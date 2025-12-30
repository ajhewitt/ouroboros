"""
Geometry Engine (Plan B)
Analyzes the positional alignment of the Cold Spot relative to Solar Geometry.
"""

import numpy as np
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u
from ouroboros import config

def find_cold_spot(map_data, smooth_fwhm_deg=5.0):
    """
    Locates the center of the Cold Spot (Eridanus Supervoid).
    
    Args:
        map_data (array): Temperature map (Kelvin).
        smooth_fwhm_deg (float): Smoothing kernel size (default 5 deg for large voids).
        
    Returns:
        tuple: (RA, Dec) of the coldest point in degrees.
    """
    # 1. Smooth the map to suppress point sources/noise
    # We use sigma = FWHM / 2.355
    sigma_rad = np.radians(smooth_fwhm_deg) / 2.355
    smoothed_map = hp.smoothing(map_data, sigma=sigma_rad)
    
    # 2. Mask the Galactic Plane and Northern Hemisphere
    # The Cold Spot is known to be in the Southern Galactic Hemisphere.
    nside = hp.npix2nside(len(map_data))
    npix = len(map_data)
    
    # Get coordinates for all pixels
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    b_rad = np.pi/2 - theta # Galactic Latitude (if input is Galactic)
    
    # We assume input map is in GALACTIC coordinates (Standard for Planck)
    # Mask everything above b = -20 degrees to focus on Eridanus region
    mask_indices = np.where(np.degrees(b_rad) > -20)
    
    # Create a working copy
    search_map = smoothed_map.copy()
    search_map[mask_indices] = np.inf # Exclude from minimum search
    
    # 3. Find Minimum
    min_idx = np.argmin(search_map)
    
    # 4. Convert to RA/Dec (ICRS)
    # Input is Galactic, so we must convert out
    gl_rad, gb_rad = hp.pix2ang(nside, min_idx, lonlat=True)
    
    coord = SkyCoord(l=gl_rad*u.deg, b=gb_rad*u.deg, frame='galactic')
    icrs = coord.transform_to('icrs')
    
    return icrs.ra.degree, icrs.dec.degree

def check_nodal_alignment(ra_spot, dec_spot):
    """
    Calculates distances from the Cold Spot to key Solar System nodes.
    
    Nodes checked:
    1. North Ecliptic Pole (NEP)
    2. South Ecliptic Pole (SEP)
    3. Vernal Equinox (Intersection of Ecliptic & Equator)
    """
    # Define Nodes
    spot_vec = hp.ang2vec(ra_spot, dec_spot, lonlat=True)
    
    # Ecliptic Poles
    nep_vec = hp.ang2vec(config.ECLIPTIC_POLE_RA_DEG, config.ECLIPTIC_POLE_DEC_DEG, lonlat=True)
    # South Pole is opposite
    sep_vec = -1 * nep_vec 
    
    # Vernal Equinox (defined as RA=0, Dec=0 in J2000/ICRS)
    # This is the "Zero Point" of the celestial sphere, tied to Earth's orbit/spin intersection.
    equinox_vec = hp.ang2vec(0.0, 0.0, lonlat=True)
    
    # Calculate Distances (Geodesic angle)
    def get_angle(v1, v2):
        dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
        return np.degrees(np.arccos(dot))
        
    return {
        "dist_nep": get_angle(spot_vec, nep_vec),
        "dist_sep": get_angle(spot_vec, sep_vec),
        "dist_equinox": get_angle(spot_vec, equinox_vec)
    }
