"""
Geometry Engine (Plan B) - Updated
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
    """
    # 1. Smooth the map
    sigma_rad = np.radians(smooth_fwhm_deg) / 2.355
    # Remove verbose=False for newer healpy versions
    smoothed_map = hp.smoothing(map_data, sigma=sigma_rad)
    
    # 2. Mask the Galactic Plane and Northern Hemisphere
    nside = hp.npix2nside(len(map_data))
    npix = len(map_data)
    
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    b_rad = np.pi/2 - theta 
    
    # Mask everything above b = -20 degrees
    mask_indices = np.where(np.degrees(b_rad) > -20)
    
    search_map = smoothed_map.copy()
    search_map[mask_indices] = np.inf 
    
    # 3. Find Minimum
    min_idx = np.argmin(search_map)
    
    # 4. Convert to RA/Dec (ICRS)
    gl_rad, gb_rad = hp.pix2ang(nside, min_idx, lonlat=True)
    coord = SkyCoord(l=gl_rad*u.deg, b=gb_rad*u.deg, frame='galactic')
    icrs = coord.transform_to('icrs')
    
    return icrs.ra.degree, icrs.dec.degree

def check_nodal_alignment(ra_spot, dec_spot):
    """
    Calculates distances from the Cold Spot to ALL Agency Nodes.
    """
    spot_vec = hp.ang2vec(ra_spot, dec_spot, lonlat=True)
    
    # 1. Ecliptic Poles (Orbit)
    nep_vec = hp.ang2vec(config.ECLIPTIC_POLE_RA_DEG, config.ECLIPTIC_POLE_DEC_DEG, lonlat=True)
    sep_vec = -1 * nep_vec 
    
    # 2. Equinox (Intersection of Orbit/Equator)
    equinox_vec = hp.ang2vec(0.0, 0.0, lonlat=True)
    
    # 3. Solar Spin (Helicity) - Using raw constants to ensure accuracy
    solar_vec = hp.ang2vec(config.SOLAR_NORTH_POLE_RA_DEG, config.SOLAR_NORTH_POLE_DEC_DEG, lonlat=True)
    
    # 4. CMB Dipole (Motion)
    dipole_vec = hp.ang2vec(config.CMB_DIPOLE_RA_DEG, config.CMB_DIPOLE_DEC_DEG, lonlat=True)
    
    def get_angle(v1, v2):
        dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
        return np.degrees(np.arccos(dot))
        
    return {
        "dist_nep": get_angle(spot_vec, nep_vec),
        "dist_sep": get_angle(spot_vec, sep_vec),
        "dist_equinox": get_angle(spot_vec, equinox_vec),
        "dist_solar_spin": get_angle(spot_vec, solar_vec),
        "dist_dipole": get_angle(spot_vec, dipole_vec)
    }
