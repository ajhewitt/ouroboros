"""
Parity Engine (Plan A)
Calculates the Point-Parity statistic P(n) for CMB maps.
"""

import numpy as np
import healpy as hp
from ouroboros import config

def get_parity_modes(alm, lmax):
    """
    Decomposes a set of alm coefficients into Even and Odd parity power.
    
    Args:
        alm (array): Spherical harmonic coefficients (complex).
        lmax (int): Maximum multipole to consider.
        
    Returns:
        tuple: (P_plus, P_minus) represents the total power in Even/Odd modes.
    """
    # Initialize power accumulators
    P_plus = 0.0
    P_minus = 0.0
    
    # Iterate through multipoles from l_min to lmax
    for l in range(config.L_MIN, lmax + 1):
        # Extract indices for this l
        # hp.Alm.getidx generates indices for all m (0 to l) at this l
        indices = hp.Alm.getidx(lmax, l, np.arange(l + 1))
        
        # Calculate C_l (Power Spectrum) for this specific multipole
        # Sum |alm|^2 over all m. Factor of 2 handles m < 0 symmetries.
        # Note: m=0 term is counted once, m>0 terms are doubled.
        cl_m0 = np.abs(alm[indices[0]])**2
        cl_m_rest = np.sum(np.abs(alm[indices[1:]])**2)
        cl_total = cl_m0 + 2 * cl_m_rest
        
        # Determine Parity
        # Even Parity: l(l+1) is Even -> l is Even (2, 4, 6...)
        # Odd Parity: l is Odd (3, 5, 7...)
        if l % 2 == 0:
            P_plus += l * (l + 1) * cl_total / (2 * np.pi)
        else:
            P_minus += l * (l + 1) * cl_total / (2 * np.pi)
            
    return P_plus, P_minus

def calculate_point_parity(map_data, lmax=config.L_MAX):
    """
    Computes the Point-Parity statistic P(n) for a given map.
    """
    # 1. Kinematic Isolation Check
    # FIXED: Removed deprecated 'verbose' argument
    map_cleaned = hp.remove_dipole(map_data)
    
    # 2. Spherical Harmonic Transform
    alm = hp.map2alm(map_cleaned, lmax=lmax)
    
    # 3. Calculate Power in Parity Modes
    P_plus, P_minus = get_parity_modes(alm, lmax)
    
    # 4. Compute Statistic
    if (P_plus + P_minus) == 0:
        return 0.0
        
    return (P_plus - P_minus) / (P_plus + P_minus)

def scan_parity_directions(map_data, nside_scan=8):
    """
    Rotates the map to align with every pixel center in a low-res HEALPix grid
    and computes the P(n) statistic for that orientation.
    
    Args:
        map_data (array): The input CMB map.
        nside_scan (int): Resolution of the scan grid (default 8 -> 768 pixels).
        
    Returns:
        tuple: (pixel_indices, parity_values)
    """
    npix = hp.nside2npix(nside_scan)
    results = np.zeros(npix)
    
    # Pre-calculate rotation angles for all pixels
    theta, phi = hp.pix2ang(nside_scan, np.arange(npix))
    
    # Scan loop (conceptually moving the observer's "Up" vector)
    # Note: Rotating the map is inverse to rotating the frame.
    # To test direction (theta, phi), we rotate the map so that (theta, phi)
    # becomes the new North Pole.
    
    # Create a rotator for each direction
    # This is computationally expensive; in production we might parallelize this.
    for i in range(npix):
        # Define rotation: (phi, theta, psi) 
        # We rotate by -phi around Z, then -theta around Y to bring point to Pole.
        r = hp.Rotator(rot=[np.degrees(phi[i]), np.degrees(theta[i])], deg=True, inv=True)
        
        # Rotate the map
        rotated_map = r.rotate_map_pixel(map_data)
        
        # Compute Parity
        results[i] = calculate_point_parity(rotated_map)
        
    return np.arange(npix), results

def _worker_scan_single_direction(args):
    """
    Helper function for multiprocessing.
    Args must be packed as a tuple: (map_data, pixel_index, nside_scan)
    """
    map_data, pix_idx, nside_scan = args
    
    # Reconstruct direction from index
    theta, phi = hp.pix2ang(nside_scan, pix_idx)
    
    # Create Rotator (Inverse rotation to bring pixel to North Pole)
    r = hp.Rotator(rot=[np.degrees(phi), np.degrees(theta)], deg=True, inv=True)
    
    # Rotate Map
    rotated_map = r.rotate_map_pixel(map_data)
    
    # Calculate Score
    return calculate_point_parity(rotated_map)
