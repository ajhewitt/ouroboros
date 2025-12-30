"""
Parity Engine (Phase III Update)
Calculates the Point-Parity statistic P(n) with tomographic slicing (l_min, l_max).
"""

import numpy as np
import healpy as hp
from ouroboros import config

def get_parity_modes(alm, lmax, lmin=2):
    """
    Decomposes a set of alm coefficients into Even and Odd parity power.
    Now supports lower-bound slicing.
    """
    P_plus = 0.0
    P_minus = 0.0
    
    # Iterate through multipoles from lmin to lmax
    # We enforce lmin >= 2 to avoid Dipole/Monopole
    start_l = max(lmin, 2)
    
    for l in range(start_l, lmax + 1):
        indices = hp.Alm.getidx(lmax, l, np.arange(l + 1))
        
        cl_m0 = np.abs(alm[indices[0]])**2
        cl_m_rest = np.sum(np.abs(alm[indices[1:]])**2)
        cl_total = cl_m0 + 2 * cl_m_rest
        
        if l % 2 == 0:
            P_plus += l * (l + 1) * cl_total / (2 * np.pi)
        else:
            P_minus += l * (l + 1) * cl_total / (2 * np.pi)
            
    return P_plus, P_minus

def calculate_point_parity(map_data, lmax=config.L_MAX, lmin=config.L_MIN):
    """
    Computes P(n) for a specific multipole range.
    """
    # 1. Kinematic Isolation
    map_cleaned = hp.remove_dipole(map_data)
    
    # 2. Transform
    alm = hp.map2alm(map_cleaned, lmax=lmax)
    
    # 3. Calculate Power with slicing
    P_plus, P_minus = get_parity_modes(alm, lmax, lmin=lmin)
    
    if (P_plus + P_minus) == 0:
        return 0.0
        
    return (P_plus - P_minus) / (P_plus + P_minus)

def scan_parity_directions(map_data, nside_scan=8, lmax=config.L_MAX, lmin=config.L_MIN):
    """
    Rotates map and computes P(n) for specific range.
    """
    npix = hp.nside2npix(nside_scan)
    results = np.zeros(npix)
    theta, phi = hp.pix2ang(nside_scan, np.arange(npix))
    
    # Optimize: Pre-compute ALMs once if we were doing pure rotation algebra,
    # but for pixel-based scanning (to preserve mask behavior), we loop.
    for i in range(npix):
        r = hp.Rotator(rot=[np.degrees(phi[i]), np.degrees(theta[i])], deg=True, inv=True)
        rotated_map = r.rotate_map_pixel(map_data)
        results[i] = calculate_point_parity(rotated_map, lmax=lmax, lmin=lmin)
        
    return np.arange(npix), results

# Helper for multiprocessing
def _worker_scan_single_direction(args):
    """
    Unpacks arguments including lmin/lmax.
    Args: (map_data, pixel_index, nside_scan, lmin, lmax)
    """
    map_data, pix_idx, nside_scan, lmin, lmax = args
    theta, phi = hp.pix2ang(nside_scan, pix_idx)
    r = hp.Rotator(rot=[np.degrees(phi), np.degrees(theta)], deg=True, inv=True)
    rotated_map = r.rotate_map_pixel(map_data)
    return calculate_point_parity(rotated_map, lmax=lmax, lmin=lmin)
