"""
Data Ingestion Layer
Handles loading of Planck FITS files and generation of synthetic mock data.
"""

import numpy as np
import healpy as hp
from ouroboros import config

def load_map(filepath, field=0):
    """
    Loads a CMB map from a FITS file and standardizes resolution.
    
    Args:
        filepath (str): Path to the .fits file.
        field (int): The FITS column to read (0=I, 1=Q, 2=U). Default 0 (Temperature).
        
    Returns:
        array: HEALPix map standardized to config.NSIDE.
    """
    print(f"Loading map from {filepath}...")
    # Load map (read-only to save memory)
    raw_map = hp.read_map(filepath, field=field, verbose=False)
    
    # Get current resolution
    nside_in = hp.npix2nside(len(raw_map))
    
    # Downgrade if necessary (Pipeline runs at NSIDE=64 for speed/LSS focus)
    if nside_in > config.NSIDE:
        print(f"Downgrading map from NSIDE {nside_in} to {config.NSIDE}...")
        return hp.ud_grade(raw_map, config.NSIDE)
    elif nside_in < config.NSIDE:
        raise ValueError(f"Input map resolution ({nside_in}) is lower than pipeline requirement ({config.NSIDE}).")
        
    return raw_map

def get_mock_map(mode='random'):
    """
    Generates a synthetic CMB map for pipeline testing.
    
    Args:
        mode (str): 'random' for noise, 'signal' to inject a fake Solar alignment.
    """
    npix = hp.nside2npix(config.NSIDE)
    
    if mode == 'random':
        # Standard Gaussian random field
        return np.random.randn(npix)
    
    elif mode == 'signal':
        # Inject a strong Odd-Parity signal aligned with the Solar Pole
        # 1. Create Odd-Parity Alm (l=3, m=0)
        alm_size = hp.Alm.getsize(config.L_MAX)
        alm = np.zeros(alm_size, dtype=np.complex128)
        idx = hp.Alm.getidx(config.L_MAX, 3, 0)
        alm[idx] = 50.0 # Strong signal
        
        # 2. Create Map
        base_map = hp.alm2map(alm, config.NSIDE)
        
        # 3. Rotate to align with Solar Pole
        # We rotate the MAP, which moves the features.
        solar_vec = config.get_solar_vector()
        # Convert vector to angles
        theta, phi = hp.vec2ang(solar_vec)
        
        # Rotator: (phi, theta, psi)
        r = hp.Rotator(rot=[np.degrees(phi), np.degrees(theta)], deg=True, inv=True)
        return r.rotate_map_pixel(base_map)
