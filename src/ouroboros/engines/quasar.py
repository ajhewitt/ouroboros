"""
Quasar Engine (Plan C) - Updated for DR16Q
Analyzes alignment of High-Z Quasar separation vectors with the CMB Quadrupole.
"""

import numpy as np
import healpy as hp
from astropy.table import Table
from ouroboros import config

def load_quasar_catalog(filepath, z_min=2.5):
    """
    Ingests SDSS-IV/DR16Q Quasar catalog and filters for High-Z.
    """
    print(f"Loading Quasar catalog from {filepath}...")
    # Load the FITS table
    t = Table.read(filepath, format='fits')
    
    # Standardize Column Names (Handle DR12/DR16/DESI differences)
    colnames = [c.upper() for c in t.colnames]
    
    if 'RA' in colnames:
        ra = t['RA']
        dec = t['DEC']
    elif 'TARGET_RA' in colnames:
        ra = t['TARGET_RA']
        dec = t['TARGET_DEC']
    else:
        # Fallback for some specific value-added catalogs
        print("Warning: Standard RA/DEC not found. Searching for alternates...")
        if 'RA_ICRS' in colnames:
            ra = t['RA_ICRS']
            dec = t['DEC_ICRS']
        else:
            raise ValueError(f"Could not find RA/DEC columns. Available: {t.colnames[:5]}...")
        
    if 'Z' not in colnames:
        if 'Z_PIPE' in colnames:
            z = t['Z_PIPE']
        else:
            raise ValueError("Could not find Redshift (Z) column.")
    else:
        z = t['Z']
    
    # Filter for High Redshift (The "First Light" Era)
    mask = (z > z_min)
    print(f" -> Filtering z > {z_min}: Keeping {np.sum(mask)} / {len(t)} objects.")
    
    return np.array(ra[mask]), np.array(dec[mask]), np.array(z[mask])

def get_separation_vectors(ra, dec, r_comoving):
    """
    Computes 3D separation vectors for catalog pairs.
    """
    # 1. Convert to Cartesian
    theta = np.radians(90.0 - dec)
    phi = np.radians(ra)
    
    x = r_comoving * np.sin(theta) * np.cos(phi)
    y = r_comoving * np.sin(theta) * np.sin(phi)
    z = r_comoving * np.cos(theta)
    
    coords = np.vstack((x, y, z)).T 

    # 2. Compute Differences
    # LIMITATION: For massive catalogs, we can't do N^2.
    # We restrict to N=2000 random samples if the catalog is huge.
    n_obj = len(coords)
    if n_obj > 2000:
        print(f"Warning: Downsampling {n_obj} -> 2000 for vector pair analysis.")
        rng = np.random.default_rng(42)
        indices = rng.choice(n_obj, 2000, replace=False)
        coords = coords[indices]
        n_obj = 2000
        
    diffs = coords[:, None, :] - coords[None, :, :]
    indices_i, indices_j = np.triu_indices(n_obj, k=1)
    vectors = diffs[indices_i, indices_j]
    
    norms = np.linalg.norm(vectors, axis=1)
    
    # Filter zero-length (duplicates) or NaNs
    valid = (norms > 0) & (np.isfinite(norms))
    
    return vectors[valid] / norms[valid, None]

def correlate_with_axis(vectors, axis_vector):
    """
    Measures the alignment of a set of vectors with a target axis.
    """
    dots = np.dot(vectors, axis_vector)
    return np.mean(np.abs(dots))
