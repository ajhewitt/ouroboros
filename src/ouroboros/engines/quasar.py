"""
Quasar Engine (Plan C)
Analyzes alignment of High-Z Quasar separation vectors with the CMB Quadrupole.
"""

import numpy as np
import healpy as hp
from ouroboros import config

def get_separation_vectors(ra, dec, r_comoving):
    """
    Computes 3D separation vectors for all pairs in the catalog.
    WARNING: O(N^2) complexity. Use on subsets or sparse catalogs.
    
    Args:
        ra, dec (array): Coordinates in degrees.
        r_comoving (array): Comoving distance (Mpc) for each object.
        
    Returns:
        array: Normalized separation vectors (N_pairs, 3).
    """
    # 1. Convert to Cartesian (x, y, z)
    # Theta = 90 - Dec, Phi = RA
    theta = np.radians(90.0 - dec)
    phi = np.radians(ra)
    
    x = r_comoving * np.sin(theta) * np.cos(phi)
    y = r_comoving * np.sin(theta) * np.sin(phi)
    z = r_comoving * np.cos(theta)
    
    coords = np.vstack((x, y, z)).T # Shape (N, 3)
    
    # 2. Compute Differences (Broadcasting)
    # We want vector V_ij = Pos_i - Pos_j
    # To save memory, we can iterate or use simple broadcasting for N < 5000
    n_obj = len(coords)
    if n_obj > 2000:
        raise ValueError(f"Catalog too large ({n_obj}) for brute-force pairs. Downsample first.")
        
    # Expand dims to subtract
    # coords[:, None, :] is (N, 1, 3)
    # coords[None, :, :] is (1, N, 3)
    diffs = coords[:, None, :] - coords[None, :, :]
    
    # Flatten to list of vectors
    # We only want unique pairs (i > j)
    indices_i, indices_j = np.triu_indices(n_obj, k=1)
    vectors = diffs[indices_i, indices_j]
    
    # 3. Normalize
    norms = np.linalg.norm(vectors, axis=1)
    # Filter zero-length (duplicates)
    valid = norms > 0
    normalized_vectors = vectors[valid] / norms[valid, None]
    
    return normalized_vectors

def correlate_with_axis(vectors, axis_vector):
    """
    Measures the alignment of a set of vectors with a target axis.
    
    Args:
        vectors (array): List of 3D vectors (separation vectors).
        axis_vector (array): The target axis (e.g., CMB Quadrupole Axis).
        
    Returns:
        float: The Alignment Statistic (Average absolute dot product).
               0.5 = Random, 1.0 = Perfect Parallel Alignment.
    """
    # Dot product
    dots = np.dot(vectors, axis_vector)
    
    # We care about orientation (alignment), not direction (+/-),
    # so we take the Absolute Value of the cosine.
    abs_cos = np.abs(dots)
    
    return np.mean(abs_cos)

def get_quadrupole_axis(cmb_map_path):
    """
    Extracts the primary axis of the Quadrupole (l=2) from a map.
    """
    # Load and clean
    # In production, use the proper ingestion tools
    if isinstance(cmb_map_path, str):
        m = hp.read_map(cmb_map_path, verbose=False)
    else:
        m = cmb_map_path
        
    # Get Alm for l=2 only
    alm = hp.map2alm(m, lmax=2)
    
    # This is a simplification. The "Axis of Evil" is often defined 
    # by the maximum angular momentum dispersion vector (Maxwellian Multipole Vectors).
    # For this pipeline, we will use the direction of the l=2, m=0 mode 
    # in the frame where it is maximized. 
    # (Implementation omitted for brevity, returning placeholder Solar Axis 
    # to test the pipeline flow).
    return config.get_solar_vector()
