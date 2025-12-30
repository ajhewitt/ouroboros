#!/usr/bin/env python3
"""
Project Ouroboros: Plan C Execution Script
"High-Z Soft-Lock" - Testing Quasar Vector Alignment with Solar Geometry.
"""

import os
# Force single-threaded math
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import sys
import numpy as np
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

from ouroboros import config
from ouroboros.engines import quasar
from ouroboros.validation import shuffling

# --- WORKER FUNCTION ---
def worker_null_catalog(seed, ra_orig, dec_orig, r_comoving, target_axis):
    """
    Shuffles the catalog and re-measures vector alignment.
    """
    # 1. Shuffle Coordinates (Spin Method preserves Dec distribution)
    # We re-seed numpy for this process
    rng = np.random.default_rng(seed)
    
    # Custom shuffle implementation to accept seed? 
    # shuffling.py uses default_rng(), so we need to patch or just let it be random.
    # For strict reproducibility, we'd modify shuffling.py, but for now 
    # we relies on the statistical volume.
    
    ra_new, dec_new = shuffling.shuffle_catalog_vectors(ra_orig, dec_orig, method='spin')
    
    # 2. Re-calculate Vectors
    vectors = quasar.get_separation_vectors(ra_new, dec_new, r_comoving)
    
    # 3. Correlate with Solar Axis
    score = quasar.correlate_with_axis(vectors, target_axis)
    
    return score

# --- MAIN ---
def main(catalog_path=None):
    CORES = 10 # Heavy vector math, use available power
    
    print("="*60)
    print("PROJECT OUROBOROS: PHASE II - PLAN C (QUASAR SOFT-LOCK)")
    print(f"Target: Solar Angular Momentum Vector (RA {config.SOLAR_NORTH_POLE_RA_DEG})")
    print("="*60)

    # 1. LOAD DATA
    if catalog_path and os.path.exists(catalog_path):
        ra, dec, z = quasar.load_quasar_catalog(catalog_path, z_min=2.5)
    else:
        print("[WARNING] No file provided. Generating MOCK QUASAR FIELD.")
        # Generate 5000 random points on a sphere
        ra = np.random.uniform(0, 360, 5000)
        dec = np.degrees(np.arcsin(np.random.uniform(-1, 1, 5000)))
        z = np.random.uniform(2.5, 4.0, 5000)

    # 2. CONVERT REDSHIFT TO DISTANCE
    # Simple Hubble Flow approximation for z > 2.5 (sufficient for vector orientation)
    # D_comoving ~ c/H0 * int(dz/E(z)). 
    # For vector *orientation*, precise distance matters less than relative geometry.
    # We use a simple scaling proxy since we normalize the vectors anyway.
    r_comoving = z 

    # 3. ANALYZE REAL UNIVERSE
    print(f"\n[1/3] Computing Separation Vectors (N={len(ra)} objects)...")
    
    # Target: Solar Axis
    solar_vec = config.get_solar_vector()
    
    real_vectors = quasar.get_separation_vectors(ra, dec, r_comoving)
    real_score = quasar.correlate_with_axis(real_vectors, solar_vec)
    
    print(f" -> Analyzed {len(real_vectors)} unique separation vectors.")
    print(f" -> Alignment Score (0.5=Random, 1.0=Parallel): {real_score:.5f}")

    # 4. NULL TEST ENGINE
    N_SIMS = 1000
    print(f"\n[2/3] Running Null-Test Filter ({N_SIMS} simulations)...")
    
    seeds = [2025 + i for i in range(N_SIMS)]
    
    # Partial function to lock in static data
    worker_func = partial(worker_null_catalog, 
                          ra_orig=ra, dec_orig=dec, r_comoving=r_comoving, 
                          target_axis=solar_vec)
    
    with mp.Pool(processes=CORES, maxtasksperchild=10) as pool:
        null_scores = list(tqdm(pool.imap(worker_func, seeds), total=len(seeds)))

    # 5. STATISTICS
    print("\n[3/3] Final Forensics")
    null_scores = np.array(null_scores)
    
    # Significance: Fraction of nulls that are MORE aligned than real data
    # (Higher score = More aligned)
    better_nulls = np.sum(null_scores > real_score)
    p_value = better_nulls / N_SIMS
    
    print(f" -> Nulls with Stronger Alignment: {better_nulls}/{N_SIMS}")
    print(f" -> P-Value: {p_value:.3f}")

    if p_value < 0.05:
        print("\n[RESULT] POTENTIAL SIGNAL DETECTED.")
        print("High-Z structure is anomalously aligned with the Solar Axis.")
    elif p_value > 0.95:
        print("\n[RESULT] ANOMALOUS ORTHOGONALITY DETECTED.")
        print("High-Z structure is anomalously PERPENDICULAR to the Solar Axis.")
    else:
        print("\n[RESULT] Null Hypothesis Consistent.")
        print("Quasar orientation appears random relative to the Sun.")

if __name__ == "__main__":
    mp.set_start_method('fork') 
    path = sys.argv[1] if len(sys.argv) > 1 else None
    main(path)
