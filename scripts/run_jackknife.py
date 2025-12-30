#!/usr/bin/env python3
"""
Phase IV-D: The Galactic Jackknife
Re-runs the Ecliptic Variance test, but strictly ignores the Galactic Plane (|b| < 40).
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import numpy as np
import healpy as hp
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

# CORRECTED IMPORTS
from ouroboros import ingestion
from ouroboros.validation import nulling 

# Global for efficiency
ECLIPTIC_INDICES = None
GALACTIC_LATS = None

def worker_jackknife(seed, original_map):
    # 1. Generate Null Map
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    # 2. Extract Ecliptic Strip Pixels
    strip_pixels = rotated_map[ECLIPTIC_INDICES]
    
    # 3. Apply Pre-computed Galactic Cut (|b| > 40)
    # GALACTIC_LATS corresponds index-to-index with ECLIPTIC_INDICES
    clean_mask = np.abs(GALACTIC_LATS) > 40
    clean_pixels = strip_pixels[clean_mask]
    
    # 4. Remove Masked/Bad Data
    valid = clean_pixels[clean_pixels > -1e20]
    
    if len(valid) == 0: 
        return 0.0
        
    return np.var(valid)

def main(path):
    global ECLIPTIC_INDICES, GALACTIC_LATS
    
    print("="*60)
    print("PHASE IV-D: GALACTIC JACKKNIFE")
    print("Checking if the Variance Anomaly survives a |b| > 40 cut.")
    print("="*60)
    
    # 1. Load Data
    data_map = ingestion.load_map(path)
    nside = hp.npix2nside(len(data_map))
    
    # 2. Precompute Ecliptic Indices (The Strip)
    print(" -> Precomputing Ecliptic/Galactic geometry...")
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    
    r = hp.Rotator(coord=['G', 'E'])
    theta_ecl, _ = r(theta, phi)
    b_ecl = 90.0 - np.degrees(theta_ecl)
    
    # Define Ecliptic Belt (+/- 20 deg)
    ECLIPTIC_INDICES = np.where(np.abs(b_ecl) < 20)[0]
    
    # 3. Precompute Galactic Lats for those specific pixels
    # theta is Galactic Co-Latitude
    theta_gal = theta[ECLIPTIC_INDICES]
    GALACTIC_LATS = 90.0 - np.degrees(theta_gal)
    
    print(f" -> Ecliptic Belt Pixels: {len(ECLIPTIC_INDICES)}")
    print(f" -> After Galactic Cut (|b|>40): {np.sum(np.abs(GALACTIC_LATS) > 40)}")
    
    # 4. Real Score
    real_var = worker_jackknife(0, data_map) # Seed 0 preserves map
    print(f" -> Real Variance (High-Latitude Only): {real_var:.2f} uK^2")
    
    # 5. Null Test
    N_SIMS = 1000
    seeds = [2025 + i for i in range(N_SIMS)]
    worker = partial(worker_jackknife, original_map=data_map)
    
    print(f"\n[2] Running {N_SIMS} Jackknife Nulls...")
    with mp.Pool(10) as pool:
        nulls = list(tqdm(pool.imap(worker, seeds), total=len(seeds)))
        
    # 6. Stats
    nulls = np.array(nulls)
    louder_nulls = np.sum(nulls > real_var)
    p_value = louder_nulls / N_SIMS
    
    print(f"\n[RESULT] Nulls with Higher Variance: {louder_nulls}/{N_SIMS}")
    print(f" -> P-Value: {p_value:.3f}")
    
    if p_value < 0.05:
        print("\n[CONCLUSION] SIGNAL PERSISTS.")
        print("The Ecliptic is anomalously noisy even without the Galaxy.")
    elif p_value > 0.95:
        print("\n[CONCLUSION] SIGNAL INVERTED.")
        print("The Ecliptic is anomalously QUIET (Shielding).")
    else:
        print("\n[CONCLUSION] Signal Vanished (Likely Galactic Contamination).")

if __name__ == "__main__":
    mp.set_start_method('fork')
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
