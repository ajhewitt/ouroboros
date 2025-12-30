#!/usr/bin/env python3
"""
Phase IV-B: Null Test for Harmonic Alignment
Checks if the l=2/3 alignment with the Ecliptic Pole is a statistical fluke.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import sys
import numpy as np
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
import healpy as hp

from ouroboros import config, ingestion
from ouroboros.engines import harmonics
from ouroboros.validation import nulling

def worker_harmonic_null(seed, original_map):
    # 1. Rotate Map
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    # 2. Extract Axes
    # We suppress output to keep the progress bar clean
    # Note: harmonics.analyze_axis_of_evil performs the extraction
    results = harmonics.analyze_axis_of_evil(rotated_map)
    
    # We care about the alignment to the Ecliptic Pole (fixed in lab frame)
    # Since we rotated the MAP, the extracted axis moves.
    # The Ecliptic Pole stays at the reference coordinates.
    
    return results['angle_3_ecl']

def main(map_path):
    CORES = 10
    N_SIMS = 1000
    
    print("="*60)
    print("PHASE IV-B: AXIS OF EVIL NULL TEST")
    print("Hypothesis: The Octopole (l=3) is locked to the Ecliptic Pole.")
    print("="*60)

    # 1. LOAD DATA
    data_map = ingestion.load_map(map_path)
    
    # 2. REAL SCORE
    print("\n[1] Measuring Real Alignment...")
    real_results = harmonics.analyze_axis_of_evil(data_map)
    real_angle = real_results['angle_3_ecl']
    print(f" -> Real Angle to Ecliptic Pole: {real_angle:.2f} deg")

    # 3. NULL TEST
    print(f"\n[2] Running {N_SIMS} Null Simulations...")
    worker_func = partial(worker_harmonic_null, original_map=data_map)
    seeds = [2025 + i for i in range(N_SIMS)]
    
    with mp.Pool(processes=CORES) as pool:
        null_angles = list(tqdm(pool.imap(worker_func, seeds), total=len(seeds)))
        
    # 4. STATS
    null_angles = np.array(null_angles)
    
    # We are testing for "Close Alignment" (Small angles)
    # How many nulls are CLOSER to the pole than the real data?
    better_nulls = np.sum(null_angles < real_angle)
    p_value = better_nulls / N_SIMS
    
    print("\n[3] FINAL FORENSICS")
    print(f" -> Nulls with closer alignment: {better_nulls}/{N_SIMS}")
    print(f" -> P-Value: {p_value:.3f}")
    
    if p_value < 0.05:
        print("\n[RESULT] SIGNIFICANT SIGNAL DETECTED.")
        print("The alignment of the Octopole with the Ecliptic is non-random.")
    else:
        print("\n[RESULT] Null Hypothesis Consistent.")

if __name__ == "__main__":
    mp.set_start_method('fork') 
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
