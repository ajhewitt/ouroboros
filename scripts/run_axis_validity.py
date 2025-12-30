#!/usr/bin/env python3
"""
Test: Axis of Evil Validity
Checks if the intrinsic alignment (l=2 vs l=3) is statistically anomalous.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import numpy as np
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from ouroboros import ingestion, config
from ouroboros.engines import harmonics
from ouroboros.validation import nulling

def worker_internal_alignment(seed, original_map):
    # 1. Generate a Null Universe (Random Phases)
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    # 2. Extract Axes
    # harmonics.analyze_axis_of_evil returns the angle between l=2 and l=3
    results = harmonics.analyze_axis_of_evil(rotated_map)
    
    return results['angle_23']

def main(path):
    print("="*60)
    print("AXIS OF EVIL VALIDITY CHECK")
    print("Hypothesis: l=2 and l=3 are coupled (Parallel Planes).")
    print("="*60)

    # 1. LOAD DATA
    data_map = ingestion.load_map(path)
    
    # 2. REAL SCORE
    real_results = harmonics.analyze_axis_of_evil(data_map)
    real_angle = real_results['angle_23']
    print(f"\n[1] Real Sky Separation (l=2 vs l=3): {real_angle:.2f} deg")

    # 3. NULL SIMULATION
    N_SIMS = 2000 # Let's do 2000 to be sure (target p < 0.003 for 3-sigma)
    print(f"[2] Simulating {N_SIMS} Random Universes...")
    
    worker = partial(worker_internal_alignment, original_map=data_map)
    seeds = [3000 + i for i in range(N_SIMS)]
    
    with mp.Pool(10) as pool:
        null_angles = list(tqdm(pool.imap(worker, seeds), total=len(seeds)))
        
    # 4. FORENSICS
    null_angles = np.array(null_angles)
    
    # How many random universes are MORE aligned (smaller angle) than ours?
    better_nulls = np.sum(null_angles < real_angle)
    p_value = better_nulls / N_SIMS
    
    print("\n[3] FINAL VERDICT")
    print(f" -> Random Universes with tighter alignment: {better_nulls}/{N_SIMS}")
    print(f" -> P-Value: {p_value:.4f}")
    
    if p_value < 0.01:
        print("\n[RESULT] THE AXIS IS REAL.")
        print("The l=2 and l=3 modes are physically coupled. This requires New Physics.")
    elif p_value < 0.05:
        print("\n[RESULT] TANTALIZING HINT.")
        print("Statistically significant (2-sigma), but not definitive.")
    else:
        print("\n[RESULT] THE AXIS IS NOISE.")
        print("A 9-degree separation happens frequently by chance. Mystery Solved.")

if __name__ == "__main__":
    mp.set_start_method('fork')
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
