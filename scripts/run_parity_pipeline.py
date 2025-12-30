#!/usr/bin/env python3
"""
Project Ouroboros: Plan A Execution Script (Memory-Safe)
"""

import os
# Force internal math libraries to be single-threaded.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import numpy as np
import healpy as hp
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

from ouroboros import config, ingestion
from ouroboros.engines import parity
from ouroboros.validation import nulling

# --- WORKER FUNCTIONS ---
def worker_null_simulation(seed, original_map, solar_vec):
    """
    Runs a single Null Simulation in a separate process.
    """
    # 1. Generate Rotated Map
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    # 2. Run the Scan on this rotated universe (Low Res for Speed)
    # Using nside_scan=4 (192 directions) instead of 8 (3072) for the nulls
    # saves massive time/memory while preserving statistical validity.
    _, scores = parity.scan_parity_directions(rotated_map, nside_scan=4)
    
    # 3. Find Max Direction
    max_idx = np.argmax(scores)
    theta, phi = hp.pix2ang(4, max_idx)
    max_vec = hp.ang2vec(theta, phi)
    
    # 4. Measure alignment
    return get_alignment_angle(max_vec, solar_vec)

def get_alignment_angle(vec_a, vec_b):
    dot = np.dot(vec_a, vec_b)
    dot = np.clip(dot, -1.0, 1.0)
    return np.degrees(np.arccos(dot))

# --- MAIN EXECUTION ---
def main(map_path=None):
    # USE 11 CORES (Leave 1 for OS)
    CORES = 11 
    print("="*60)
    print(f"PROJECT OUROBOROS: PARALLEL MODE ({CORES} Workers)")
    print("="*60)

    # 1. LOAD DATA
    if map_path:
        data_map = ingestion.load_map(map_path)
    else:
        print("[WARNING] No file. Generating MOCK SIGNAL...")
        data_map = ingestion.get_mock_map(mode='signal')

    # 2. PARALLEL PARITY SCAN (Real Universe)
    print(f"\n[1/3] Scanning 3,072 Directions (Pool={CORES})...")
    
    # We use a chunksize to batch tasks, which reduces overhead
    nside_scan = 8
    npix = hp.nside2npix(nside_scan)
    scan_args = [(data_map, i, nside_scan) for i in range(npix)]
    
    # maxtasksperchild=10 ensures workers clean up periodically during the long scan
    with mp.Pool(processes=CORES, maxtasksperchild=10) as pool:
        results = pool.map(parity._worker_scan_single_direction, scan_args)
        
    parity_scores = np.array(results)
    
    # Analyze Real Results
    max_idx = np.argmax(parity_scores)
    max_val = parity_scores[max_idx]
    theta, phi = hp.pix2ang(nside_scan, max_idx)
    max_vec = hp.ang2vec(theta, phi)
    solar_vec = config.get_solar_vector()
    alignment_angle = get_alignment_angle(max_vec, solar_vec)
    
    print(f" -> Max Score: {max_val:.4f}")
    print(f" -> Distance from Solar Pole: {alignment_angle:.2f} deg")

    # 3. PARALLEL NULL TEST ENGINE (Memory Critical Section)
    print(f"\n[2/3] Running {config.N_SIMS} Null Simulations (Pool={CORES})...")
    
    seeds = [2025 + i for i in range(config.N_SIMS)]
    worker_func = partial(worker_null_simulation, original_map=data_map, solar_vec=solar_vec)
    
    # CRITICAL FIX: maxtasksperchild=1
    # This forces the worker process to DIE and Restart after every single simulation.
    # This guarantees that any memory leaked by healpy/numpy is returned to the OS immediately.
    with mp.Pool(processes=CORES, maxtasksperchild=1) as pool:
        null_angles = list(tqdm(pool.imap(worker_func, seeds), total=len(seeds)))

    # 4. STATISTICS
    print("\n[3/3] Final Forensics")
    null_angles = np.array(null_angles)
    better_nulls = np.sum(null_angles < alignment_angle)
    p_value = better_nulls / config.N_SIMS
    
    print(f" -> Nulls Closer to Sun: {better_nulls}/{config.N_SIMS}")
    print(f" -> P-Value: {p_value:.3f}")
    
    if p_value < 0.05:
        print("\n[RESULT] POTENTIAL SIGNAL DETECTED.")
    else:
        print("\n[RESULT] Null Hypothesis Consistent.")

if __name__ == "__main__":
    mp.set_start_method('fork') 
    path = sys.argv[1] if len(sys.argv) > 1 else None
    main(path)
