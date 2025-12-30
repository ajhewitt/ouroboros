#!/usr/bin/env python3
"""
Project Ouroboros: Phase III - Tomography
"The Multipole Ladder" - Testing Kinematic Alignment across Scales.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

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
def worker_null_tomography(seed, original_map, lmin, lmax, target_vec):
    """
    Runs a single Null Simulation for a specific L-range.
    Returns alignment angle to target.
    """
    # 1. Rotate Map
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    # 2. Fast Scan (NSIDE=4 is sufficient for Nulls)
    _, scores = parity.scan_parity_directions(
        rotated_map, nside_scan=4, lmin=lmin, lmax=lmax
    )
    
    # 3. Find Max Direction
    max_idx = np.argmax(scores)
    theta, phi = hp.pix2ang(4, max_idx)
    max_vec = hp.ang2vec(theta, phi)
    
    # 4. Measure alignment to Target
    dot = np.dot(max_vec, target_vec)
    return np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))

def main(map_path=None):
    CORES = 10 
    
    # DEFINING THE SLICES
    # We slice the multipole range to see where the signal lives.
    # Bin 1: Large Scales (The "Axis of Evil" domain)
    # Bin 2: Intermediate
    # Bin 3: Smaller scales (Sound Horizon)
    BINS = [
        (2, 30, "Large Scale (L=2-30)"),
        (31, 60, "Mid Scale   (L=31-60)"),
        (61, 100, "Fine Scale  (L=61-100)")
    ]
    
    print("="*60)
    print("PROJECT OUROBOROS: PHASE III (TOMOGRAPHY)")
    print("Target: Kinematic CMB Dipole (Motion)")
    print("="*60)

    # 1. LOAD DATA
    if map_path:
        data_map = ingestion.load_map(map_path)
    else:
        print("[WARNING] No file. Using MOCK data.")
        data_map = ingestion.get_mock_map(mode='random')

    dipole_vec = config.get_dipole_vector()
    
    # 2. ITERATE THROUGH BINS
    for lmin, lmax, label in BINS:
        print(f"\n--- Analyzing {label} ---")
        
        # A. REAL SCAN
        print(f" -> Scanning Real Universe (Pool={CORES})...")
        nside_scan = 8
        npix = hp.nside2npix(nside_scan)
        
        # Update args to include lmin/lmax
        scan_args = [(data_map, i, nside_scan, lmin, lmax) for i in range(npix)]
        
        with mp.Pool(processes=CORES, maxtasksperchild=10) as pool:
            results = pool.map(parity._worker_scan_single_direction, scan_args)
        
        parity_scores = np.array(results)
        max_idx = np.argmax(parity_scores)
        theta, phi = hp.pix2ang(nside_scan, max_idx)
        max_vec = hp.ang2vec(theta, phi)
        
        real_angle = np.degrees(np.arccos(np.clip(np.dot(max_vec, dipole_vec), -1.0, 1.0)))
        print(f" -> Real Alignment to Dipole: {real_angle:.2f} deg")
        
        # B. NULL TEST
        N_SIMS = 200 # Higher N not strictly needed for scoping
        print(f" -> Running {N_SIMS} Nulls...")
        
        seeds = [2025 + i for i in range(N_SIMS)]
        worker_func = partial(
            worker_null_tomography, 
            original_map=data_map, 
            lmin=lmin, lmax=lmax, 
            target_vec=dipole_vec
        )
        
        with mp.Pool(processes=CORES, maxtasksperchild=1) as pool:
            null_angles = list(tqdm(pool.imap(worker_func, seeds), total=len(seeds)))
            
        # C. STATISTICS
        null_angles = np.array(null_angles)
        better_nulls = np.sum(null_angles < real_angle)
        p_value = better_nulls / N_SIMS
        
        print(f" -> P-Value: {p_value:.3f}")
        
        if p_value < 0.05:
            print(f" -> [RESULT] SIGNIFICANT SIGNAL in {label}")
        else:
            print(f" -> [RESULT] Null Consistent.")

if __name__ == "__main__":
    mp.set_start_method('fork') 
    path = sys.argv[1] if len(sys.argv) > 1 else None
    main(path)
