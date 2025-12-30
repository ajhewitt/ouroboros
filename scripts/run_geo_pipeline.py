#!/usr/bin/env python3
"""
Project Ouroboros: Plan B Execution Script (Omni-Check)
Checks Cold Spot alignment against Ecliptic, Solar Spin, and Dipole.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import sys
import numpy as np
import multiprocessing as mp
from functools import partial
from tqdm import tqdm

from ouroboros import config, ingestion
from ouroboros.engines import geo
from ouroboros.validation import nulling

def worker_null_spot(seed, original_map):
    """Returns the distance specifically to the Solar Spin Axis for the null."""
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    ra, dec = geo.find_cold_spot(rotated_map, smooth_fwhm_deg=5.0)
    metrics = geo.check_nodal_alignment(ra, dec)
    
    # We return the Solar Spin distance specifically for the P-value calculation
    return metrics['dist_solar_spin']

def main(map_path=None):
    CORES = 6
    print("="*60)
    print("PROJECT OUROBOROS: PLAN B (OMNI-CHECK)")
    print("Targeting: Ecliptic, Solar Spin, and Dipole")
    print("="*60)

    # 1. LOAD DATA
    if map_path:
        data_map = ingestion.load_map(map_path)
    else:
        print("[WARNING] No file. Using MOCK data.")
        data_map = ingestion.get_mock_map()

    # 2. ANALYZE REAL UNIVERSE
    print(f"\n[1/3] Locating Eridanus Supervoid...")
    real_ra, real_dec = geo.find_cold_spot(data_map, smooth_fwhm_deg=5.0)
    real_metrics = geo.check_nodal_alignment(real_ra, real_dec)
    
    print(f" -> Found Cold Spot at: RA {real_ra:.2f}, Dec {real_dec:.2f}")
    print("-" * 40)
    print(f" -> Dist to N. Ecliptic Pole: {real_metrics['dist_nep']:.2f} deg")
    print(f" -> Dist to S. Ecliptic Pole: {real_metrics['dist_sep']:.2f} deg")
    print(f" -> Dist to Vernal Equinox:   {real_metrics['dist_equinox']:.2f} deg")
    print("-" * 40)
    print(f" -> Dist to CMB Dipole:       {real_metrics['dist_dipole']:.2f} deg")
    print(f" -> Dist to Solar Spin Axis:  {real_metrics['dist_solar_spin']:.2f} deg")
    
    # We focus significance testing on the Solar Spin, as Ecliptic was already tested (p=0.614)
    target_metric = real_metrics['dist_solar_spin']

    # 3. NULL TEST ENGINE
    N_SIMS = 1000
    print(f"\n[2/3] Running Null-Test (Focus: Solar Spin, N={N_SIMS})...")
    
    seeds = [2025 + i for i in range(N_SIMS)]
    worker_func = partial(worker_null_spot, original_map=data_map)
    
    with mp.Pool(processes=CORES, maxtasksperchild=10) as pool:
        null_dists = list(tqdm(pool.imap(worker_func, seeds), total=len(seeds)))

    # 4. STATISTICS
    print("\n[3/3] Final Forensics (Solar Spin Alignment)")
    null_dists = np.array(null_dists)
    
    # How many random cold spots are CLOSER to the Sun's axis than the real one?
    better_nulls = np.sum(null_dists < target_metric)
    p_value = better_nulls / N_SIMS
    
    print(f" -> Nulls Closer to Solar Spin: {better_nulls}/{N_SIMS}")
    print(f" -> P-Value: {p_value:.3f}")

    if p_value < 0.05:
        print("\n[RESULT] POTENTIAL SIGNAL DETECTED.")
    else:
        print("\n[RESULT] Null Hypothesis Consistent.")

if __name__ == "__main__":
    mp.set_start_method('fork') 
    path = sys.argv[1] if len(sys.argv) > 1 else None
    main(path)
