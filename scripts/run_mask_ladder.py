#!/usr/bin/env python3
"""
Phase I: The Mask Sensitivity Ladder (Corrected)
Hypothesis: Does a Galactic Mask force random noise to align?
Method: Generates fresh Gaussian Random Fields (synfast) to test purely geometric artifacts.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import numpy as np
import healpy as hp
import multiprocessing as mp
from functools import partial

# We only need harmonics for the measurement
from ouroboros.engines import harmonics
from ouroboros import ingestion

# -----------------------------------------------------------------------------
# 1. UTILS
# -----------------------------------------------------------------------------
def generate_galactic_mask(nside, lat_cut_deg):
    """Generates a binary mask (0 = cut, 1 = keep)."""
    if lat_cut_deg == 0:
        return np.ones(hp.nside2npix(nside))
        
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    b_val = 90.0 - np.degrees(theta)
    
    mask = np.ones(npix)
    mask[np.abs(b_val) < lat_cut_deg] = 0.0
    return mask

# -----------------------------------------------------------------------------
# 2. WORKER
# -----------------------------------------------------------------------------
def worker_pure_noise(seed, cls, nside, mask):
    """
    Generates a FRESH random universe (Gaussian) and tests alignment.
    """
    # 1. Generate Random Map from Power Spectrum (synfast)
    # We use lmax=32 because we only care about l=2,3
    np.random.seed(seed)
    
    # Generate random alm
    # This ensures an isotropic universe (No intrinsic axis)
    random_map = hp.synfast(cls, nside, lmax=64, verbose=False, new=True)
    
    # 2. Apply Mask (Hard Zeroing)
    # This is the "geometric damage"
    masked_map = random_map * mask
    
    # 3. Measure Alignment
    try:
        results = harmonics.analyze_axis_of_evil(masked_map)
        return results['angle_23']
    except Exception:
        return 90.0 # Fail safe

# -----------------------------------------------------------------------------
# 3. RUNNER
# -----------------------------------------------------------------------------
def run_ladder_step(cls, nside, lat_cut, n_sims=500):
    # 1. Create Mask
    mask = generate_galactic_mask(nside, lat_cut)
    f_sky = np.sum(mask) / len(mask)
    
    print(f"--- Mask |b| < {lat_cut} deg (f_sky={f_sky:.1%}) ---")
    
    # 2. Run Nulls
    worker = partial(worker_pure_noise, cls=cls, nside=nside, mask=mask)
    seeds = [10000 + (int(lat_cut)*10000) + i for i in range(n_sims)]
    
    with mp.Pool(10) as pool:
        angles = list(pool.imap(worker, seeds))
        
    angles = np.array(angles)
    
    # 3. Stats
    # We define "Evil" as alignment < 10 degrees (Standard definition is ~10)
    THRESHOLD = 10.0
    n_evil = np.sum(angles < THRESHOLD)
    prob = n_evil / n_sims
    
    return f_sky, prob

def main(path):
    print("="*60)
    print("PHASE I: THE MASK SENSITIVITY LADDER (PURE NOISE)")
    print("Does a mask create an 'Axis of Evil' out of nothing?")
    print("="*60)
    
    # 1. Load Data to get Power Spectrum (Cl)
    print(" -> extracting C_l from real sky...")
    data_map = ingestion.load_map(path)
    nside = hp.npix2nside(len(data_map))
    
    # Remove monopole/dipole before Cl extraction
    map_clean = hp.remove_dipole(data_map, verbose=False)
    cls = hp.anafast(map_clean, lmax=64)
    
    # 2. The Ladder
    cuts = [0, 10, 20, 30, 40]
    results = []
    
    for cut in cuts:
        f_sky, prob = run_ladder_step(cls, nside, cut, n_sims=500)
        results.append((cut, f_sky, prob))
        print(f"    -> P(Align < 10 deg) = {prob:.3f}")

    print("\n" + "="*60)
    print(f"{'Lat Cut':<10} | {'f_sky':<10} | {'P(Evil)':<10}")
    print("-" * 50)
    
    for cut, f_sky, prob in results:
        print(f"{cut:<10} | {f_sky:.1%}   | {prob:.3f}")
    print("="*60)

if __name__ == "__main__":
    mp.set_start_method('fork')
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
