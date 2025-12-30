#!/usr/bin/env python3
"""
Phase II: The "Anatomy" Check (Spectrum Dependency)
Hypothesis: Does the "Low Quadrupole" of the real sky make alignment more likely?
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import numpy as np
import healpy as hp
import multiprocessing as mp
from functools import partial

from ouroboros import ingestion
from ouroboros.engines import harmonics

# -----------------------------------------------------------------------------
# 1. SPECTRUM GENERATORS
# -----------------------------------------------------------------------------
def get_real_spectrum(map_path, lmax=64):
    """Extracts the actual C_l from the provided map."""
    print(f" -> Extracting Real C_l from {map_path}...")
    data_map = ingestion.load_map(map_path)
    # Remove monopole/dipole to get clean Cl
    map_clean = hp.remove_dipole(data_map, verbose=False)
    cls = hp.anafast(map_clean, lmax=lmax)
    return cls

def get_theoretical_spectrum(lmax=64):
    """
    Generates a Scale-Invariant (Harrison-Zel'dovich) spectrum.
    D_l = l(l+1)C_l / 2pi = Constant
    Therefore C_l ~ 1 / (l * (l+1))
    """
    print(" -> Generating Theoretical C_l (Scale Invariant)...")
    ls = np.arange(lmax + 1)
    cls = np.zeros(lmax + 1)
    
    # Avoid division by zero at l=0
    cls[1:] = 1.0 / (ls[1:] * (ls[1:] + 1.0))
    
    # Normalize to match the amplitude of the Real Spectrum at l=10
    # (Just so the map units aren't wildly different, though angle is scale-independent)
    return cls

# -----------------------------------------------------------------------------
# 2. WORKER
# -----------------------------------------------------------------------------
def worker_spectrum(seed, cls, nside):
    """Generates a random universe from a specific C_l."""
    np.random.seed(seed)
    
    # Generate map from Cl (Isotropic phases)
    random_map = hp.synfast(cls, nside, lmax=64, verbose=False, new=True)
    
    # Check Alignment
    try:
        results = harmonics.analyze_axis_of_evil(random_map)
        return results['angle_23']
    except Exception:
        return 90.0

# -----------------------------------------------------------------------------
# 3. RUNNER
# -----------------------------------------------------------------------------
def run_track(track_name, cls, nside, n_sims=1000):
    print(f"\n--- Running Track: {track_name} ---")
    
    worker = partial(worker_spectrum, cls=cls, nside=nside)
    seeds = [50000 + i for i in range(n_sims)]
    
    with mp.Pool(10) as pool:
        angles = list(pool.imap(worker, seeds))
        
    angles = np.array(angles)
    
    # Define "Evil" as < 10 degrees
    THRESHOLD = 10.0
    n_evil = np.sum(angles < THRESHOLD)
    prob = n_evil / n_sims
    
    return prob

def main(path):
    print("="*60)
    print("PHASE II: SPECTRUM DEPENDENCY CHECK")
    print("Does a suppressed Quadrupole force alignment?")
    print("="*60)
    
    # 1. Setup
    real_cls = get_real_spectrum(path)
    flat_cls = get_theoretical_spectrum()
    
    # Normalize Flat to Real for sanity (match l=2 power roughly)
    # Actually, we WANT Flat to have "Normal" power, so we normalize to l=10
    norm_factor = real_cls[10] / flat_cls[10]
    flat_cls *= norm_factor
    
    # Print Quadrupole Comparison
    print(f" -> Real C_2 (Quadrupole): {real_cls[2]:.2e}")
    print(f" -> Flat C_2 (Quadrupole): {flat_cls[2]:.2e}")
    ratio = real_cls[2] / flat_cls[2]
    print(f" -> Ratio (Real/Flat):     {ratio:.2f} (Real is {ratio:.0%} of expected)")
    
    nside = 64 # Sufficient for low-l
    
    # 2. Run Tracks
    prob_real = run_track("A. Real Spectrum (Low Quadrupole)", real_cls, nside)
    prob_flat = run_track("B. Flat Spectrum (Normal Quadrupole)", flat_cls, nside)
    
    # 3. Verdict
    print("\n" + "="*60)
    print("FINAL RESULTS")
    print(f"Track A (Real C_l): P(Evil) = {prob_real:.3f}")
    print(f"Track B (Flat C_l): P(Evil) = {prob_flat:.3f}")
    print("-" * 30)
    
    if abs(prob_real - prob_flat) < 0.005:
        print("[VERDICT] INDEPENDENT.")
        print("The alignment probability is NOT affected by the Quadrupole power.")
        print("The Axis of Evil is a pure Phase Anomaly.")
    elif prob_real > prob_flat:
        print("[VERDICT] DEPENDENT.")
        print("The Low Quadrupole INCREASES the chance of alignment.")
    else:
        print("[VERDICT] INVERSE.")
        print("A Normal Quadrupole makes alignment MORE likely.")
    print("="*60)

if __name__ == "__main__":
    mp.set_start_method('fork')
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
