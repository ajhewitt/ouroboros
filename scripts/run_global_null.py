#!/usr/bin/env python3
"""
Phase III: The Global Null (Look-Elsewhere Effect)
Hypothesis: Is the l=2/3 alignment unique, or do all low-l modes align?
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import numpy as np
import healpy as hp
from ouroboros import ingestion
from ouroboros.engines import harmonics

def get_pair_angle(alm, l_a, l_b):
    """Calculates the angle between the principal axes of two multipoles."""
    axis_a = harmonics.get_principal_axis(alm, l_a)
    axis_b = harmonics.get_principal_axis(alm, l_b)
    
    dot = np.abs(np.dot(axis_a, axis_b))
    angle = np.degrees(np.arccos(dot))
    # Alignment is mod 90 or 180? Principal axes are headless.
    # Angle is between 0 and 90.
    if angle > 90: angle = 180 - angle
    return angle

def main(path):
    print("="*60)
    print("PHASE III: THE GLOBAL NULL")
    print("Checking alignment of multipole pairs up to l=10.")
    print("="*60)
    
    # 1. Load Real Sky
    data_map = ingestion.load_map(path)
    map_clean = hp.remove_dipole(data_map, verbose=False)
    alm = hp.map2alm(map_clean, lmax=16)
    
    # 2. Measure Real Angles
    print("\n[1] Real Sky Alignments:")
    print(f"{'Pair':<10} | {'Angle (deg)':<12} | {'Status'}")
    print("-" * 40)
    
    real_angles = {}
    
    # We check pairs 2-3, 3-4, ... 9-10
    for l in range(2, 10):
        l_next = l + 1
        angle = get_pair_angle(alm, l, l_next)
        real_angles[l] = angle
        
        # Simple flag for "Close" (< 20 deg)
        status = "ALIGNED" if angle < 20 else "random"
        if angle < 10: status = "** EVIL **"
        
        print(f"l={l} vs {l_next} | {angle:.2f}        | {status}")

    # 3. Interpretation
    # We already know from Phase I/II that the random expectation for 'Evil' (<10deg) is ~1.9%
    # So finding ONE is rare (1/50). Finding TWO would be impossible (1/2500).
    
    evil_count = sum(1 for a in real_angles.values() if a < 10)
    
    print("\n[2] FORENSICS")
    if evil_count == 1 and real_angles[2] < 10:
        print(" -> RESULT: The l=2/3 alignment is UNIQUE.")
        print("    No other low-l pairs show this behavior.")
        print("    This increases the statistical significance of the Axis of Evil.")
    elif evil_count > 1:
        print(" -> RESULT: Multiple pairs are aligned.")
        print("    This suggests a systematic effect or broad structural anisotropy.")
    else:
        print(" -> RESULT: No pairs are strongly aligned (Wait, did 2-3 disappear?)")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
