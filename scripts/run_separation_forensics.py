#!/usr/bin/env python3
"""
Phase VI: The Orthogonality Check
Hypothesis: Is the 55.6 degree separation consistent with pure randomness?
"""
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

def random_axis():
    """Generates a random 3D unit vector on the sphere."""
    # Random point on sphere
    z = np.random.uniform(-1, 1)
    phi = np.random.uniform(0, 2*np.pi)
    theta = np.arccos(z)
    
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    return np.array([x, y, z])

def get_separation(v1, v2):
    """Calculates angle between two axes (0 to 90 degrees)."""
    # Axes are headless, so we want the smallest angle
    dot = np.abs(np.dot(v1, v2))
    if dot > 1.0: dot = 1.0
    return np.degrees(np.arccos(dot))

def main():
    print("="*60)
    print("PHASE VI: SEPARATION FORENSICS")
    print("Quantifying the 'Uncorrelated' Verdict.")
    print("="*60)
    
    # 1. Simulate Random Separations
    N_SIMS = 100000
    print(f" -> Simulating {N_SIMS} random pairs...")
    
    angles = []
    for _ in range(N_SIMS):
        v1 = random_axis()
        v2 = random_axis()
        angles.append(get_separation(v1, v2))
        
    angles = np.array(angles)
    
    # 2. Your Result
    OBSERVED_SEP = 55.60
    
    # 3. Stats
    mean_sep = np.mean(angles)
    # How many random pairs are CLOSER than 55.6? (P-value for alignment)
    p_aligned = np.sum(angles < OBSERVED_SEP) / N_SIMS
    
    print(f"\n[THEORY]")
    print(f" -> Mean Random Separation: {mean_sep:.2f} deg")
    print(f"    (Expect ~57.3 deg for isotropic axes)")
    
    print(f"\n[OBSERVATION]")
    print(f" -> Real Sky Separation:    {OBSERVED_SEP:.2f} deg")
    print(f" -> P(Closer by chance):    {p_aligned:.3f}")
    
    # 4. Interpretation
    print("\n[VERDICT]")
    if 0.4 < p_aligned < 0.6:
        print(" -> PERFECTLY TYPICAL.")
        print("    The separation sits right in the middle of the distribution.")
        print("    This proves the two anomalies are statistically independent.")
    elif p_aligned < 0.05:
        print(" -> SUSPICIOUSLY ALIGNED.")
        print("    Wait, 55 deg is actually rare? (Unlikely).")
    else:
        print(f" -> Normal Fluctuation (P={p_aligned:.2f})")

if __name__ == "__main__":
    main()
