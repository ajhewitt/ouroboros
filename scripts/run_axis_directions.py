#!/usr/bin/env python3
"""
Phase V: The Directional Tie-Breaker
Hypothesis: If l=2/3 and l=5/6 are systematic, they should point to the same sky location.
"""
import sys
import numpy as np
import healpy as hp
from ouroboros import ingestion
from ouroboros.engines import harmonics

def vec2dir(vec):
    """Converts 3D cartesian vector to RA/Dec."""
    # Healpy vec2ang returns theta (colat), phi (lon)
    theta, phi = hp.vec2ang(vec)
    
    # Force conversion to simple scalar floats (extract from numpy array if needed)
    dec = 90.0 - np.degrees(np.array(theta).item())
    ra = np.degrees(np.array(phi).item())
    
    return float(ra), float(dec)

def main(path):
    print("="*60)
    print("PHASE V: DIRECTIONAL CONSISTENCY CHECK")
    print("Do the primary and secondary axes point to the same location?")
    print("="*60)
    
    # 1. Load Data
    data_map = ingestion.load_map(path)
    map_clean = hp.remove_dipole(data_map, verbose=False)
    alm = hp.map2alm(map_clean, lmax=10)
    
    # 2. Extract Vectors
    axes = {}
    for l in [2, 3, 5, 6]:
        vec = harmonics.get_principal_axis(alm, l)
        axes[l] = vec
        ra, dec = vec2dir(vec)
        print(f" -> l={l} Axis: RA={ra:.1f}, Dec={dec:.1f}")

    # 3. Compare Pair Averages
    # Since l=2 and l=3 are aligned, we can average them to get the "Evil Axis 1"
    # Note: Axes are headless, so we ensure dot product is positive before averaging
    if np.dot(axes[2], axes[3]) < 0: axes[3] *= -1
    axis_evil_1 = (axes[2] + axes[3])
    axis_evil_1 /= np.linalg.norm(axis_evil_1)
    
    if np.dot(axes[5], axes[6]) < 0: axes[6] *= -1
    axis_evil_2 = (axes[5] + axes[6])
    axis_evil_2 /= np.linalg.norm(axis_evil_2)
    
    # 4. Measure Separation
    dot = np.abs(np.dot(axis_evil_1, axis_evil_2))
    separation = np.degrees(np.arccos(dot))
    
    ra1, dec1 = vec2dir(axis_evil_1)
    ra2, dec2 = vec2dir(axis_evil_2)
    
    print("\n[RESULTS]")
    print(f"Primary Axis (l=2,3):   RA={ra1:.1f}, Dec={dec1:.1f}")
    print(f"Secondary Axis (l=5,6): RA={ra2:.1f}, Dec={dec2:.1f}")
    print(f"Separation Angle:       {separation:.2f} deg")
    
    print("\n[VERDICT]")
    if separation < 20:
        print(" -> SYSTEMATIC ERROR CONFIRMED.")
        print("    Both anomalies point to the same patch of sky.")
        print("    (Likely Ecliptic Scan or Galactic Mask artifact).")
    else:
        print(" -> UNCORRELATED FLUKES.")
        print("    The anomalies point in totally different directions.")
        print("    The 'Axis of Evil' is just random noise clustering.")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
