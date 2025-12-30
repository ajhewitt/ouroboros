#!/usr/bin/env python3
import sys
import numpy as np
import healpy as hp
from ouroboros import ingestion, config
from ouroboros.engines import harmonics
from ouroboros.validation import nulling

def main():
    print("DEBUG: Loading Map...")
    path = "data/raw/planck/smica.fits"
    data_map = ingestion.load_map(path)
    
    print("DEBUG: Generating Single Null...")
    gen = nulling.NullGenerator(data_map, seed=12345)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    print("DEBUG: Applying Hard Mask (|b| < 20)...")
    nside = hp.npix2nside(len(rotated_map))
    npix = hp.nside2npix(nside)
    theta, _ = hp.pix2ang(nside, np.arange(npix))
    b_val = 90.0 - np.degrees(theta)
    
    mask = np.ones(npix)
    mask[np.abs(b_val) < 20] = 0.0
    
    masked_map = rotated_map * mask
    
    # Check for NaNs before processing
    if np.any(np.isnan(masked_map)):
        print(" -> WARNING: Map contains NaNs!")
        # Healpy hates NaNs. We must replace them.
        masked_map[np.isnan(masked_map)] = 0.0
        print(" -> Fixed NaNs.")
        
    print("DEBUG: Running Harmonic Analysis (No Safety Block)...")
    
    # This is where it likely crashes
    results = harmonics.analyze_axis_of_evil(masked_map)
    
    print("\n[SUCCESS]")
    print(f"Angle 2-3: {results['angle_23']:.4f} deg")

if __name__ == "__main__":
    main()
