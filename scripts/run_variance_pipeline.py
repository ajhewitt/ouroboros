#!/usr/bin/env python3
"""
Phase IV-C: Ecliptic Variance (The Shielding Test)
Checks if the Ecliptic Plane is anomalously "Quiet" (Low Variance).
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
from astropy.coordinates import SkyCoord
import astropy.units as u

from ouroboros import ingestion
from ouroboros.validation import nulling

def get_ecliptic_variance(map_data):
    """
    Calculates the variance of CMB temperature pixels within 
    +/- 20 degrees of the Ecliptic Plane.
    """
    nside = hp.npix2nside(len(map_data))
    
    # 1. Define the Ecliptic Strip Mask
    # Get pixel coordinates
    ipix = np.arange(len(map_data))
    theta, phi = hp.pix2ang(nside, ipix)
    
    # Convert to Galactic then Ecliptic to find latitude
    # (Doing this via astropy is accurate but slow for full maps, 
    #  so we use a rotation matrix approximation for speed in the loop, 
    #  OR we pre-compute the mask indices ONCE since the strip is fixed relative to the map?
    #  Actually, for Null tests, we rotate the MAP, so the Strip stays fixed.)
    
    # We need to know which pixels are in the Ecliptic Strip in the LAB FRAME.
    # This calculation only needs to happen once.
    pass 

def _precompute_ecliptic_indices(nside):
    """Generates indices for the Ecliptic Belt."""
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    
    # Healpy rotates to Galactic by default. We need to check Ecliptic Lat.
    r = hp.Rotator(coord=['G', 'E']) # Galactic to Ecliptic
    theta_ecl, _ = r(theta, phi)
    
    # Ecliptic Latitude b_ecl = 90 - theta_ecl
    b_ecl = 90.0 - np.degrees(theta_ecl)
    
    # Define Belt: +/- 20 degrees
    indices = np.where(np.abs(b_ecl) < 20)[0]
    return indices

# Global variable for the worker to access (avoids passing huge arrays)
ECLIPTIC_INDICES = None

def worker_variance(seed, original_map):
    # 1. Rotate Map
    gen = nulling.NullGenerator(original_map, seed=seed)
    _, rotated_map = next(gen.generate_nulls(n_sims=1))
    
    # 2. Extract Ecliptic Indices (Global)
    pixels = rotated_map[ECLIPTIC_INDICES]
    
    # 3. GET GALACTIC LATITUDE FOR THESE PIXELS (To filter the Galaxy)
    # We need the coordinates of the *Strip* pixels.
    # Since the Strip is fixed, we should have pre-computed their galactic lat too.
    # (Simplified for script update: Just calculate it here, it's fast enough for 1000 sims)
    
    npix = len(rotated_map)
    nside = hp.npix2nside(npix)
    theta, phi = hp.pix2ang(nside, ECLIPTIC_INDICES)
    # Theta is co-latitude (0 at North Pole). Galactic Lat b = 90 - theta_deg
    b_gal = 90.0 - np.degrees(theta)
    
    # 4. JACKKNIFE: Harsh Galactic Cut
    # Only keep pixels where Galactic Latitude is > 40 or < -40
    # This removes the Milky Way entirely.
    clean_mask = np.abs(b_gal) > 40
    
    clean_pixels = pixels[clean_mask]
    
    # Filter Healpy mask
    valid_pixels = clean_pixels[clean_pixels > -1e20]
    
    if len(valid_pixels) == 0:
        return 0.0
        
    return np.var(valid_pixels)

def main(map_path):
    global ECLIPTIC_INDICES
    CORES = 10
    N_SIMS = 1000
    
    print("="*60)
    print("PHASE IV-C: THE SHIELDING HYPOTHESIS")
    print("Testing if the Ecliptic Plane is anomalously quiet.")
    print("="*60)

    # 1. LOAD DATA
    data_map = ingestion.load_map(map_path)
    nside = hp.npix2nside(len(data_map))
    
    # 2. DEFINE STRIP
    print(" -> Pre-computing Ecliptic Belt indices...")
    ECLIPTIC_INDICES = _precompute_ecliptic_indices(nside)
    
    # 3. REAL VARIANCE
    real_pixels = data_map[ECLIPTIC_INDICES]
    # Filter out mask (assuming Healpy standard unseen pixel)
    real_valid = real_pixels[real_pixels > -1e20]
    real_variance = np.var(real_valid)
    
    print(f" -> Real Ecliptic Variance: {real_variance:.2f} uK^2")

    # 4. NULL TEST
    print(f"\n[2] Running {N_SIMS} Null Simulations...")
    worker_func = partial(worker_variance, original_map=data_map)
    seeds = [2025 + i for i in range(N_SIMS)]
    
    with mp.Pool(processes=CORES) as pool:
        null_vars = list(tqdm(pool.imap(worker_func, seeds), total=len(seeds)))
        
    # 5. STATS
    null_vars = np.array(null_vars)
    
    # Hypothesis: We expect Real Variance to be LOWER than random (Shielding)
    lower_nulls = np.sum(null_vars < real_variance)
    p_value = lower_nulls / N_SIMS
    
    print("\n[3] FINAL FORENSICS")
    print(f" -> Nulls with LOWER variance: {lower_nulls}/{N_SIMS}")
    print(f" -> P-Value (Low Variance): {p_value:.3f}")
    
    if p_value < 0.05:
        print("\n[RESULT] SIGNAL DETECTED.")
        print("The Ecliptic Plane is anomalously quiet.")
    elif p_value > 0.95:
        print("\n[RESULT] INVERSE SIGNAL.")
        print("The Ecliptic Plane is anomalously NOISY.")
    else:
        print("\n[RESULT] Null Hypothesis Consistent.")

if __name__ == "__main__":
    mp.set_start_method('fork') 
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
