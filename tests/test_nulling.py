"""
Null-Test Engine Audit
Verifies that the Monte Carlo generator preserves physical properties.
"""

import pytest
import numpy as np
import healpy as hp
from ouroboros.validation import nulling
from ouroboros import config

def test_rotation_conservation():
    """
    Audit 1: Physics Conservation.
    Rotating the Universe should not change its total energy (Power Spectrum).
    """
    # 1. Create a synthetic map with known power (e.g., l=10)
    nside = 64
    npix = hp.nside2npix(nside)
    idx = hp.Alm.getidx(config.L_MAX, 10, 0)
    alm = np.zeros(hp.Alm.getsize(config.L_MAX), dtype=np.complex128)
    alm[idx] = 10.0
    original_map = hp.alm2map(alm, nside)
    
    # Calculate original Power Spectrum (Cl)
    cl_orig = hp.anafast(original_map, lmax=config.L_MAX)
    
    # 2. Generate a Rotated Null
    generator = nulling.NullGenerator(original_map, seed=42)
    _, rotated_map = next(generator.generate_nulls(n_sims=1))
    
    # Calculate new Power Spectrum
    cl_new = hp.anafast(rotated_map, lmax=config.L_MAX)
    
    # 3. Assert Conservation
    # We check l=10 specifically. 
    # Note: Pixel rotation involves interpolation, so small leakage is expected.
    # We allow 5% tolerance for pixelation effects at NSIDE=64.
    power_diff = np.abs(cl_new[10] - cl_orig[10]) / cl_orig[10]
    
    assert power_diff < 0.05, \
        f"Rotation corrupted the Power Spectrum! Diff: {power_diff:.2%}"

def test_randomness():
    """
    Audit 2: Ensure rotations are actually changing the map.
    """
    nside = 64
    original_map = np.random.randn(hp.nside2npix(nside))
    
    generator = nulling.NullGenerator(original_map, seed=42)
    _, rotated_map = next(generator.generate_nulls(n_sims=1))
    
    # If the map is identical, the rotator is broken (Identity matrix).
    correlation = np.corrcoef(original_map, rotated_map)[0, 1]
    
    assert correlation < 0.99, "Rotated map is identical to original! (Rotation failed)"
