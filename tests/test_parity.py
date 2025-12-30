"""
Parity Logic Audit
Verifies that the engine correctly separates Even and Odd parity modes.
"""

import pytest
import numpy as np
import healpy as hp
from ouroboros.engines import parity
from ouroboros import config

def generate_pure_parity_map(l_target, nside=64):
    """Helper to create a map with power ONLY at a specific multipole l."""
    npix = hp.nside2npix(nside)
    alm_size = hp.Alm.getsize(config.L_MAX)
    alm = np.zeros(alm_size, dtype=np.complex128)
    
    # Set power for l=l_target, m=0 (simplest mode)
    idx = hp.Alm.getidx(config.L_MAX, l_target, 0)
    alm[idx] = 1.0 + 0j
    
    return hp.alm2map(alm, nside)

def test_pure_even_parity():
    """
    Audit 1: Pure Even Input (l=2) must yield P(n) = +1.0
    """
    # Create map with only Quadrupole (l=2)
    map_even = generate_pure_parity_map(2)
    
    score = parity.calculate_point_parity(map_even)
    
    # Allow small float tolerance due to map2alm/alm2map transforms
    assert np.isclose(score, 1.0, atol=1e-3), \
        f"Pure Even map (l=2) scored {score}, expected +1.0"

def test_pure_odd_parity():
    """
    Audit 2: Pure Odd Input (l=3) must yield P(n) = -1.0
    """
    # Create map with only Octupole (l=3)
    map_odd = generate_pure_parity_map(3)
    
    score = parity.calculate_point_parity(map_odd)
    
    assert np.isclose(score, -1.0, atol=1e-3), \
        f"Pure Odd map (l=3) scored {score}, expected -1.0"

def test_mixed_parity():
    """
    Audit 3: Equal mixture of l=2 and l=3 should yield near 0.0
    """
    map_even = generate_pure_parity_map(2)
    map_odd = generate_pure_parity_map(3)
    map_mixed = map_even + map_odd
    
    score = parity.calculate_point_parity(map_mixed)
    
    # It won't be exactly 0 because power scales with l(l+1)
    # l=2 power factor ~ 6, l=3 power factor ~ 12.
    # So odd dominates slightly. We just check it's not +1 or -1.
    assert -0.8 < score < 0.8, \
        f"Mixed map scored {score}, expected intermediate value"

def test_kinematic_isolation():
    """
    Audit 4: Ensure Dipole (l=1) does not contaminate the result.
    """
    # Create pure Even map
    map_base = generate_pure_parity_map(2)
    
    # Add massive Dipole
    map_dirty = map_base + 100 * generate_pure_parity_map(1)
    
    # The engine should clean the dipole. Result should still be +1.0 (Even).
    score = parity.calculate_point_parity(map_dirty)
    
    assert np.isclose(score, 1.0, atol=1e-3), \
        "Dipole contamination detected! Engine failed Kinematic Isolation."
