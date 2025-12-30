"""
Quasar Engine Audit
Verifies vector separation logic and alignment statistics.
"""

import pytest
import numpy as np
from ouroboros.engines import quasar
from ouroboros.validation import shuffling

def test_separation_logic():
    """
    Audit 1: Geometry Check.
    If we have two objects on the X-axis, their separation vector must be X.
    """
    # Object A at x=100, y=0, z=0 (RA=0, Dec=0, r=100)
    # Object B at x=200, y=0, z=0 (RA=0, Dec=0, r=200)
    ra = np.array([0.0, 0.0])
    dec = np.array([0.0, 0.0])
    r = np.array([100.0, 200.0])

    vecs = quasar.get_separation_vectors(ra, dec, r)

    # Result should be aligned with X-axis.
    # The engine computes (A - B) -> (100 - 200) = -100.
    # So we expect [-1, 0, 0].
    expected = np.array([-1.0, 0.0, 0.0])

    assert np.allclose(vecs[0], expected, atol=1e-5), \
        f"Separation vector {vecs[0]} does not match expected X-axis alignment."

def test_perfect_alignment():
    """
    Audit 2: Statistic Check.
    If all vectors align with the axis, score should be 1.0.
    """
    # Create vectors exactly along Z
    vecs = np.array([
        [0, 0, 1],
        [0, 0, -1] # Opposite direction should also count as aligned (abs value)
    ])
    axis = np.array([0, 0, 1])

    score = quasar.correlate_with_axis(vecs, axis)

    assert score == 1.0, f"Perfect alignment scored {score} (expected 1.0)"

def test_random_shuffling():
    """
    Audit 3: Shuffling Check.
    Ensure the 'spin' method actually changes coordinates.
    """
    ra = np.array([10.0, 20.0, 30.0])
    dec = np.array([0.0, 0.0, 0.0])

    ra_new, dec_new = shuffling.shuffle_catalog_vectors(ra, dec, method='spin')

    # Dec should be identical
    assert np.allclose(dec, dec_new)

    # RA should be different (shifted)
    assert not np.allclose(ra, ra_new), "Shuffling failed to scramble RA!"
