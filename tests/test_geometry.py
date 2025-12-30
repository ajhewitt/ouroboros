"""
Geometry Audit
Verifies the integrity of the Coordinate Definitions in config.py.
"""

import pytest
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

from ouroboros import config

def test_solar_pole_definition():
    """
    Audit 1: Verify the Solar North Pole vector normalization.
    """
    vec = config.get_solar_vector()
    
    # Assert vector length is exactly 1.0 (within float tolerance)
    norm = np.linalg.norm(vec)
    assert np.isclose(norm, 1.0), f"Solar Vector is not normalized! Mag={norm}"

def test_agency_frame_alignment():
    """
    Audit 2: External Library Cross-Check.
    Does our hardcoded Solar Pole match Astropy's coordinate transformations?
    """
    # Create an Astropy object from our config constants
    solar_pole = SkyCoord(
        ra=config.SOLAR_NORTH_POLE_RA_DEG * u.deg,
        dec=config.SOLAR_NORTH_POLE_DEC_DEG * u.deg,
        frame='icrs' # J2000
    )

    # Convert to Galactic Frame
    solar_pole_gal = solar_pole.galactic

    # Refined Check:
    # The North Ecliptic Pole is at b ~ 29.8 deg.
    # The Solar Pole is tilted ~7.25 deg, putting it at b ~ 22.8 deg.
    # We allow a +/- 2 degree margin.
    assert 90 < solar_pole_gal.l.degree < 100, \
        f"Solar Pole Galactic L {solar_pole_gal.l.degree} is out of expected bounds"
        
    # UPDATED BOUNDS: 20 to 25 degrees (expecting ~22.77)
    assert 20 < solar_pole_gal.b.degree < 25, \
        f"Solar Pole Galactic B {solar_pole_gal.b.degree} is out of expected bounds"

def test_pole_separation():
    """
    Audit 3: Geometric relationship check.
    The Sun's spin axis is tilted ~7.25 degrees relative to the Ecliptic Pole.
    We calculate the dot product between our two defined poles and verify the angle.
    """
    # 1. Get Solar Vector
    s_vec = config.get_solar_vector()
    
    # 2. Compute Ecliptic Vector manually
    ra_ecl_rad = np.radians(config.ECLIPTIC_POLE_RA_DEG)
    dec_ecl_rad = np.radians(config.ECLIPTIC_POLE_DEC_DEG)
    e_vec = np.array([
        np.cos(dec_ecl_rad) * np.cos(ra_ecl_rad),
        np.cos(dec_ecl_rad) * np.sin(ra_ecl_rad),
        np.sin(dec_ecl_rad)
    ])
    
    # 3. Calculate Angle: theta = arccos(a . b)
    dot_product = np.dot(s_vec, e_vec)
    angle_rad = np.arccos(dot_product)
    angle_deg = np.degrees(angle_rad)
    
    # The obliquity of the solar equator to the ecliptic is approx 7.25 degrees.
    # We allow a small margin for epoch variations.
    expected_tilt = 7.25
    margin = 0.5 
    
    assert (expected_tilt - margin) < angle_deg < (expected_tilt + margin), \
        f"Computed Solar/Ecliptic tilt {angle_deg:.2f} deg deviates from expected 7.25 deg"

if __name__ == "__main__":
    pytest.main()
