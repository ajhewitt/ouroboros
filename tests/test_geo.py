"""
Geometry Engine Audit
Verifies Cold Spot localization and nodal distance calculations.
"""

import pytest
import numpy as np
import healpy as hp
from ouroboros.engines import geo
from ouroboros import config

def test_cold_spot_localization():
    """
    Audit 1: Can the engine find a fake Cold Spot injected into a map?
    """
    nside = 64
    npix = hp.nside2npix(nside)
    # Create blank map (warm)
    map_data = np.zeros(npix) + 100.0
    
    # Inject a 'Cold Spot' at a known location
    # Location: Galactic l=200, b=-50 (Southern Hemisphere)
    target_l, target_b = 200.0, -50.0
    target_idx = hp.ang2pix(nside, target_l, target_b, lonlat=True)
    
    # Make it very cold
    map_data[target_idx] = -100.0
    
    # Run Finder
    ra_found, dec_found = geo.find_cold_spot(map_data, smooth_fwhm_deg=0.0)
    
    # Convert found RA/Dec back to Galactic to compare
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra=ra_found*u.deg, dec=dec_found*u.deg, frame='icrs')
    gal = c.galactic
    
    # Check proximity (within pixel size ~1 deg)
    assert np.isclose(gal.l.degree, target_l, atol=2.0)
    assert np.isclose(gal.b.degree, target_b, atol=2.0)

def test_nodal_math():
    """
    Audit 2: Verify distance calculations to Ecliptic Poles.
    """
    # Hypothetical Spot exactly at the South Ecliptic Pole
    # The SEP is opposite the NEP.
    # NEP RA ~ 270 (-90), Dec ~ 66.5
    # SEP RA ~ 90, Dec ~ -66.5
    
    # Let's just use the NEP for simplicity of math
    # If the spot is AT the NEP, dist_nep should be 0.
    ra_test = config.ECLIPTIC_POLE_RA_DEG
    dec_test = config.ECLIPTIC_POLE_DEC_DEG
    
    metrics = geo.check_nodal_alignment(ra_test, dec_test)
    
    assert np.isclose(metrics['dist_nep'], 0.0, atol=0.1)
    assert np.isclose(metrics['dist_sep'], 180.0, atol=0.1)
