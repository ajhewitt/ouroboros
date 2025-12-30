"""
Ouroboros Configuration Registry
Defines the 'Agency Frame' geometry and pipeline constants.
"""

import numpy as np

# --- PIPELINE SETTINGS ---
# HEALPix resolution. NSIDE=64 is sufficient for Large Scale Structure (L < 100).
NSIDE = 64

# The multipole limit for Parity and Correlation analysis.
# Lower limit L_MIN=2 excludes the Dipole (Kinematic Isolation).
L_MIN = 2
L_MAX = 100

# Monte Carlo parameters for the Null-Test Engine
N_SIMS = 100  # Minimum required for significance testing


# --- AGENCY FRAME GEOMETRY (J2000 Epoch) ---
# Coordinates of the Solar Angular Momentum Vector (The Sun's North Pole).
# Source: Archinal et al. 2011 (IAU) / Carrington Elements
# RA: 286.13 deg, Dec: +63.87 deg
SOLAR_NORTH_POLE_RA_DEG = 286.13
SOLAR_NORTH_POLE_DEC_DEG = 63.87

# Coordinates of the North Ecliptic Pole (NEP).
# RA: 18h 00m 00s (270.0 deg), Dec: +66d 33m 38.55s (66.56 deg)
ECLIPTIC_POLE_RA_DEG = 270.00
ECLIPTIC_POLE_DEC_DEG = 66.56

# The Galactic North Pole (for reference transformations)
GALACTIC_NORTH_POLE_RA_DEG = 192.85948
GALACTIC_NORTH_POLE_DEC_DEG = 27.12825

def get_solar_vector():
    # TEMPORARY SWAP: Testing Kinematic Alignment
    return get_dipole_vector()

    # TEMPORARY SWAP: Testing Ecliptic Alignment
    return get_ecliptic_vector()

    """Returns the Solar North Pole as a normalized Cartesian vector [x, y, z]."""
    ra_rad = np.radians(SOLAR_NORTH_POLE_RA_DEG)
    dec_rad = np.radians(SOLAR_NORTH_POLE_DEC_DEG)

    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.array([x, y, z])

# --- ADDITIONAL VECTORS ---

# The direction of the Solar System's motion relative to the CMB (The Dipole).
# Source: Planck 2018 Solar Dipole
# RA: 167.942 deg, Dec: -6.944 deg (Galactic l=263.99, b=48.26)
CMB_DIPOLE_RA_DEG = 167.942
CMB_DIPOLE_DEC_DEG = -6.944

def get_dipole_vector():
    """Returns the CMB Dipole direction as a normalized vector."""
    ra_rad = np.radians(CMB_DIPOLE_RA_DEG)
    dec_rad = np.radians(CMB_DIPOLE_DEC_DEG)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.array([x, y, z])

def get_ecliptic_vector():
    """Returns the North Ecliptic Pole as a normalized vector."""
    ra_rad = np.radians(ECLIPTIC_POLE_RA_DEG)
    dec_rad = np.radians(ECLIPTIC_POLE_DEC_DEG)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.array([x, y, z])
