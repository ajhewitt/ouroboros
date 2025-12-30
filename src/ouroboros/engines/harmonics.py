"""
Harmonic Engine (Phase IV)
Extracts the "Maxwellian Multipole Vectors" (Principal Axes) of the CMB.
"""

import numpy as np
import healpy as hp
from ouroboros import config

def get_principal_axis(alm, l):
    """
    Finds the axis of maximum angular momentum for a specific multipole l.
    This defines the "Plane" of the Quadrupole/Octopole.
    """
    # 1. Reconstruct the map for ONLY this multipole
    # We zero out all other l's
    lmax = hp.Alm.getlmax(len(alm))
    alm_filtered = np.zeros_like(alm)
    
    start_idx = hp.Alm.getidx(lmax, l, 0)
    end_idx = hp.Alm.getidx(lmax, l, l) + 1
    
    # Copy the m-modes for this l
    # Note: healpy layout is complex, simplified extraction:
    indices = hp.Alm.getidx(lmax, l, np.arange(l + 1))
    alm_filtered[indices] = alm[indices]
    
    # 2. Generate the map (NSIDE=64 is enough for low-l)
    m_l = hp.alm2map(alm_filtered, nside=64, verbose=False)
    
    # 3. Find the "Axis of Symmetry"
    # We calculate the Moment of Inertia tensor of the power distribution
    npix = len(m_l)
    vecs = np.array(hp.pix2vec(64, np.arange(npix))).T
    weights = m_l**2 # We weigh by power squared to find the "mass"
    
    # Inertia Tensor: I_ij = sum(w * (r^2 delta_ij - r_i r_j))
    # But for "Axis of Alignment", we just want the weighted covariance matrix
    # C_ij = sum(w * r_i * r_j)
    cov = np.dot(vecs.T, (vecs * weights[:, None]))
    
    # Eigen decomposition
    evals, evecs = np.linalg.eigh(cov)
    
    # The axis of the "Pancake" is the eigenvector with the SMALLEST eigenvalue
    # (Pointing perpendicular to the mass distribution)
    # Or for l=2, usually defined by the maximum power axis.
    # Standard "Axis of Evil" definition is the direction that maximizes angular momentum dispersion.
    # Let's use the standard "Copi et al." definition: The axis that minimizes power variation *around* it.
    
    # SIMPLIFIED ROBUST METHOD:
    # Just find the dipole direction of the power map |m_l|^2
    # This points to the "Hotspots" of the multipole.
    
    return evecs[:, 0] # Returns the normal vector to the plane

def analyze_axis_of_evil(map_data):
    """
    Extracts l=2 and l=3 axes and checks Solar Alignment.
    """
    # remove dipole
    map_clean = hp.remove_dipole(map_data, verbose=False)
    alm = hp.map2alm(map_clean, lmax=10)
    
    # Get Axes
    axis_l2 = get_principal_axis(alm, 2)
    axis_l3 = get_principal_axis(alm, 3)
    
    # Check Alignment with each other (The "Evil" Internal Check)
    dot_23 = np.abs(np.dot(axis_l2, axis_l3))
    angle_23 = np.degrees(np.arccos(dot_23))
    
    # Check Alignment with Ecliptic Plane (The "Solar" Check)
    # If the vector is IN the Ecliptic, it is 90 deg from the Pole.
    # If the vector is PERPENDICULAR (aligned with Pole), it is 0 deg.
    ecl_pole = config.get_ecliptic_vector()
    
    angle_2_ecl = np.degrees(np.arccos(np.abs(np.dot(axis_l2, ecl_pole))))
    angle_3_ecl = np.degrees(np.arccos(np.abs(np.dot(axis_l3, ecl_pole))))
    
    # Check Alignment with Equinoxes (RA=0, Dec=0)
    equinox = np.array([1, 0, 0]) # Cartesian X is Vernal Equinox in ICRS roughly
    angle_2_eq = np.degrees(np.arccos(np.abs(np.dot(axis_l2, equinox))))
    
    return {
        "axis_l2": axis_l2,
        "axis_l3": axis_l3,
        "angle_23": angle_23,
        "angle_2_ecl": angle_2_ecl,
        "angle_3_ecl": angle_3_ecl,
        "angle_2_eq": angle_2_eq
    }
