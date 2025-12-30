#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from ouroboros import ingestion
from ouroboros.engines import harmonics

def plot_axes_map(path_map, out_path):
    """Plot l=2,3 (Primary) and l=5,6 (Secondary) axes on the sky."""
    print(f"Generating Sky Map: {out_path}...")
    
    # 1. Load & Process Map
    data_map = ingestion.load_map(path_map)
    map_clean = hp.remove_dipole(data_map, verbose=False)
    alm = hp.map2alm(map_clean, lmax=10)

    # 2. Extract Vectors
    vecs = {}
    for l in [2, 3, 5, 6]:
        vecs[l] = harmonics.get_principal_axis(alm, l)

    # 3. Setup Plot
    plt.figure(figsize=(10, 6))
    
    # Background Map (RdBu_r: Red=Hot, Blue=Cold)
    hp.mollview(map_clean, title="CMB Principal Axes Locations", 
                unit=r"$\Delta T$ [$\mu K$]", cmap='RdBu_r', hold=True)
    
    # 4. Plot Axes
    
    # Primary Axis (l=2, 3) -> BRIGHT GREEN (High contrast on Red/Blue)
    # (Red was okay, but Green pops better against the Red regions of the map)
    # actually let's stick to RED for Primary (classic "Evil") and BLACK for Secondary.
    
    theta2, phi2 = hp.vec2ang(vecs[2])
    theta3, phi3 = hp.vec2ang(vecs[3])
    hp.projscatter(theta2, phi2, marker='X', color='red', s=150, label='l=2 Axis (Quadrupole)')
    hp.projscatter(theta3, phi3, marker='o', facecolors='none', edgecolors='red', s=150, lw=2.5, label='l=3 Axis (Octopole)')

    # Secondary Axis (l=5, 6) -> BLACK (Maximum Contrast)
    theta5, phi5 = hp.vec2ang(vecs[5])
    theta6, phi6 = hp.vec2ang(vecs[6])
    hp.projscatter(theta5, phi5, marker='X', color='black', s=150, label='l=5 Axis')
    hp.projscatter(theta6, phi6, marker='o', facecolors='none', edgecolors='black', s=150, lw=2.5, label='l=6 Axis')

    plt.legend(loc='lower right', frameon=True, framealpha=0.9)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_separation_hist(out_path):
    """Plot the PDF of random axis separations."""
    print(f"Generating Histogram: {out_path}...")
    
    N_SIMS = 100000
    angles = []
    for _ in range(N_SIMS):
        v1 = np.random.randn(3); v1 /= np.linalg.norm(v1)
        v2 = np.random.randn(3); v2 /= np.linalg.norm(v2)
        dot = np.abs(np.dot(v1, v2))
        angles.append(np.degrees(np.arccos(dot)))
    
    angles = np.array(angles)
    
    plt.figure(figsize=(8, 5))
    count, bins, ignored = plt.hist(angles, 50, density=True, alpha=0.6, color='gray', label='Random Isotropic PDF')
    
    # Add Signal Line
    obs_val = 55.60
    plt.axvline(obs_val, color='r', linestyle='--', linewidth=2, label=f'Observed ($55.6^\circ$)')
    
    # Theory Curve
    x = np.linspace(0, 90, 100)
    pdf = np.sin(np.radians(x)) * (np.pi / 180.0)
    plt.plot(x, pdf, 'k-', lw=1.5, label='Theory (sin $\\theta$)') 

    plt.xlabel('Separation Angle (deg)')
    plt.ylabel('Probability Density')
    plt.title('Directional Coherence Check (l=2/3 vs l=5/6)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(out_path, dpi=300)
    plt.close()

if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    plot_axes_map(path, "docs/tex/fig_axes.png")
    plot_separation_hist("docs/tex/fig_separation.png")
