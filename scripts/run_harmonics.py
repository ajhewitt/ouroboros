#!/usr/bin/env python3
import os
import sys
import numpy as np
from ouroboros import ingestion
from ouroboros.engines import harmonics

def main(map_path):
    print("="*60)
    print("PHASE IV: THE AXIS OF EVIL FORENSICS")
    print("="*60)
    
    # Load Real Map
    data_map = ingestion.load_map(map_path)
    
    # Analyze
    results = harmonics.analyze_axis_of_evil(data_map)
    
    print("\n[1] INTERNAL ALIGNMENT (The 'Evil' Test)")
    print(f" -> Angle between Quadrupole (l=2) and Octopole (l=3): {results['angle_23']:.2f} deg")
    print("    (Random expectation is ~60 deg. < 20 deg is anomalous.)")
    
    print("\n[2] SOLAR ALIGNMENT (The 'PbC' Test)")
    print(f" -> Quadrupole vs Ecliptic Pole: {results['angle_2_ecl']:.2f} deg")
    print(f" -> Octopole   vs Ecliptic Pole: {results['angle_3_ecl']:.2f} deg")
    print(f" -> Quadrupole vs Equinoxes:     {results['angle_2_eq']:.2f} deg")
    
    # Interpretation
    print("\n[3] INTERPRETATION")
    if results['angle_23'] < 20:
        print(" -> [CONFIRMED] The Axis of Evil exists (l=2 and l=3 are aligned).")
    else:
        print(" -> [FALSIFIED] The Axis of Evil appears to be a statistical fluke.")
        
    if abs(results['angle_2_ecl'] - 90) < 10:
        print(" -> [ALERT] The Axis lies IN the Ecliptic Plane.")
    elif results['angle_2_ecl'] < 10:
        print(" -> [ALERT] The Axis points AT the Ecliptic Pole.")
    else:
        print(" -> No obvious Ecliptic correlation.")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data/raw/planck/smica.fits"
    main(path)
