![Project Ouroboros](docs/ouroboros.jpg)

**Principal Investigator:** A. Hewitt  
**Status:** COMPLETE (Hypothesis Falsified)  
**Date:** December 30th, 2025

## 1. Project Overview

Project Ouroboros was a computational investigation into the **Parochial by Construction (PbC)** hypothesis. This framework posited that the apparent "fine-tuning" of the Universe's history might be a retro-causal artifact selected to be compatible with the local observer's reference frame (the Solar System).

The pipeline performed a rigorous **Null-Test Audit (N >= 1000)** across three cosmological epochs:
1.  **The CMB (Phase I-IV):** Testing geometric alignments of temperature anomalies.
2.  **The Dark Ages:** Testing kinematic coupling of the Dipole.
3.  **The Reionization Era:** Testing structural alignment of High-Z Quasars.

## 2. Final Results
The investigation systematically ruled out the hypothesis across all domains.

| Investigation | Target Vector | P-Value | Conclusion |
| :--- | :--- | :--- | :--- |
| **A. Parity Mirror** | Solar Spin Axis | **0.890** | **Null** (Counter-Aligned) |
| **B. Cold Spot** | Ecliptic Nodes | **0.614** | **Null** (Random Placement) |
| **C. Quasar Structure** | Solar Spin Axis | **0.221** | **Null** (Random Orientation) |
| **D. Ecliptic Shielding**| Variance | **0.747** | **Null** (Signal vanished after Galactic Cut) |

### The "Axis of Evil" Finding (Phase IV)
We investigated the anomalous alignment between the CMB Quadrupole (l=2) and Octopole (l=3), which are observed to be parallel within ~9 degrees.
* **Solar Alignment:** We found no statistically significant link between this axis and the Solar System (p=0.459).
* **Intrinsic Validity:** A null test (N=2000) revealed that **60.0%** of random universes produce an alignment tighter than 9 degrees when subjected to a standard Galactic Mask.
* **Verdict:** The "Axis of Evil" is confirmed to be a **geometric artifact** of the mask, not a cosmological structure.

---

## 3. Repository Structure

### `/scripts` (Execution Entry Points)
These scripts run specific scientific audits.

#### **Phase I: Initial Ledger Audit**

#### **Phase II: Geometric Audits**
* `run_parity_pipeline.py`: Calculates Point-Parity Asymmetry and checks alignment with Solar Spin/Orbit.
* `run_geo_pipeline.py`: Checks if the Eridanus Supervoid (Cold Spot) is geometrically locked to Ecliptic Nodes.
* `run_quasar_pipeline.py`: (Plan C) Analyzes separation vectors of High-Z Quasars (z > 2.5) relative to Solar Angular Momentum.

#### **Phase III: Kinematic Audits**
* `run_phase3_tomography.py`: Slices the Parity signal into multipole bins (l=2-30, 31-60, 61-100) to check for directional stability (The "Drift" Test).

#### **Phase IV: Harmonic & Variance Audits (The "Axis of Evil")**
* `run_harmonics.py`: Extracts the Principal Axes ("Pancakes") of the l=2 and l=3 modes.
* `run_harmonics_null.py`: Tests if the Octopole (l=3) is statistically locked to the Ecliptic Pole.
* `run_axis_validity.py`: **The "Axis of Evil" Test.** Checks if the internal l=2/3 alignment is anomalous vs. a masked random sky.
* `run_variance_pipeline.py`: Tests the "Shielding Hypothesis" (Is the Ecliptic plane anomalously quiet?).
* `run_jackknife.py`: Re-runs the variance test with a harsh Galactic Cut (|b| > 40 deg) to rule out foregrounds.

#### **Utilities**
* `inspect_columns.py`: Helper to view FITS headers for catalog debugging.

---

### `/src/ouroboros/engines` (Physics Logic)
Core mathematical modules implementing the specific tests.

* `parity.py`: Implements P(n) = C(n) + C(-n) symmetry calculations on Healpix maps.
* `geo.py`: Geodesic distance calculators for void/node alignment.
* `quasar.py`: 3D vector algebra for quasar distribution analysis (includes memory-mapped loading for DESI).
* `harmonics.py`: Eigenvalue decomposition for extracting multipole orientation vectors.

### `/src/ouroboros/validation` (Statistical Controls)
Generators for the Null Hypothesis.

* `nulling.py`: Generates random realizations of the CMB by shuffling a_lm phases while preserving the Power Spectrum (C_l).
* `shuffling.py`: Generates random Quasar catalogs by rotating positions on the sphere.

### `/src/ouroboros` (Core)
* `ingestion.py`: Standardized loading for Planck `SMICA`/`NPIPE` maps (downgrading to NSIDE=64).
* `config.py`: Central repository of Solar System vector constants (Solar Spin, Ecliptic Pole, Dipole).

---

## 4. Usage

To reproduce the primary "Axis of Evil" masking artifact result:

```bash
# 1. Download Planck SMICA map to data/raw/planck/smica.fits
# 2. Run the internal validity check
python scripts/run_axis_validity.py data/raw/planck/smica.fits
```

---

## 5. Acknowledgments

Data provided by the Planck Collaboration (2018 Release) and DESI (Data Release 1). Analysis performed using `healpy` and `astropy`.

---

## 6. License

MIT
