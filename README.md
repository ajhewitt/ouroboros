![Project Ouroboros](docs/ouroboros.jpg)

**Cosmological History Selection by Local Agency**
*Principal Investigator: A. Hewitt | The Ante Institute*

## Overview
This repository contains the computational pipeline for **Phase II** of the Universal Ledger investigation. It implements a rigor-first approach to testing the **Biological Selection Hypothesis**—the proposal that the "Axis of Evil" (CMB Quadrupole/Octupole alignment) and other cosmological anomalies are not measurement errors, but selection artifacts soft-locked to the observer's local frame (The Solar System).

## Repo Layout

```text
ouroboros/
├── data/                   # IGNORED by git (local storage for FITS/catalogs)
│   ├── raw/                # Planck NPIPE, SMICA maps, SDSS catalogs
│   └── processed/          # Intermediate outputs (cached masks, filtered tables)
│
├── docs/                   # PDF documentation
│   └── tex/                # LaTeX documentation
│       └── research_plan.tex # The "Project Ouroboros: Phase II" manuscript
│
├── notebooks/              # Jupyter notebooks for exploratory analysis & visualization
│   ├── 01_parity_scan.ipynb
│   ├── 02_cold_spot_geo.ipynb
│   └── 03_quasar_corr.ipynb
│
├── src/                    # Source code package
│   └── ouroboros/
│       ├── __init__.py
│       ├── config.py       # Hardcoded constants (Solar Poles, coordinate frames)
│       ├── ingestion.py    # Map loading (healpy) and catalog parsing (astropy)
│       ├── tools/          # Shared mathematical utilities
│       │   ├── math.py     # Spherical geometry helpers
│       │   └── stats.py    # Statistical functions (P(n), significance)
│       ├── engines/        # The three core investigation tracks
│       │   ├── parity.py   # Plan A: Point-Parity logic [cite: 20]
│       │   ├── geo.py      # Plan B: Geodesic distance calculators [cite: 26]
│       │   └── quasar.py   # Plan C: Vector alignment logic [cite: 35]
│       └── validation/     # The Null-Test Filter
│           ├── nulling.py  # Monte Carlo generators (rotations) [cite: 15]
│           └── shuffling.py# Catalog shuffling logic [cite: 16]
│
├── scripts/                # Executable CLI scripts to run full pipelines
│   ├── run_parity_pipeline.py
│   └── run_null_tests.py
│
├── tests/                  # Unit tests (pytest)
│   ├── test_geometry.py
│   └── test_ingestion.py
│
├── .gitignore              # Standard Python + Data exclusions
├── pyproject.toml          # Dependency management & build config
└── README.md               # Documentation & Manifest
```

## Pipeline Architecture

The codebase is structured around three investigation tracks backed by a shared "Null-Test Engine":

### 1. The Parity Mirror (Plan A)
* **Goal:** Test if the anomalous CMB Parity Asymmetry aligns with the **Solar Angular Momentum Vector** rather than the Ecliptic.
* **Module:** `src/ouroboros/engines/parity.py`
* **Methodology:** Directional scanning of the Point-Parity statistic ($P(n)$) for multipoles $l=2$ to $100$.

### 2. The Cold Spot Shadow (Plan B)
* **Goal:** Determine if the **Cold Spot (Eridanus Supervoid)** lies at a harmonic node relative to Solar geometry.
* **Module:** `src/ouroboros/engines/geo.py`
* **Methodology:** Geodesic distance calculation from the Cold Spot center to Ecliptic Poles and Equinox nodes.

### 3. High-Z Quasar Soft-Lock (Plan C)
* **Goal:** Search for "Soft-Lock" alignments in the Reionization era ($z > 6$) that may have faded in lower-redshift structures.
* **Module:** `src/ouroboros/engines/quasar.py`
* **Methodology:** Correlating Quasar Polarization/Separation Vectors with the CMB Quadrupole axis.

## Validation Protocol: The Null-Test Engine
To rule out False Positives (specifically Zodiacal light or window-function aliasing), all findings must pass the **Forensic Filter**:

1.  **Kinematic Isolation:** Explicit removal of Monopole ($l=0$) and Dipole ($l=1$) modes.
2.  **Foreground Veto:** Signals must persist in both Raw and Component-Separated (SMICA/SEVEM) maps.
3.  **Geometric Nulling (Monte Carlo):**
    * **Rotation:** Significance is determined by rotating the sky $N \ge 100$ times relative to the Solar System frame.
    * **Shuffling:** For catalogs, positions are randomized while preserving mask density.

## Setup & Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/ante-institute/ouroboros.git](https://github.com/ante-institute/ouroboros.git)
    cd ouroboros
    ```

2.  **Install dependencies:**
    ```bash
    pip install .
    # OR for development
    pip install -e .[dev]
    ```

3.  **Data Ingestion:**
    * Place Planck `.fits` files (NPIPE/SMICA) in `data/raw/`.
    * Place SDSS/eBOSS catalogs in `data/raw/catalogs/`.

## License
MIT
