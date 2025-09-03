# Density Calculation of Protein Nanodroplets Using SPACEBALL


[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

---

## Overview

Here is a computational tool designed for the density calculation of protein nanodroplets by determining their volume. It analyzes PDB structures to compute the nanodroplet volume, which can then be used to estimate density accurately. The previous version [SPACEBALL v. 2.0](http://info.ifpan.edu.pl/~chwastyk/spaceball/) was primarily focused on calculating the volume of cavities within proteins.  In this release, we introduce a **new method to determine the optimal probe radius**, enabling more precise density estimation of protein nanodroplets. By using this optimized probe radius, researchers can obtain **more reliable volume measurements**,  supporting studies in **biophysics, phase separation, and protein nanostructure analysis**.

## Repository structure

```bash
---
SPACEBALL_DENSITY/
├── LICENSE                 # License file (MIT)
├── README.md               # Documentation and usage guide
│
│
├── scripts/                # Python helper scripts
│   ├── probe.py        # Calculates optimal probe radius & prepares input file for SPACEBALL calculation
│
└── simulation_clusters               # Example PDB input file
            
```

## Requirements

- **Python 3.x**  
  Install Python 3 from [python.org](https://www.python.org/) or via package manager:
  ```bash
  # Ubuntu/Debian
  sudo apt install python3 python3-pip

  # macOS
  brew install python

  # Install required packages using pip
  pip install math
  pip install argparse

       

## Installation

Clone the repository:

```bash
git clone https://github.com/midhunanila/SPACEBALL_DENSITY.git
```


## Calculation

### Probe Radius

The Python scripts in this repository calculate the optimal probe radius for a protein nanodroplet. The input file should be in standard PDB format.

---

## Steps to Run

### 1. calculation of optimal probe radius with `probe.py`

1. Navigate to the `scripts` folder:

```bash
cd scripts
```

2. Run the script:

```bash
python probe.py <pdb_file> --chain_size <chainsize> --surface_fractions <surface_fraction>
```

- `--chain_size`: Set according to your protein  
- `--surface_fractions`: Fraction of the surface (default: 0.5)  

This script gives the corresponding value of optimal probe radius in **Ångströms**

Use this value on the Spaceball [website](http://info.ifpan.edu.pl/~chwastyk/spaceball/)

### 2. Input parameter description

The input parameters for the SPACEBALL calculation. Each line corresponds to a key parameter:

- **`pdb_structure: <pdb_file>`**  
  The input structure file in PDB format contains the 3D coordinates of all atoms in the system.

- **`wall_probe_radius_[A]: probe radius `**  
  Radius of the wall probe in **Ångströms**. Defines the outer boundary for SPACEBALL calculations.

- **`grid_X_[A]: 3`**  
  Grid spacing along the X-axis in **Ångströms**. Determines the resolution of the computational grid. The default value of the lattice constant, a = 3 **Ångströms**, corresponds roughly to the water molecule diameter

- **`grid_Y_[A]: 3`**  
  Grid spacing along the Y-axis in **Ångströms**.

- **`grid_Z_[A]: 3`**  
  Grid spacing along the Z-axis in **Ångströms**.


- **`number_of_rotations: 1`**  
  Number of rotational orientations applied to the probe during the calculation.


---

### 4. Density calculation

From the volume obtained from the given platform, we can compute the density of the protein nanodroplet

---

### 5. Example workflow

```bash
# Generate probe radius
cd scripts
python probe.py droplet.pdb --chain_size 140 --surface_fractions 0.5
```

---

## Notes

- Ensure Python and its necessary packages are installed.
-  The input PDB file must be in standard PDB format.
- Adjust `chain_size` according to your protein.  
- Verify the parameters for the SPACEBALL calculation are correct.  
- Paths are relative; adjust if your folder structure differs.

## Developers

| Name                     | Role / Contribution                  |
|--------------------------|---------------------------------------|
| Midhun Mohan Anila       | Python scripts and algorithm design   |
| Bartosz Rozycki          | algorithm design, corrections and validation |
| Michał Wojciechowski     | corrections and validation   |
| Mateusz Chwastyk         | developer of the SPACEBALL method and website, algorithm design, corrections and validation|
| Jan Malinowski           | developer of the SPACEBALL website |

## Publication (Under Review)

Our methodology and the SPACEBALL tool are described in a recent manuscript currently **Accepted**:

- **Title:** Theoretical Methods for Assessing the Density of Protein Nanodroplets 
- **Authors:** Midhun Mohan Anila, Michał Wojciechowski, Mateusz Chwastyk, Bartosz Rozycki  
- **Status:** Accepted 


