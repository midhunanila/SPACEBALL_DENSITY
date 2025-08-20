# SPACEBALL v.3.0


**Density Calculation of Protein Nanodroplets Using SPACEBALL**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

---

## Overview

**SPACEBALL v. 3.0** is a computational tool for **density calculation of protein nanodroplets**.  It analyzes PDB structures and computes the nanodroplet volume.  The previous version of [SPACEBALL v. 2.0](http://info.ifpan.edu.pl/~chwastyk/spaceball/) was primarily designed for cavity detection.  In this release, we introduce a new method to determine the **optimal probe radius**, enabling more accurate density estimation of protein nanodroplets. 

## Repository structure

```bash
---
SPACEBALL-v3.0/
├── LICENSE                 # License file (MIT)
├── README.md               # Documentation and usage guide
│
├── spaceball/              # Fortran source code and executable
│   ├── Makefile            # Run `make` here to compile
│   ├── *.f / *.f90         # Fortran source files
│   └── prot_hole_paral3.xg # Compiled executable (created after `make`)
│
├── scripts/                # Python helper scripts
│   ├── probe.py        # Calculates optimal probe radius & prepares input file for SPACEBALL calculation
│   ├── density.py     # Computes density from PDB + SPACEBALL output file
│   └── (other utilities if added)
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


- **gfortran (Fortran compiler)**
  ```bash
  # Ubuntu/Debian
  sudo apt install gfortran
  
  # macOS (Homebrew)
  brew install gcc  # includes gfortran


## Installation

Clone the repository:

```bash
git clone https://github.com/midhunanila/SPACEBALL-v3.0.git
cd SPACEBALL-v3.0/spaceball
```

Compile the Spaceball executable:

```bash
make
```

This will create the executable file named `prot_hole_paral3.xg` for SPACEBALL calculations.

---

## Calculation

### Probe Radius

The Python scripts in this repository calculate the optimal probe radius for a protein nanodroplet. The input file should be in standard PDB format.

---

## Steps to Run

### 1. Generate Probe Radius with `probe.py`

1. Navigate to the `scripts` folder:

```bash
cd ../scripts
```

2. Run the script:

```bash
python probe.py <pdb_file> --chain_size <chainsize> --surface_fractions <surface_fraction>
```

- `--chain_size`: Set according to your protein  
- `--surface_fractions`: Fraction of the surface (default: 0.5)  

This script sets the probe radius and generates the input file for SPACEBALL.

### 2. Input File Description

The input file specifies parameters for running SPACEBALL. Each line corresponds to a key parameter:

- **`pdb_structure: <pdb_file>`**  
  The input structure file in PDB format containing the 3D coordinates of all atoms in the system.

- **`wall_probe_radius_[A]: proberadius `**  
  Radius of the wall probe in **Ångströms**. Defines the outer boundary for SPACEBALL calculations.

- **`water_probe_radius_[A]: 1.42`**  
  Radius of the water probe in **Ångströms**. Represents the effective size of a water molecule when probing the structure. The default value is 1.42 **Ångströms**

- **`grid_X_[A]: 3`**  
  Grid spacing along the X-axis in **Ångströms**. Determines the resolution of the computational grid. The default value of the lattice constant, a = 3 **Ångströms**, corresponds roughly to the water molecule diameter

- **`grid_Y_[A]: 3`**  
  Grid spacing along the Y-axis in **Ångströms**.

- **`grid_Z_[A]: 3`**  
  Grid spacing along the Z-axis in **Ångströms**.

- **`number_of_clusters_written_to_the_output: 1`**  
  Number of clusters to save in the output file. Only the largest or first cluster will be written.

- **`machine_[1_PC/2_CLUSTER]: 4`**  
  Machine configuration or identifier. Could indicate parallelization (1 = single PC, 2 = cluster). Here, `4` may specify number of cores or threads.

- **`number_of_rotations: 1`**  
  Number of rotational orientations applied to the probe during the calculation.

- **`EOF`**  
  End-of-file marker indicating the input file ends here.


---

### 3. Run SPACEBALL Calculation

1. Navigate to the `spaceball` folder if not already there:

```bash
cd ../spaceball
```

2. Run the executable:

```bash
./prot_hole_paral3.xg inputfile
```

- Use the probe radius from `probe.py` in the input file.  
- Output will be saved as:

```bash
pdbfilename.out
```

*(replace `pdbfilename` with your PDB file name).*

---
### 4. density calculation

1. Navigate to the `spaceball` folder if not already there:

```bash
cd ../scripts
```

2. Run the script:
   
```bash
python density.py <pdb_file> <out_file>
```
The result gives the density of the droplet  in units of **res/nm³**

### 5. Example Workflow

```bash
# Generate probe radius
cd scripts
python spaceball.py droplet.pdb --chain_size 140 --surface_fractions 0.5

# Compile Spaceball
cd ../spaceball
make

# Run SPACEBALL calculation
./prot_hole_paral3.xg inputfile

# Check output
ls -l pdbfilename.out
# Density calculation
python density_calc.py droplet.pdb pdbfilename.out
```

---

## Notes

- Ensure Python and a Fortran compiler are installed.  
- Adjust `chain_size` according to your protein.  
- Verify `inputfile` contains the correct probe radius values.  
- Paths are relative; adjust if your folder structure differs.
