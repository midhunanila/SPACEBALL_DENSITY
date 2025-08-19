# SPACEBALL-v3.0


**Density Calculation of Protein Nanodroplets Using SPACEBALL**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

---

## Overview

**SPACEBALL v. 3.0** is a computational tool for **density calculation of protein nanodroplets**.  It analyzes PDB structures and computes the nanodroplet volume.  The previous version of [SPACEBALL v. 2.0](http://info.ifpan.edu.pl/~chwastyk/spaceball/) was primarily designed for cavity detection.  In this release, we introduce a new method to determine the **optimal probe radius**, enabling more accurate density estimation of protein nanodroplets. 

---

## Features

- Density calculations for protein nanodroplets  
- Supports multiple simulation file formats (XYZ, PDB, etc.)  
- High-performance computation for large datasets  

---

## Requirements

- **Python 3.x**  
  Install Python 3 from [python.org](https://www.python.org/) or via package manager:
  ```bash
  # Ubuntu/Debian
  sudo apt install python3 python3-pip

  # macOS
  brew install python

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
cd SPACEBALL-v3.0
