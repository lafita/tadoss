# TADOSS - TAndem DOmain Swap Stability prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

TADOSS estimates the stability of tandem domain swap conformations using an alchemical approximation based on coarse-grained (Go-like) simulation models from the three-dimensional structure of a protein domain.
The stability is defined as the relative free energy difference (ΔΔG) between the native and swap conformations of a pair of tandem identical domains.

## Installation

### Requirements

- Python 2.7: https://www.python.org/download/releases/2.7
- BioPython: http://biopython.org/wiki/Download
- One of the following two:
  - GROMACS 5.0 or higher: http://manual.gromacs.org/documentation
  - Reduce: http://kinemage.biochem.duke.edu/software/reduce.php

Additional
- R version 3.4 or higher https://www.r-project.org, including:
  - ggplot2
  - optparse
- Molecular 3D viewer:
  - PyMOL: https://pymol.org/2
  - Chimera: https://www.cgl.ucsf.edu/chimera/download.html

### Setup

1. Make sure you have installed the required software and packages (see above).
2. Clone this repository to your desired location.
```bash
git clone git@github.com:lafita/tadoss.git ~/tadoss
```
3. Optional: make tadoss accessible in your `PATH` to access it from anywhere.
```bash
ln -s /usr/bin/tadoss ~/tadoss
```

## Usage

The method consist in four steps starting from a `PDB` file of a protein domain structure (`e1shgA1.pdb` as example).

### Bundle script

In order to simplify the usage of the method, a `Bash` script that bundles the four steps steps is provided with a simple interface to the user.

```bash
tadoss -d ubq -f e1shgA1.pdb
```

### Individual steps

Each of the steps can also be run separately to allow more flexibility and control to users.

1. Download and prepare the domain structure (ECOD domain `e1shgA1` as example):
```bash
# Under the example folder in this repository
cd example
mkdir tmp # temporary intermediate files will be stored here
```
```bash
# Select only protein residues from the structure
python2.7 ~/tadoss/trim_nonprotein.py e1shgA1.pdb tmp/sh3_protein.pdb
```

2. Add missing hydrogens to the structure:
```bash
# Using GROMACS
gmx pdb2gmx -f tmp/sh3_protein.pdb -o tmp/sh3_hadded.pdb -ignh <<< $'14\n3'
```
```bash
# Using Reduce (use the path to the "reduce" script)
reduce tmp/sh3_protein.pdb > tmp/sh3_hadeed.pdb
```

3. Generate native contact energy (Go model) from the structure:
```
python2.7 ~/go_builder/go_builder.py -k sh3 --gmx tmp/sh3_hadded.pdb
```

The generated native contact energies are listed in file `go_sh3/go_sh3_gomodel_golist.dat`.
Other files are not needed for domain swap stability estimation.

The first two columns of the file contain the pair of residues forming the contact (i,j), the third column contains the distance between them in Amstrongs and the last column is the absolute GO energy of the interaction in kcal/mol.

4. Alchemical domain swap stability estimation from the GO model contact energies:
```
python ~/tadoss.py sh3 tmp/sh3_hadded.pdb
```

### Results

Results are split into two files: [sh3-dG_cut.tsv](example/sh3-dG_cut.tsv) and [sh3-dG_join.tsv](example/sh3-dG_join.tsv).

The final free energy difference (ΔΔG) can be obtained by summing up the contributions in each file:
```
ΔΔG = ΔGjoin + ΔGcut
```

More information about the output can be found here.

### Visualization

In order to visualize the mapping of the profile ΔGc onto the native domain structure, a new PDB file `ubq-dG_cut.pdb` is also generated with the ΔGc values in the B-factor column. Using a molecular viewer like Pymol, residues can be colored according to their ΔGc values.

In addition, two R scripts are provided to visualize the results.

## Publications

- Manuscript in preparation.
- Original description by Tian and Best (2016): https://doi.org/10.1371/journal.pcbi.1004933
