# Conformer generation and analysis

This repository contains code to generate and analyze molecular conformers.
The code uses RDKit to create conformers from SMILES, and returns a list of ASE Atoms.
Optionally, the code allows optimizing them at the CCSD(T)-level with a pre-trained neural network potential (ANI-1ccx implemented in [TorchANI](https://github.com/aiqm/torchani)).

## Installation

To use this package, first clone the repository and install it under your virtual environment using pip:

```
git clone git@github.com:dskoda/conformers.git
cd conformers
pip install -e .
```

## Usage

The `conformers` package contains an API-like interface that can be imported in other code, and scripts that enable computing the conformers using the command line interface.
These scripts can be accessed using the `conformers` command that is installed with the repository.
For example, to generate an ASE Atoms object based on a single conformer:

```python
from conformers.confgen import ConformerGenerator
from conformers.convert import mol_to_atoms

cgen = ConformerGenerator(smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
mol, energies = cgen.run()
atoms = mol_to_atoms(mol)
```
