# FLASK: Conformer generation and analysis

This repository contains code to generate and analyze molecular conformers.
The code uses RDKit to create conformers from SMILES, optimizes them at the CCSD(T)-level with a pre-trained neural network potential (ANI-1ccx implemented in [TorchANI](https://github.com/aiqm/torchani)), and analyzes the conformers according to heuristics.

## Installation

To use this package, first clone the repository and install it under your virtual environment using pip:

```
git clone ssh://git@czgitlab.llnl.gov:7999/dskoda/flask-conformers.git
cd flask_conformers
pip install -e .
```

## Usage

The `flask_conformers` package contains an API-like interface that can be imported in other code, and scripts that enable computing the conformers using the command line interface.
These scripts can be accessed using the `flask_conf` command that is installed with the repository.
For example, if you wanted to compute the estimated density of a molecule using the heuristic approach:

```
flask_conf density "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --optimize
```

The result will be something like:

```
CN1C=NC2=C1C(=O)N(C(=O)N2C)C, 1.390 g/cm3
```

If you want to process a list of SMILES, you can use the `process_list` command:

```
flask_conf process_list smiles.txt --optimize --output results.json
```

This command loops over all SMILES in the file (one per line) and saves the densities, optimized conformers, and energies, as predicted by the ANI-1ccx model.
If the `--optimize` flag is not used, the RDKit conformers are used instead (and no energies are saved).
This saves a bit of time as the structural optimization is not performed, but the quality of the final conformers may not be as good.

## Authors and acknowledgment

The code for this repository was written by Daniel Schwalbe-Koda (schwalbekoda1), who acknowledges data and helpful comments from Evan Antoniuk (antoniuk1).
The conformer generator code was copied from [mkite_conformer](https://github.com/mkite-group/mkite_conformer).

## License

This project is not open-sourced and is available for internal use of the FLASK team only.
This work was produced under the auspices of the U.S. Department of
Energy by Lawrence Livermore National Laboratory under Contract
DE-AC52-07NA27344.
