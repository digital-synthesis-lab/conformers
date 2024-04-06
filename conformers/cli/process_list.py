import json
import sys
import os

import click
import tqdm
from conformers.confgen import ConformerGenerator
from conformers.convert import mol_to_atoms
from conformers.density import ConformerAnalyzer
from conformers.json import NumpyEncoder


@click.command("process_list")
@click.argument("smiles_list", type=str)
@click.option("--optimize", is_flag=True, default=False)
@click.option("--fmax", type=float, default=0.005)
@click.option("-o", "--output", type=str, default="density.json")
def process_list(smiles_list, optimize, fmax, output):
    assert os.path.exists(smiles_list), f"Path {smiles_list} does not exist."

    if os.path.exists(output):
        print(f"Output {output} exists. Exiting...")
        sys.exit()

    with open(smiles_list, "r") as f:
        all_smiles = [s.strip() for s in f.readlines()]

    if optimize:
        from conformers.optimize import Optimizer

        opt = Optimizer(fmax=fmax)

    results = {}
    for smiles in tqdm.tqdm(all_smiles):
        cgen = ConformerGenerator(smiles)
        mol, _ = cgen.run()
        atoms = mol_to_atoms(mol)

        extra = {}
        if optimize:
            atoms = opt.run(atoms, logfile=None)
            extra["energy"] = atoms.get_potential_energy()

        analyzer = ConformerAnalyzer(atoms, smiles)
        density = analyzer.get_density()

        results[smiles] = {"density": density, "conformer": atoms.todict(), **extra}

    with open(output, "w") as f:
        json.dump(results, f, cls=NumpyEncoder)
