import click
from conformers.confgen import ConformerGenerator
from conformers.convert import mol_to_atoms
from conformers.density import ConformerAnalyzer


@click.command("density")
@click.argument("smiles", type=str)
@click.option("-o", "--optimize", is_flag=True, default=False)
@click.option("--fmax", type=float, default=0.005)
def density(smiles, optimize, fmax):
    cgen = ConformerGenerator(smiles)
    mol, _ = cgen.run()
    atoms = mol_to_atoms(mol)

    if optimize:
        from conformers.optimize import Optimizer

        opt = Optimizer(fmax=fmax)
        atoms = opt.run(atoms, logfile=None)

    analyzer = ConformerAnalyzer(atoms, smiles)
    density = analyzer.get_density()
    print(f"{smiles}, {density:.3f} g/cm3")
