from ase import Atoms
from ase.optimize import BFGS
from torchani.models import ANI1ccx


class Optimizer:
    def __init__(self, fmax: float = 0.005):
        """Creates an optimizer using TorchANI (model ANI-1ccx)"""
        self.fmax = fmax

    def get_model(self, **kwargs):
        return ANI1ccx(**kwargs)

    def get_calculator(self):
        model = self.get_model()
        return model.ase()

    def run(self, atoms: Atoms, **kwargs):
        """Optimizes an Atoms object using ANI-1ccx"""
        atoms = atoms.copy()
        calc = self.get_calculator()
        atoms.set_calculator(calc)
        opt = BFGS(atoms, **kwargs)
        opt.run(fmax=self.fmax)

        return atoms
