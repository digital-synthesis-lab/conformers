from typing import Dict

import numpy as np
from ase import Atoms
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors
from scipy.spatial.distance import cdist


N_AVOGADRO = 6.0221408e23  # mol^(-1)
CM3 = 1e-24  # Å^3
DEFAULT_RADII = {"H": 1.496, "C": 2.013, "O": 1.780, "N": 1.858}


class ConformerAnalyzer:
    def __init__(
        self,
        atoms: Atoms,
        smiles: str,
        radii: Dict[str, float] = DEFAULT_RADII,
        grid_size: float = 0.1,
        pad: float = None,
    ):
        self.atoms = atoms
        self.smiles = smiles
        self.radii = radii
        self.grid_size = grid_size

        if pad is None:
            pad = max(radii.values())

        self.pad = pad

    def get_volume(self) -> float:
        """Computes the volume (in Å^3) of a molecule specified by `atoms`.
        The volume is computed using a voxel approach that uses cutoffs specified
            by `radii` and a grid size (in Å).

        Arguments:
            atoms (Atoms): ASE Atoms containing the information about the molecule
            radii (dict): radii of the atoms to consider the region as part of the
                molecule or not.
            grid_size (float): length of the voxel.
            pad (float): padding around the molecule. If not specified, will
                consider the maximum radius in `radii` as the padding.

        Returns:
            volume (float): the volume of the molecule in Å^3
        """
        positions = self.atoms.get_positions()
        pos_radii = np.array(
            [self.radii.get(sym, 1.0) for sym in self.atoms.get_chemical_symbols()]
        )

        # making the voxels
        x_min, y_min, z_min = positions.min(axis=0) - self.grid_size - self.pad
        x_max, y_max, z_max = positions.max(axis=0) + self.grid_size + self.pad
        x, y, z = np.mgrid[
            x_min : x_max : self.grid_size,
            y_min : y_max : self.grid_size,
            z_min : z_max : self.grid_size,
        ]
        voxels = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
        grid_vol = self.grid_size**3

        # computing the voxel using vector operations
        dm = cdist(voxels, positions)
        diff = dm - pos_radii.reshape(1, -1)
        isin = np.min(diff, axis=1) <= 0
        volume = isin.sum() * grid_vol

        return volume

    def get_mass(self):
        """Compute the mass (in g/mol) of a given molecule
        using its SMILES string"""
        mol = Chem.MolFromSmiles(self.smiles)
        mol = Chem.AddHs(mol, addCoords=True)
        mass = Descriptors.ExactMolWt(mol)
        return mass

    def get_density(self):
        """Computes the density (in g/cm^3) of a molecule specified by `atoms`
            and a SMILES string.
        The volume is computed using a voxel approach that uses cutoffs specified
            by `radii` and a grid size (in Å).

        Arguments:
            atoms (Atoms): ASE Atoms containing the information about the molecule
            radii (dict): radii of the atoms to consider the region as part of the
                molecule or not.
            grid_size (float): length of the voxel.
            pad (float): padding around the molecule. If not specified, will
                consider the maximum radius in `radii` as the padding.

        Returns:
            volume (float): the volume of the molecule in Å^3
        """
        mass = self.get_mass()
        vol = self.get_volume() * CM3 * N_AVOGADRO

        return mass / vol
