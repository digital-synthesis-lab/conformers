from ase import Atoms
from rdkit.Chem import AllChem as Chem


def canonical_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def atoms_to_mol(atoms: Atoms, smiles: str = None):
    if smiles is None:
        mol = Chem.RWMol()
    else:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol, addCoords=True)

    for z in atoms.get_atomic_numbers():
        mol.AddAtom(Chem.Atom(z))

    conformer = Chem.Conformer(len(atoms))
    for i, pos in enumerate(atoms.positions):
        conformer.SetAtomPosition(i, pos)

    mol.AddConformer(conformer)
    return mol


def mol_to_atoms(mol: Chem.Mol, conformer_id: int = -1):
    symbols = [at.GetSymbol() for at in mol.GetAtoms()]
    conf = mol.GetConformer(conformer_id)
    positions = conf.GetPositions().tolist()
    return Atoms(symbols=symbols, positions=positions)
