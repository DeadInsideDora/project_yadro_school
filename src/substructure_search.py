from typing import List
from rdkit import Chem


def substructure_search(
        molecules: List[str],
        pattern_molecule: str) -> List[str]:
    pattern = Chem.MolFromSmiles(pattern_molecule)
    if not pattern:
        raise ValueError(f"'{pattern_molecule}' isn't smiles string")
    return [molecule for molecule in molecules
            if Chem.MolFromSmiles(molecule).HasSubstructMatch(pattern)]
