from typing import List
from rdkit import Chem
from logging import Logger


def substructure_search(
        molecules: List[str],
        pattern_molecule: str,
        logger: Logger = None) -> List[str]:
    pattern = Chem.MolFromSmiles(pattern_molecule)
    if not pattern:
        if logger:
            logger.info(f"'{pattern_molecule}' isn't smiles string")
        raise ValueError(f"'{pattern_molecule}' isn't smiles string")
    return [molecule for molecule in molecules
            if Chem.MolFromSmiles(molecule).HasSubstructMatch(pattern)]
