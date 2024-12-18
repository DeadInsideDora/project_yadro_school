import pytest
from src.substructure_search import substructure_search


def test_correct_data():
    res = substructure_search(
        ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
    assert res == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


def test_empty_result():
    res = substructure_search(["CN=C=O", "[OH-]", "[NH4+]", "c1ccccc1"], "CCO")
    assert not res


def test_invalid_pattern():
    molecules = ["CN=C=O", "CC(=O)Oc1ccccc1C(=O)O", "CCC(=O)O", "[OH-]"]
    pattern = "not smile sting"
    with pytest.raises(ValueError) as e:
        substructure_search(molecules, pattern)

    assert str(e.value) == f"'{pattern}' isn't smiles string"


def test_empty_molecules_list():
    res = substructure_search([], "CN=C=O")
    assert not res
