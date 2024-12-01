from os import getenv

from fastapi import FastAPI, HTTPException, UploadFile, File
from rdkit import Chem
from uuid import uuid4
from typing import List


smile_molecules = dict()

app = FastAPI()


def substructure_search(molecules: List[str], pattern_molecule: str) -> List[str]:
    pattern = Chem.MolFromSmiles(pattern_molecule)
    if not pattern:
        raise ValueError(f"'{pattern_molecule}' isn't smiles string")
    return [molecule for molecule in molecules
            if Chem.MolFromSmiles(molecule).HasSubstructMatch(pattern)]


def add_base(smile_name: str):
    if Chem.MolFromSmiles(smile_name) is None:
        raise HTTPException(status_code=400, detail=f"{smile_name} isn't smiles string")
    if smile_name in smile_molecules.values():
        raise HTTPException(status_code=409, detail=f"Molecule {smile_name} already added")
    unique_id = str(uuid4())
    smile_molecules[unique_id] = smile_name
    return unique_id


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add")
def add_molecule(smile_name: str) -> dict:
    return {"id": add_base(smile_name)}


@app.get("/get_by_id/{uuid}")
def get_molecule_by_id(uuid: str) -> str:
    if uuid not in smile_molecules.keys():
        raise HTTPException(status_code=404, detail="Identifier not found")
    return smile_molecules[uuid]


@app.patch("/update_by_id/{uuid}")
def update_by_id(uuid: str, new_smile_name: str) -> None:
    if uuid not in smile_molecules.keys():
        raise HTTPException(status_code=404, detail="Identifier not found")
    smile_molecules[uuid] = new_smile_name


@app.delete("/delete/{uuid}")
def delete_by_id(uuid: str) -> None:
    if uuid not in smile_molecules.keys():
        raise HTTPException(status_code=404, detail="Identifier not found")
    del smile_molecules[uuid]


@app.get("/get_all")
def get_molecules_list():
    return {"smiles_name": list(smile_molecules.values())}


@app.get("/sub_search")
def get_search_result(pattern: str):
    return {"sub_search_result": substructure_search(list(smile_molecules.values()), pattern)}


@app.post("/add_from_file")
def upload_file(file: UploadFile = File(...)):
    smiles_to_add = file.file.read().decode().split('\n')
    if all(smile not in smile_molecules.values() for smile in smiles_to_add):
        added_molecules = []
        for line in smiles_to_add:
            added_molecules.append({"smile_name": line, "id": add_base(line)})
        return {"added": added_molecules}
    raise HTTPException(status_code=409, detail=f"File contains one or more already added molecules")


