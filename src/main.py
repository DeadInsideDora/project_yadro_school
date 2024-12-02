from os import getenv

from fastapi import FastAPI, HTTPException, UploadFile, File
from rdkit import Chem
from uuid import uuid4
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from database_meta import Base, Molecules
from substructure_search import substructure_search

app = FastAPI()

engine = create_engine(f"postgresql://{getenv('POSTGRES_USER')}:"
                       f"{getenv('POSTGRES_PASSWORD')}@"
                       f"{getenv('POSTGRES_DB')}:5432/"
                       f"{getenv('POSTGRES_DB')}")

Base.metadata.create_all(bind=engine)


def validate_smile(smile_name: str) -> None:
    if Chem.MolFromSmiles(smile_name) is None:
        raise HTTPException(status_code=400,
                            detail=f"{smile_name} isn't smiles string")


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add")
def add_molecule(smile_name: str) -> dict:
    validate_smile(smile_name)
    with Session(autoflush=False, bind=engine) as db:
        if db.query(Molecules).filter_by(smile_name=smile_name).first():
            raise HTTPException(status_code=409,
                                detail=f"Molecule {smile_name} already added")
        new_molecule = Molecules(id=str(uuid4()), smile_name=smile_name)
        db.add(new_molecule)
        db.commit()

        return {"id": new_molecule.id}


@app.get("/get_by_id/{uuid}")
def get_molecule_by_id(uuid: str) -> str:
    with Session(autoflush=False, bind=engine) as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            raise HTTPException(status_code=404, detail="Identifier not found")
    return molecule.smile_name


@app.patch("/update_by_id/{uuid}")
def update_by_id(uuid: str, new_smile_name: str) -> None:
    validate_smile(new_smile_name)
    with Session(autoflush=False, bind=engine) as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            raise HTTPException(status_code=404, detail="Identifier not found")
        if db.query(Molecules).filter(
                Molecules.smile_name.in_([new_smile_name])).first():
            raise HTTPException(
                status_code=409,
                detail=f"Molecule {new_smile_name} already added, cant update"
            )
        molecule.smile_name = new_smile_name
        db.commit()


@app.delete("/delete/{uuid}")
def delete_by_id(uuid: str) -> None:
    with Session(autoflush=False, bind=engine) as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            raise HTTPException(status_code=404, detail="Identifier not found")
        db.delete(molecule)
        db.commit()


@app.get("/get_all")
def get_molecules_list():
    with Session(autoflush=False, bind=engine) as db:
        return {
            "smile_name": [
                molecule.smile_name for molecule in db.query(Molecules).all()]}


@app.get("/sub_search")
def get_search_result(pattern: str):
    with Session(autoflush=False, bind=engine) as db:
        return {"sub_search_result": substructure_search(
            [molecule.smile_name for molecule in db.query(Molecules).all()],
            pattern)
        }


@app.post("/add_from_file")
def upload_file(file: UploadFile = File(...)):
    smiles_to_add = file.file.read().decode().split('\n')
    with Session(autoflush=False, bind=engine) as db:
        added_molecules = []
        for line in smiles_to_add:
            validate_smile(line)
            if db.query(Molecules).filter_by(smile_name=line).first():
                raise HTTPException(
                    status_code=409,
                    detail=f"Molecule {line} already added")
            added_molecules.append(Molecules(id=str(uuid4()), smile_name=line))
        db.add_all(added_molecules)
        db.commit()
        return {"added": [{"smile_name": molecule.smile_name,
                           "id": molecule.id} for molecule in added_molecules]}
