import logging

from typing import Optional
from uuid import uuid4
from os import getenv

from fastapi import FastAPI, HTTPException, UploadFile, File
from rdkit import Chem
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from schema import Base, Molecules
from substructure_search import substructure_search


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
log_handler = logging.FileHandler(f'{__name__}.log', mode='w')
log_formatter = logging.Formatter("%(asctime)s %(name)s %(funcName)s %(levelname)s %(message)s")
log_handler.setFormatter(log_formatter)
logger.addHandler(log_handler)

app = FastAPI()

engine = create_engine(f"postgresql://{getenv('POSTGRES_USER')}:"
                       f"{getenv('POSTGRES_PASSWORD')}@"
                       f"{getenv('POSTGRES_DB')}:5432/"
                       f"{getenv('POSTGRES_DB')}")

Base.metadata.create_all(bind=engine)


def validate_smile(smile_name: str) -> None:
    if Chem.MolFromSmiles(smile_name) is None:
        logger.error(f"{smile_name} isn't smiles string")
        raise HTTPException(status_code=400,
                            detail=f"{smile_name} isn't smiles string")


@app.get("/")
def get_server():
    logger.info('Request: GET /')
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add")
def add_molecule(smile_name: str) -> dict:
    logger.info(f'Request: POST /add?smile_name={smile_name}')
    validate_smile(smile_name)
    with Session(autoflush=False, bind=engine) as db:
        if db.query(Molecules).filter_by(smile_name=smile_name).first():
            logger.error(f"Molecule {smile_name} already added")
            raise HTTPException(status_code=409,
                                detail=f"Molecule {smile_name} already added")
        new_molecule = Molecules(id=str(uuid4()), smile_name=smile_name)
        db.add(new_molecule)
        db.commit()
        logger.info(f'New molecule {smile_name} was successfully added with id: {new_molecule.id}')

        return {"id": new_molecule.id}


@app.get("/get_by_id/{uuid}")
def get_molecule_by_id(uuid: str) -> str:
    logger.info(f'Request: GET //get_by_id/{uuid}')
    with Session(autoflush=False, bind=engine) as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            logger.error(f'Identifier {uuid} not found')
            raise HTTPException(status_code=404, detail="Identifier not found")
    logger.info(f'Molecule was successfully gotten by id= {uuid}: {molecule.smile_name}')
    return molecule.smile_name


@app.patch("/update_by_id/{uuid}")
def update_by_id(uuid: str, new_smile_name: str) -> None:
    logger.info(f'Request: PATCH /update_by_id/{uuid}')
    validate_smile(new_smile_name)
    with Session(autoflush=False, bind=engine) as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            logger.error(f'Identifier {uuid} not found')
            raise HTTPException(status_code=404, detail="Identifier not found")
        if db.query(Molecules).filter(
                Molecules.smile_name.in_([new_smile_name])).first():
            logger.error(f'Molecule {new_smile_name} already added, cant update')
            raise HTTPException(
                status_code=409,
                detail=f"Molecule {new_smile_name} already added, cant update"
            )
        molecule.smile_name = new_smile_name
        db.commit()
        logger.info(f'Molecule was successfully updated by id= {uuid}, new smile_name= {new_smile_name}')


@app.delete("/delete/{uuid}")
def delete_by_id(uuid: str) -> None:
    logger.info(f'Request: DELETE /delete/{uuid}')
    with Session(autoflush=False, bind=engine) as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            logger.error(f'Identifier {uuid} not found')
            raise HTTPException(status_code=404, detail="Identifier not found")
        db.delete(molecule)
        db.commit()
        logger.info(f'Molecule was successfully deleted by id={uuid}')


@app.get("/get_all")
def get_molecules_list(limit: Optional[int] = None):
    logger.info('Request: GET /get_all' + (f'?limit={limit}' if limit else ''))
    if limit is not None and limit <= 0:
        logger.error('limit value must be positive integer')
        raise HTTPException(status_code=400,
                            detail='limit value must be positive integer')
    with Session(autoflush=False, bind=engine) as db:
        query = db.query(Molecules)
        if limit:
            query = query.limit(limit)
        return {
            "smile_name": [
                molecule.smile_name for molecule in query.all()]}


@app.get("/sub_search")
def get_search_result(pattern: str):
    logger.info('Request: GET /sub_search?pattern={pattern}')
    with Session(autoflush=False, bind=engine) as db:
        return {"sub_search_result": substructure_search(
            [molecule.smile_name for molecule in db.query(Molecules).all()],
            pattern, logger)
        }


@app.post("/add_from_file")
def upload_file(file: UploadFile = File(...)):
    logger.info(f'Request: POST /add_from_file -F file={file}')
    smiles_to_add = file.file.read().decode().split('\n')
    with Session(autoflush=False, bind=engine) as db:
        added_molecules = []
        for line in smiles_to_add:
            validate_smile(line)
            if db.query(Molecules).filter_by(smile_name=line).first():
                logger.info(f'Molecule {line} already added')
                raise HTTPException(
                    status_code=409,
                    detail=f'Molecule {line} already added')
            added_molecules.append(Molecules(id=str(uuid4()), smile_name=line))
        db.add_all(added_molecules)
        db.commit()
        logger.info(f'All molecules from file={file} were successfully added')
        return {"added": [{"smile_name": molecule.smile_name,
                           "id": molecule.id} for molecule in added_molecules]}
