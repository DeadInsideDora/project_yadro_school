import logging
from typing import Optional
from os import getenv
from celery.result import AsyncResult
from fastapi import FastAPI, HTTPException, UploadFile, File
from redis import Redis
from rdkit import Chem
from uuid import uuid4

import database
from tasks import get_by_id_task, sub_search_task, celery
from models import Base, Molecules

logging.basicConfig(
    filename=f'{__name__}.log',
    level=logging.INFO,
    format="%(asctime)s %(name)s %(funcName)s %(levelname)s %(message)s")

logger = logging.getLogger(__name__)

Base.metadata.create_all(bind=database.engine)

app = FastAPI()

redis_client = Redis(host='redis', port=6379)


def del_sub_search():
    for key in redis_client.keys(pattern='sub_search:*'):
        redis_client.delete(key.decode('utf-8'))


@app.get("/get_by_id/{uuid}")
def get_molecule_by_id(uuid: str) -> dict:
    logger.info(f'Request: GET //get_by_id/{uuid}')
    cache_key = f"get_by_id:{uuid}"
    cached_result = redis_client.get(cache_key)

    if cached_result:
        return {"result": cached_result.decode('utf-8'), "status": "completed"}
    else:
        task_id = redis_client.get(f"task:{cache_key}")
        if not task_id:
            task = get_by_id_task.delay(uuid)
            task_id = task.id
            redis_client.set(f"task:{cache_key}", task_id)
            return {"task_id": task_id, "status": "processing",
                    "message": "Calculation started. Please check back later."}
        else:
            return {
                "task_id": task_id,
                "status": "processing",
                "message": "Calculation is in progress. "
                           "Please check back later."}


def validate_smile(smile_name: str) -> None:
    if Chem.MolFromSmiles(smile_name) is None:
        logger.error(f"{smile_name} isn't smiles string")
        raise HTTPException(status_code=400,
                            detail=f"{smile_name} isn't smiles string")


def get_molecule(uuid: str):
    logger.info(f'get_by_id({uuid})')
    cache_key = f'get_by_id:{uuid}'
    result = redis_client.get(cache_key)
    if result:
        logger.info(
            f'Molecule was successfully gotten by id={uuid} from cache')
        return result.decode()

    with database.SessionLocal() as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            logger.error(f'Identifier {uuid} not found')
            raise HTTPException(status_code=404, detail="Identifier not found")
    logger.info(
        f'Molecule was successfully gotten by id={uuid}:'
        f'{molecule.smile_name}')
    redis_client.set(cache_key, molecule.smile_name)
    return molecule.smile_name


@app.get("/")
def get_server():
    logger.info('Request: GET /')
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/add")
def add_molecule(smile_name: str) -> dict:
    logger.info(f'Request: POST /add?smile_name={smile_name}')
    validate_smile(smile_name)
    with database.SessionLocal() as db:
        if db.query(Molecules).filter_by(smile_name=smile_name).first():
            logger.error(f"Molecule {smile_name} already added")
            raise HTTPException(status_code=409,
                                detail=f"Molecule {smile_name} already added")
        new_molecule = Molecules(id=str(uuid4()), smile_name=smile_name)
        db.add(new_molecule)
        db.commit()
        del_sub_search()
        logger.info(
            f'New molecule {smile_name} was successfully added with id:'
            f'{new_molecule.id}')

        return {"id": new_molecule.id}


@app.patch("/update_by_id/{uuid}")
def update_by_id(uuid: str, new_smile_name: str) -> None:
    logger.info(f'Request: PATCH /update_by_id/{uuid}')
    validate_smile(new_smile_name)
    with database.SessionLocal() as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            logger.error(f'Identifier {uuid} not found')
            raise HTTPException(status_code=404, detail="Identifier not found")
        if db.query(Molecules).filter(
                Molecules.smile_name.in_([new_smile_name])).first():
            logger.error(
                f'Molecule {new_smile_name} already added, cant update')
            raise HTTPException(
                status_code=409,
                detail=f"Molecule {new_smile_name} already added, cant update"
            )
        molecule.smile_name = new_smile_name
        db.commit()
        redis_client.delete(f'get_by_id:{uuid}')
        del_sub_search()
        logger.info(
            f'Molecule was successfully updated by id= {uuid},'
            f'new smile_name= {new_smile_name}')


@app.delete("/delete/{uuid}")
def delete_by_id(uuid: str) -> None:
    logger.info(f'Request: DELETE /delete/{uuid}')
    with database.SessionLocal() as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        if not molecule:
            logger.error(f'Identifier {uuid} not found')
            raise HTTPException(status_code=404, detail="Identifier not found")
        db.delete(molecule)
        db.commit()
        redis_client.delete(f'get_by_id:{uuid}')
        del_sub_search()
        logger.info(f'Molecule was successfully deleted by id={uuid}')


@app.get("/get_all")
def get_molecules_list(limit: Optional[int] = None):
    logger.info('Request: GET /get_all' + (f'?limit={limit}' if limit else ''))
    if limit is not None and limit <= 0:
        logger.error('limit value must be positive integer')
        raise HTTPException(status_code=400,
                            detail='limit value must be positive integer')
    with database.SessionLocal() as db:
        query = db.query(Molecules)
        if limit:
            query = query.limit(limit)
        return {
            "smile_name": [
                molecule.smile_name for molecule in query.all()]}


@app.get("/sub_search")
def get_search_result(pattern: str):
    logger.info(f'Request: GET /sub_search?pattern={pattern}')
    cache_key = f"sub_search:{pattern}"
    cached_result = redis_client.get(cache_key)

    if cached_result:
        return {"result": cached_result.decode('utf-8'), "status": "completed"}
    else:
        task_id = redis_client.get(f"task:{cache_key}")
        if not task_id:
            task = sub_search_task.delay(pattern)
            task_id = task.id
            redis_client.set(f"task:{cache_key}", task_id)
            return {"task_id": task_id, "status": "processing",
                    "message": "Calculation started. Please check back later."}
        else:
            return {"task_id": task_id, "status": "processing",
                    "message": "Calculation started. Please check back later."}


@app.get("/tasks/{task_id}")
def get_task_result(task_id: str):
    logger.info(f'Request: GET /tasks/{task_id}')
    task_result = AsyncResult(task_id, app=celery)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {
            "task_id": task_id,
            "status": "Task completed",
            "result": task_result.result}
    else:
        return {"task_id": task_id, "status": task_result.state}


@app.post("/add_from_file")
def upload_file(file: UploadFile = File(...)):
    logger.info(f'Request: POST /add_from_file -F file={file}')
    smiles_to_add = file.file.read().decode().split('\n')
    with database.SessionLocal() as db:
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
        del_sub_search()
        logger.info(f'All molecules from file={file} were successfully added')
        return {"added": [{"smile_name": molecule.smile_name,
                           "id": molecule.id} for molecule in added_molecules]}
