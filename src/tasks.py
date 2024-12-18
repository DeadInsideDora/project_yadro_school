import logging
from typing import List
from celery import Celery
from redis import Redis

import database
from models import Molecules
from substructure_search import substructure_search

celery = Celery(
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0'
)

redis_client = Redis(host='redis', port=6379, db=0)


logger = logging.getLogger(__name__)


@celery.task
def get_by_id_task(uuid: str) -> str:
    logger.info(f'get_by_id:{uuid}. create task')
    cache_key = f'get_by_id:{uuid}'
    with database.SessionLocal() as db:
        molecule = db.query(Molecules).filter_by(id=uuid).first()
        smile_name = str(molecule.smile_name) if molecule else 'None'
        redis_client.set(cache_key, smile_name)
        logger.info(f'delete task: task:{cache_key}')
        redis_client.delete(f'task:{cache_key}')
        return smile_name


@celery.task
def sub_search_task(pattern: str) -> List[str]:
    logger.info(f'sub_search_task:{pattern}. create task')
    cache_key = f'sub_search:{pattern}'
    with database.SessionLocal() as db:
        try:
            result = substructure_search(
                [molecule.smile_name for molecule
                 in db.query(Molecules).all()],
                pattern, logger)
        except ValueError as e:
            result = str(e)
    redis_client.set(cache_key, str(result))
    logger.info(f'delete task: task:{cache_key}')
    redis_client.delete(f'task:{cache_key}')
    return result
