from os import getenv
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker


DATABASE_URL = (f"postgresql://{getenv('POSTGRES_USER')}:"
                f"{getenv('POSTGRES_PASSWORD')}@"
                f"db:5432/{getenv('POSTGRES_DB')}")


engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()
