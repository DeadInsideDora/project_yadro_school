from sqlalchemy.orm import DeclarativeBase
from sqlalchemy import Column, String


class Base(DeclarativeBase):
    pass


class Molecules(Base):
    __tablename__ = "molecules"

    id = Column(String, primary_key=True, index=True)
    smile_name = Column(String)
