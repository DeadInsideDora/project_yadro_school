from sqlalchemy import Column, String
from database import Base


class Molecules(Base):
    __tablename__ = "molecules"

    id = Column(String, primary_key=True, index=True)
    smile_name = Column(String)
