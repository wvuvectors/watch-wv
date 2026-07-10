# app/models/county.py

from sqlalchemy import Column, String, Float
from app.database import Base

class County(Base):
    __tablename__ = "county"
    
    county_id = Column(String(10), primary_key=True)
    county_labcode = Column(String(5), nullable=True)
    county_fips = Column(String(10), nullable=True)
    county_name = Column(String(20), nullable=True)
    county_population = Column(Float, nullable=True)