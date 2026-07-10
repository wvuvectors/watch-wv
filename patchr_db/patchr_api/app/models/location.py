# app/models/location.py

from sqlalchemy import Column, String, Text, Float, DECIMAL
from app.database import Base

class Location(Base):
    __tablename__ = "location"
    
    location_id = Column(String(100), primary_key=True)
    sample_code_prefix = Column(String(20), nullable=True)
    location_primary_lab = Column(String(20), nullable=True)
    location_status = Column(String(50), nullable=True)
    location_common_name = Column(String(255), nullable=True)
    location_category = Column(String(100), nullable=True)
    location_group = Column(String(100), nullable=True)
    location_lng = Column(DECIMAL(10, 6), nullable=True)
    location_lat = Column(DECIMAL(10, 6), nullable=True)
    location_primary_wwtp_id = Column(String(100), nullable=True)
    location_counties_served = Column(String(255), nullable=True)
    location_population_served = Column(String(50), nullable=True)
    location_sampler_type = Column(String(100), nullable=True)
    location_collection_window_hrs = Column(Float, nullable=True)
    location_collection_pull_ml = Column(Float, nullable=True)
    location_collection_step_min = Column(Float, nullable=True)
    location_collection_basis = Column(String(100), nullable=True)
    location_collection_type = Column(String(100), nullable=True)
    location_zipcode = Column(String(20), nullable=True)
    location_comment = Column(Text, nullable=True)