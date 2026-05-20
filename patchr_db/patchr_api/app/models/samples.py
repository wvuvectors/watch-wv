# app/models/samples.py

from sqlalchemy import Column, String, Date, DateTime, Float, Text
from app.database import Base

class Samples(Base):
    __tablename__ = "samples"
    
    sample_id = Column(String(50), primary_key = True)
    sample_status = Column(String(50), nullable=True)
    location_id = Column(String(100))
    sample_event = Column(String(100), nullable=True)
    sample_qc = Column(String(100), nullable=True)
    sample_collection_start_datetime = Column(DateTime, nullable=True)
    sample_collection_end_datetime = Column(DateTime, nullable=True)
    sample_recovered_datetime = Column(DateTime, nullable=True)
    sample_flow = Column(Float, nullable=True)
    sample_received_by = Column(String(100), nullable=True)
    sample_received_date = Column(Date, nullable=True)
    sample_ph_lab = Column(Float, nullable=True)
    sample_comment = Column(Text, nullable=True)