# app/models/ebatch.py

from sqlalchemy import Column, String, Date, Datetime, Float, Text
from app.database import Base 

class eBatch(Base):
    __tablename__ = "ebatch"
    
    extraction_batch_id = Column(String(50), primary_key = True)
    extraction_date = Column(Date, nullable=True)
    extraction_input_ul = Column(Float, nullable=True)
    extraction_eluant = Column(String(100), nullable=True)
    extraction_machine = Column(String(100), nullable=True)
    extraction_method = Column(String(100), nullable=True)
    extraction_method_lot_id = Column(String(20), nullable=True)
    extraction_output_ul = Column(Float, nullable=True)
    extraction_batch_record_version = Column(String(10), nullable=True)
    extraction_run_by = Column(String(50, nullable=True)
    extraction_batch_comment = Column(String(100), nullable=True)