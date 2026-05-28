# app/models/cbatch.py

from sqlalchemy import Column, String, Date, DateTime, Float, Text
from app.database import Base

class cBatch(Base):
    __tablename__ = "cbatch"
    
    concentration_batch_id = Column(String(50), primary_key = True)
    concentration_date = Column(Date, nullable=True)
    concentration_input_ml = Column(Float, nullable=True)
    concentration_machine = Column(String(100), nullable=True)
    concentration_method = Column(String(100), nullable=True)
    concentration_method_lot_id = Column(String(20), nullable=True)
    concentration_output_ml = Column(Float, nullable=True)
    concentration_run_by = Column(String(50), nullable=True)
    concentration_batch_record_version = Column(String(10), nullable=True)
    concentration_batch_comment = Column(Text, nullable=True)