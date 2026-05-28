# app/models/abatch.py

from sqlalchemy import Column, String, Date, DateTime, Float, Text
from app.database import Base 

class aBatch(Base):
    __tablename__ = "abatch"
    
    assay_batch_id = Column(String(50), primary_key = True)
    assay_date = Column(Date, nullable=True)
    assay_reaction_ul = Column(Float, nullable=True)
    assay_machine = Column(String(50), nullable=True)
    assay_amplification_method = Column(String(100), nullable=True)
    assay_amplification_method_lot_id = Column(String(50), nullable=True)
    assay_quantification_method = Column(String(50), nullable=True)
    assay_quantification_type = Column(String(50), nullable=True)
    assay_batch_record_version = Column(String(50), nullable=True) 
    assay_method = Column(String(50), nullable=True)
    assay_method_lot_id = Column(String(50), nullable=True)
    assay_qx_manager_version = Column(String(100), nullable=True)
    assay_run_by = Column(String(50), nullable=True)
    assay_batch_comment = Column(Text, nullable=True)