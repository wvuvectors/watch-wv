# app/models/extractions.py

from sqlalchemy import Column, String, Text, ForeignKey
from app.database import Base

class Extractions(Base):
    __tablename__ = "extractions"
    
    extraction_id = Column(String(50), primary_key = True)
    concentration_id = Column(String(50), ForeignKey("concentration.concentration_id"))
    extraction_batch_id = Column(String(50))
    extraction_location_in_batch = Column(String(50))
    extraction_location_in_storage = Column(String(50))
    extraction_comment = Column(Text)