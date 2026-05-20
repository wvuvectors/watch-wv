# app/models/concentration.py

from sqlalchemy import Column, String, Text, ForeignKey
from app.database import Base

class Concentration(Base):
    __tablename__ = "concentration"
    
    concentration_id = Column(String(50), primary_key = True)
    sample_id = Column(String(50), ForeignKey("samples.sample_id"))
    concentration_batch_id = Column(String(50))
    concentration_location_in_batch = Column(String(50))
    concentration_comment = Column(Text)