# app/models/assays.py

from sqlalchemy import Column, String, Text, Integer, ForeignKey
from app.database import Base

class Assays(Base):
    __tablename__ = "assays"
    
    assay_id = Column(String(50), primary_key = True)
    extraction_id = Column(String(50), nullable=True)
    sample_id = Column(String(50), nullable=True)
    assay_batch_id = Column(String(50), nullable=True)
    assay_location_in_batch = Column(String(50), nullable=True)
    assay_input_ul = Column(Float, nullable=True)
    assay_class = Column(String(50), nullable=True)
    assay_type = Column(String(50), nullable=True)
    assay_target = Column(String(50), nullable=True)
    assay_target_genetic_locus = Column(String(50), nullable=True)
    assay_template = Column(String(50), nullable=True)
    assay_target_macromolecule = Column(String(50), nullable=True)
    assay_target_fluorophore = Column(String(50), nullable=True)
    assay_accepted_droplets = Column(Float, nullable=True)
    assay_target_predicted_copies_per_ul_reaction = Column(Float, nullable=True)
    assay_target_copies_per_ul_reaction = Column(Float, nullable=True)
    assay_comment = Column(Text, nullable=True)