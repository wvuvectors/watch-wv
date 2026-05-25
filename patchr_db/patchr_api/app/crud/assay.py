# app/crud/assay.py
# Defines helper functions to be used throughout app 

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_

from app.models.assay import Assay

# Function to retrieve a single assay by its assay_id
def get_assay_by_id(db: Session, assay_id: str):
    return (
        db.query(Assay)
        .filter(Assay.assay_id == assay_id)
        .first()
    )

# List assays
def list_assays(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(Assay)
        .offset(skip)
        .limit(limit)
        .all()
    )

# Dynamic Query
def query_assays(
    db: Session,
    *,
    assay_id: str | None = None,
    extraction_id: str | None = None,
    sample_id: str | None = None,
    assay_batch_id: str | None = None,
    assay_location_in_batch: str | None = None,
    assay_input_ul: float | None = None,
    assay_class: str | None = None,
    assay_type: str | None = None,
    assay_target: str | None = None,
    assay_target_genetic_locus: str | None = None,
    assay_template: str | None = None,
    assay_target_macromolecule: str | None = None,
    assay_target_flourophore: str | None = None,
    assay_accepted_droplets: float | None = None,
    assay_target_predicted_copies_per_ul_reaction: float | None = None,
    assay_target_copies_per_ul_reaction: float | None = None,
    assay_comment: str | None = None,
    skip: int = 0,
    limit: int = 10000,
):
    filters = []
    
    # ---- String / categorical filters ----
    if assay_id is not None:
        filters.append(func.lower(func.trim(Assay.assay_id)) == assay_id.strip().lower())
    
    if extraction_id is not None: 
        filters.append(func.lower(func.trim(Assay.extraction_id)) == extraction_id.strip().lower())
        
    if sample_id is not None: 
        filters.append(func.lower(func.trim(Assay.samples_id)) == sample_id.strip().lower())
        
    if assay_batch_id is not None: 
        filters.append(func.lower(func.trim(Assay.assay_batch_id)) == assay_batch_id.strip().lower())
    
    if assay_location_in_batch is not None: 
        filters.append(func.lower(func.trim(Assay.assay_location_in_batch)) == assay_location_in_batch.strip().lower())
    
    if assay_input_ul is not None:
        filters.append(func.lower(func.trim(Assay.assay_input_ul)) == assay_input_ul.strip().lower())
        
    if assay_class is not None:
        filters.append(func.lower(func.trim(Assay.assay_class)) == assay_class.strip().lower())
        
    if assay_type is not None: 
        filters.append(func.lower(func.trim(Assay.assay_type)) == assay_type.strip().lower())
        
    if assay_target is not None:
        filters.append(func.lower(func.trim(Assay.assay_target)) == assay_target.strip().lower())
        
    if assay_target_genetic_locus is not None:
        filters.append(func.lower(func.trim(Assay.assay_target_genetic_locus)) == assay_target_genetic_locus.strip().lower())
        
    if assay_template is not None: 
        filters.append(func.lower(func.trim(Assay.assay_template)) == assay_template.strip().lower())
        
    if assay_target_macromolecule is not None:
        filters.append(func.lower(func.trim(Assay.assay_target_macromolecule)) == assay_target_macromolecule.strip().lower())
        
    if assay_target_flourophore is not None:
        filters.append(func.lower(func.trim(Assay.assay_target_flourophore)) == assay_target_flourophore.strip().lower())
        
    if assay_accepted_droplets is not None:
        filters.append(func.lower(func.trim(Assay.assay_accepted_droplets)) == assay_accepted_droplets.strip().lower())
        
    if assay_target_predicted_copies_per_ul_reaction is not None: 
        filters.append(func.lower(func.trim(Assay.assay_target_predicted_copies_per_ul_reaction)) == assay_target_predicted_copies_per_ul_reaction.strip().lower())
        
    if assay_target_copies_per_ul_reaction is not None:
        filters.append(func.lower(func.trim(Assay.assay_target_copies_per_ul_reaction)) == assay_target_copies_per_ul_reaction.strip().lower())
        
    if assay_comment is not None:
        filters.append(func.lower(func.trim(Assay.assay_comment)) == assay_comment.strip().lower())
    
    # Build statement
    stmt = select(Assays).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result