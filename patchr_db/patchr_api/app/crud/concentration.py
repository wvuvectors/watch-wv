# app/crud/concentration.py
# Defines helper functions to be used throughout app

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_

from app.models.concentration import Concentration

# Function to retrieve a single concentration by its concentration_id
def get_concentration_by_id(db: Session, concentration_id: str):
    return (
        db.query(Concentration)
        .filter(Concentration.concentration_id == concentration_id)
        .first()
    )

# List concentrations
def list_concentrations(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(Concentration)
        .offset(skip)
        .limit(limit)
        .all()
    )

# Dynamic Query
def query_concentrations(
    db: Session, 
    *,
    concentration_id: str | None = None,
    sample_id: str | None = None,
    concentration_batch_id: str | None = None,
    concentration_location_in_batch: str | None = None,
    concentration_comment: str | None = None,
    skip: int = 0,
    limit: int = 10000,
):
    
    filters = []
    
    # ---- String / categorical filters ----
    if concentration_id is not None:
        filters.append(func.lower(func.trim(Concentration.concentration_id)) == concentration_id.strip().lower())
        
    if sample_id is not None:
        filters.append(func.lower(func.trim(Concentration.sample_id)) == sample_id.strip().lower())
        
    if concentration_batch_id is not None:
        filters.append(func.lower(func.trim(Concentration.concentration_batch_id)) == concentration_batch_id.strip().lower())
        
    if concentration_location_in_batch is not None:
        filters.append(func.lower(func.trim(Concentration.concentration_location_in_batch)) == concentration_location_in_batch.strip().lower())
        
    if concentration_comment is not None: 
        filters.append(func.lower(func.trim(Concentration.concentration_comment)) == concentration_comment.strip().lower())
        
    # Build statement
    stmt = select(Concentration).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result
        