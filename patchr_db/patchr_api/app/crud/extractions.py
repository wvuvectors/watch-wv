# app/crud/extractions.py
# Defines helper functions to be used throughout app

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_

from app.models.extractions import Extractions

# Function to retrieve a single extraction by its extraction_id
def get_extraction_by_id(db: Session, extractions_id: str):
    return (
        db.query(Extractions)
        .filter(Extractions.extraction_id == extraction_id)
        .first()
    )
    
# List extractions
def list_extractions(db: Session, skip: int = 0, limit: int = 10000):
    return (
        db.query(Extractions)
        .offset(skip)
        .limit(limit)
        .all()
    )

# Dynamic Query
def query_extractions(
    db: Session,
    *,
    extraction_id: str | None = None,
    concentration_id: str | None = None,
    extraction_batch_id: str | None = None,
    extraction_location_in_batch: str | None = None,
    extraction_location_in_storage: str | None = None,
    extraction_comment: str | None = None,
    skip: int = 0,
    limit: int = 100,
):

    filters = []
    
    # ---- String / categorical filters ----
    if extraction_id is not None:
        filters.append(func.lower(func.trim(Extractions.extraction_id)) == extraction_id.strip().lower())
        
    if concentration_id is not None:
        filters.append(func.lower(func.trim(Extractions.concentration_id)) == concentration_id.strip().lower())
        
    if extraction_batch_id is not None:
        filters.append(func.lower(func.trim(Extractions.extraction_batch_id)) == extraction_batch_id.strip().lower())
        
    if extraction_location_in_batch is not None:
        filters.append(func.lower(func.trim(Extractions.extraction_location_in_batch)) == extraction_location_in_batch.strip().lower())
        
    if extraction_location_in_storage is not None:
        filters.append(func.lower(func.trim(Extractions.extraction_location_in_storage)) == extraction_location_in_storage.strip().lower())
        
    if extraction_comment is not None:
        filters.append(func.lower(func.trim(Extractions.extraction_comment)) == extraction_comment.strip().lower())
    
    # Build statement
    stmt = select(Extractions).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result