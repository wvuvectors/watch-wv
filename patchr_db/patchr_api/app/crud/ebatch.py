# app/crud/ebatch.py
# Defines helper functions to be used throughout app 

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_
from datetime import date, datetime

from app.models.ebatch import eBatch

# Function to retrieve a single ebatch by its extraction_batch_id
def get_ebatch_by_id(db: Session, extraction_batch_id: str):
    return (
        db.query(eBatch)
        .filter(eBatch.extraction_batch_id == extraction_batch_id)
        .first()
    )

# List ebatch entries
def list_ebatch(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(eBatch)
        .offset(skip)
        .limit(limit)
        .all()
    )
    
# Dynamic Query
def query_ebatch(
    db: Session,
    *,
    extraction_batch_id: str | None = None, 
    start_date: date | None = None, 
    end_date: date | None = None,
    extraction_input_ul: float | None = None, 
    extraction_eluant: str | None = None, 
    extraction_machine: str | None = None,
    extraction_method: str | None = None, 
    extraction_method_lot_id: str | None = None, 
    min_output_ul: float | None = None, 
    max_output_ul: float | None = None,
    extraction_batch_record_version: str | None = None, 
    extraction_run_by: str | None = None, 
    skip: int = 0,
    limit: int = 10000,
):

    filters = []
    
    # ----- String / categorical filters ----
    if extraction_batch_id is not None:
        filters.append(func.lower(func.trim(eBatch.extraction_batch_id)) == extraction_batch_id.strip().lower())
        
    if extraction_eluant is not None: 
        filters.append(func.lower(func.trim(eBatch.extraction_eluant)) == extraction_eluant.strip().lower())
        
    if extraction_machine is not None: 
        filters.append(func.lower(func.trim(eBatch.extraction_machine)) == extraction_machine.strip().lower())
        
    if extraction_method is not None: 
        filters.append(func.lower(func.trim(eBatch.extraction_method)) == extraction_method.strip().lower())
       
    if extraction_method_lot_id is not None: 
        filters.append(func.lower(func.trim(eBatch.extraction_method_lot_id)) == extraction_method_lot_id.strip().lower())
        
    if extraction_batch_record_version is not None:
        filters.append(func.lower(func.trim(eBatch.extraction_batch_record_version)) == extraction_batch_record_version.strip().lower())
        
    if extraction_run_by is not None:
        filterse.append(func.lower(func.trim(eBatch.extraction_run_by)) == extraction_run_by.strip().lower())
        
    # ---- Numerical filters ----
    if extraction_input_ul is not None:
        filters.append(eBatch.extraction_input_ul == extraction_input_ul)
        
    if min_output_ul is not None:
        filters.append(eBatch.extraction_output_ul >= min_output_ul)
        
    if max_output_ul is not None: 
        filters.append(eBatch.extraction_output_ul <= max_output_ul)
        
    # ---- Date / datetime filters ----
    if start_date is not None:
        filters.append(eBatch.extraction_date >= start_date)
        
    if end_date is not None:
        filters.append(eBatch.extraction_date <= end_date)
        
    # Build statement
    stmt = select(eBatch).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result