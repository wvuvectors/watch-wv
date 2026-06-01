# app/crud/cbatch.py
# Defines helper functions to be used throughout app

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_
from datetime import date, datetime

from app.models.cbatch import cBatch

# Function to retrieve a single cbatch by its concentration_batch_id
def get_cbatch_by_id(db: Session, concentration_batch_id: str):
    return (
        db.query(cBatch)
        .filter(cBatch.concentration_batch_id == concentration_batch_id)
        .first()
    )
    
# List cbatch entries 
def list_cbatch(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(cBatch)
        .offset(skip)
        .limit(limit)
        .all()
    )
    
# Dynamic Query
def query_cbatch(
    db: Session, 
    *, 
    concentration_batch_id: str | None = None,
    start_date: date | None = None,
    end_date: date | None = None,
    concentration_input_ml: float | None = None, 
    concentration_machine: str | None = None, 
    concentration_method: str | None = None, 
    concentration_method_lot_id: str | None = None, 
    min_output_ml: float | None = None, 
    max_output_ml: float | None = None,
    concentration_run_by: str | None = None, 
    concentration_batch_record_version: str | None = None, 
    skip: int = 0,
    limit: int = 10000,
):

    filters = []
    
    # ---- String / categorical filters ----
    if concentration_batch_id is not None: 
        filters.append(func.lower(func.trim(cBatch.concentration_batch_id)) == concentration_batch_id.strip().lower())
        
    if concentration_machine is not None:
        filters.append(func.lower(func.trim(cBatch.concentration_machine)) == concentration_machine.strip().lower())
        
    if concentration_method is not None:
        filters.append(func.lower(func.trim(cBatch.concentration_method)) == concentration_method.strip().lower())
        
    if concentration_method_lot_id is not None:
        filters.append(func.lower(func.trim(cBatch.concentration_method_lot_id)) == concentration_method_lot_id.strip().lower())
        
    if concentration_batch_record_version is not None:
        filters.append(func.lower(func.trim(cBatch.concentration_batch_record_version)) == concentration_batch_record_version.strip().lower())
        
    # ---- Numeric filters ----
    if concentration_input_ml is not None:
        filters.append(cBatch.concentration_input_ml == concentration_input_ml)
        
    if min_output_ml is not None:
        filters.append(cBatch.concentration_output_ml >= min_output_ml)
        
    if max_output_ml is not None:
        filters.append(cBatch.concentration_output_ml <= max_output_ml)

    # ---- Date filters ----
    if start_date is not None:
        filters.append(cBatch.concentration_date >= start_date)
        
    if end_date is not None:
        filters.append(cBatch.concentration_date <= end_date)
     
    # Build statement
    stmt = select(cBatch).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result