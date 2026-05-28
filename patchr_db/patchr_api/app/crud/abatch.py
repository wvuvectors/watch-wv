# app/crud/abatch.py
# Defines helper functions to be used throughout app 

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_
from datetime import date, datetime

from app.models.abatch import aBatch

# Function to retrieve a single abatch by its assay_batch_id
def get_abatch_by_id(db: Session, assay_batch_id: str):
    return (
        db.query(aBatch)
        .filter(aBatch.assay_batch_id == assay_batch_id)
        .first()
    )

# List abatch entries
def list_abatch(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(aBatch)
        .offset(skip)
        .limit(limit)
        .all()
    )
    
# Dynamic Query
def query_abatch(
    db: Session,
    *,
    assay_batch_id: str | None = None,
    assay_date: date | None = None, 
    assay_reaction_ul: float | None = None, 
    assay_machine: str | None = None, 
    assay_amplification_method: str | None = None, 
    assay_amplification_method_lot_id: str | None = None, 
    assay_quantification_method: str | None = None, 
    assay_quantification_type: str | None = None, 
    assay_batch_record_version: str | None = None, 
    assay_method: str | None = None, 
    assay_method_lot_id: str | None = None, 
    assay_qx_manager_version: str | None = None, 
    assay_run_by: str | None = None, 
    skip: int = 0,
    limit: int = 10000,
):

    filter = []
    
    # ---- String / categorical filters ----
    if assay_batch_id is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_batch_id)) == assay_batch_id.strip().lower())
        
    if assay_machine is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_machine)) == assay_machine.strip().lower())
        
    if assay_amplification_method is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_amplification_method)) == assay_amplification_method.strip().lower())
        
    if assay_amplification_method_lot_id is not None:
        filters.append(func.lower(func.trim(aBatch.assay_amplification_method_lot_id)) == assay_amplification_method_lot_id.strip().lower())
        
    if assay_quantification_type is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_quantification_type)) == assay_quantification_type.strip().lower())
        
    if assay_batch_record_version is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_batch_record_version)) == assay_batch_record_version.strip().lower())
        
    if assay_method is not None:
        filters.append(func.lower(func.trim(aBatch.assay_method)) == assay_method.strip().lower())
        
    if assay_method_lot_id is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_method_lot_id)) == assay_method_lot_id.strip().lower())
        
    if assay_qx_manager_version is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_qx_manager_version)) == assay_qx_manager_version.strip().lower())
        
    if assay_run_by is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_run_by)) == assay_run_by.strip().lower())
        
    # ---- Numerical filters ----
    if assay_reaction_ul is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_reaction_ul)) == assay_reaction_ul.strip().lower())
        
    # ---- Date / datetime filters ---- 
    if assay_date is not None: 
        filters.append(func.lower(func.trim(aBatch.assay_date)) == assay_date.strip().lower())
        
    # Build statement
    stmt = select(aBatch).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt.scalars().all()
    return result 