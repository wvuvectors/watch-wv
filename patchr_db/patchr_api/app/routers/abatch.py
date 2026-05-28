# app/routers/abatch.py
# Handles all abatch endpoints

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from datetime import date, datetime

from app.database import SessionLocal, get_db
from app.schemas.abatch import aBatchSchema
from app.crud.abatch import get_abatch_by_id, list_abatch, query_abatch

# Create a router object
router = APIRouter(
    prefix="/abatch",
    tags=["abatch"]
)

# Get list of 100000 abatch entries
@router.get("/", response_model=list[aBatchSchema])
def read_abatch(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    abatch = list_abatch(db=db, skip=skip, limit=limit)
    return abatch
    
# Dynamic Querying
@router.get("/query", response_model=list[aBatchSchema])
def query_abatch_endpoint(
    assay_batch_id: str | None = Query(None),
    assay_date: date | None = Query(None), 
    assay_reaction_ul: float | None = Query(None), 
    assay_machine: str | None = Query(None), 
    assay_amplification_method: str | None = Query(None), 
    assay_amplification_method_lot_id: str | None = Query(None), 
    assay_quantification_method: str | None = Query(None), 
    assay_quantification_type: str | None = Query(None), 
    assay_batch_record_version: str | None = Query(None), 
    assay_method: str | None = Query(None), 
    assay_method_lot_id: str | None = Query(None), 
    assay_qx_manager_version: str | None = Query(None), 
    assay_run_by: str | None = Query(None), 
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query abatch table with optional filters 
    
    return query_abatch(
        db=db,
        assay_batch_id=assay_batch_id,
        assay_date=assay_date,
        assay_reaction_ul=assay_reaction_ul,
        assay_machine=assay_machine,
        assay_amplification_method=assay_amplification_method,
        assay_amplification_method_lot_id=assay_amplification_method_lot_id,
        assay_quantification_method=assay_quantification_method,
        assay_quantification_type=assay_quantification_type,
        assay_batch_record_version=assay_batch_record_version,
        assay_method=assay_method,
        assay_method_lot_id=assay_method_lot_id,
        assay_qx_manager_version=assay_qx_manager_version,
        assay_run_by=assay_run_by,
        skip=skip,
        limit=limit
    )

# Get single abatch id
@router.get("/{assay_batch_id}", response_mode=aBatchSchema)
def read_abatch(assay_batch_id: str, db: Session = Depends(get_db)):
    abatch = get_abatch_by_id(db=db, assay_batch_id=assay_batch_id)
    if not abatch: 
        raise HTTPException(status_code=404, detail="aBatch record not found")
    return abatch