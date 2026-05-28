# app/routers/cbatch.py
# Handles all cbatch endpoints

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from datetime import datetime, datetime

from app.database import SessionLocal, get_db 
from app.schemas.cbatch import cBatchSchema
from app.crud.cbatch import get_cbatch_by_id, list_cbatch, query_cbatch

# Create a router object
router = APIRouter(
    prefix="/cbatch",
    tags=["cbatch"]
)

# Get list of 100000 cbatch entries
@router.get("/", response_model=list[cBatchSchema])
def read_cbatch(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    cbatch = list_cbatch(db=db, skip=skip, limit=limit)
    return cbatch
    
# Dynamic querying
@router.get("/query", response_model=list[cBatchSchema])
def query_cbatch_endpoint(
    concentration_batch_id: str | None = Query(None),
    concentration_date: date | None = Query(None), 
    concentration_input_ml: float | None = Query(None), 
    concentration_machine: str | None = Query(None), 
    concentration_method: str | None = Query(None), 
    concentration_method_lot_id: str | None = Query(None), 
    concentration_output_ml: float | None = Query(None), 
    concentration_run_by: str | None = Query(None), 
    concentration_batch_record_version: str | None = Query(None), 
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query cbatch table with optional filters
    
    return query_cbatch(
        db=db,
        concentration_batch_id=concentration_batch_id,
        concentration_date=concentration_date,
        concentration_input_ml=concentration_input_ml,
        concentration_machine=concentration_machine,
        concentration_method=concentration_method,
        concentration_method_lot_id=concentration_method_lot_id,
        concentration_output_ml=concentration_output_ml,
        concentration_run_by=concentration_run_by,
        concentration_batch_record_version=concentration_batch_record_version,
        skip=skip,
        limit=limit
    )
    
# Get single cbatch_id
@router.get("/{concentration_batch_id}", response_model=cBatchSchema)
def read_cbatch(concentration_batch_id: str, db: Session = Depends(get_db)):
    cbatch = get_cbatch_by_id(db=db, concentration_batch_id=concentration_batch_id)
    if not cbatch: 
        raise HTTPException(status_code=404, detail="cBatch record not found")
    retrn cbatch