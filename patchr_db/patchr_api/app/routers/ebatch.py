# app/routers/ebatch.py
# Handles all ebatch endpoints

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from datetime import date, datetime

from app.database import SessionLocal, get_db
from app.schemas.ebatch import eBatchSchema
from app.crud.ebatch import get_ebatch_by_id, list_ebatch, query_ebatch

# Create a router object
router = APIRouter(
    prefix="/ebatch",
    tags=["ebatch"]
)

# Get list of 100000 ebatch entries 
@router.get("/", response_model=list[eBatchSchema])
def read_ebatch(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    ebatch = list_ebatch(db=db, skip=skip, limit=limit)
    return ebatch
    
# Dynamic Querying
@router.get("/query", response_model=list[eBatchSchema])
def query_ebatch_endpoint(
    extraction_batch_id: str | None = Query(None), 
    extraction_date: date | None = Query(None), 
    extraction_input_ml: float | None = Query(None), 
    extraction_eluant: str | None = Query(None),
    extraction_machine: str | None = Query(None), 
    extraction_method: str | None = Query(None), 
    extraction_method_lot_id: str | None = Query(None), 
    extraction_output_ml: float | None = Query(None), 
    extraction_batch_record_version: str | None = Query(None), 
    extraction_run_by: str | None = Query(None), 
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query ebatch table with optional filters
    
    return query_ebatch(
        db=db,
        extraction_batch_id=extraction_batch_id,
        extraction_date=extraction_date,
        extraction_input_ml=extraction_input_ml,
        extraction_eluant=extraction_eluant,
        extraction_machine=extraction_machine,
        extraction_method=extraction_method,
        extraction_method_lot_id=extraction_method_lot_id,
        extraction_output_ml=extraction_output_ml,
        extraction_batch_record_version=extraction_batch_record_version,
        extraction_run_by=extraction_run_by,
        skip=skip,
        limit=limit
    )
    
# Get single ebatch id
@router.get("/{extraction_batch_id}", response_model=eBatchSchema)
def read_ebatch(extraction_batch_id: str, db: Session = Depends(get_db)):
    ebatch = get_ebatch_by_id(db=db, extraction_batch_id=extraction_batch_id)
    if not ebatch: 
        raise HTTPException(status_code=404, detail="eBatch record not found")
    return ebatch