# app/routers/samples.py
# Handles all endpoints
# Queries tables from MySQL DB, serializing w Pydantic schemas, and exposes API routes

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from datetime import datetime, date

from app.database import SessionLocal, get_db # FastAPI dependency for db sessions
from app.schemas.samples import SamplesSchema  # Pydantic schema for serialization
from app.crud.samples import get_sample_by_id, list_samples, query_samples # CRUD helper functions

# Create a router object
router = APIRouter(
    prefix="/samples",
    tags=["samples"]
)

# Get list of 10000 samples 
@router.get("/", response_model=list[SamplesSchema])
def read_samples(skip: int = 0, limit: int = Query(1000, le=10000), db: Session = Depends(get_db)):
    """
    Returns a list of samples from the database.
    
    - `skip`: Number of records to skip (for pagination)
    - `limit`: Maximum number of records to return
    """
    samples = list_samples(db=db, skip=skip, limit=limit)
    return samples
    
# Dynamic querying
@router.get("/query", response_model=list[SamplesSchema])
def query_samples_endpoint(
    sample_id: str | None = Query(None),
    status: str | None = Query(None, alias="sample_status"),
    location_id: str | None = Query(None),
    sample_event: str | None = Query(None),
    qc: str | None = Query(None, alias="sample_qc"),
    collection_start: date | None = Query(None),
    collection_end: date | None = Query(None),
    recovered_start: date | None = Query(None),
    recovered_end: date | None = Query(None),
    min_flow: float | None = Query(None),
    max_flow: float | None = Query(None),
    received_by: str | None = Query(None),
    received_start_date: date | None = Query(None),
    received_end_date: date | None = Query(None),
    min_ph: float | None = Query(None),
    max_ph: float | None = Query(None),
    skip: int = Query(0),
    limit: int = Query(100),
    db: Session = Depends(get_db)
):
    """
    Query Samples table with optional filters.
    All parameters are optional. Date inputs should be in YYYY-MM-DD format.
    """

    # Convert date inputs to datetime objects for CRUD function
    collection_start_dt = datetime.combine(collection_start, datetime.min.time()) if collection_start else None
    collection_end_dt = datetime.combine(collection_end, datetime.max.time()) if collection_end else None
    recovered_start_dt = datetime.combine(recovered_start, datetime.min.time()) if recovered_start else None
    recovered_end_dt = datetime.combine(recovered_end, datetime.max.time()) if recovered_end else None

    return query_samples(
        db=db,
        sample_id=sample_id,
        status=status,
        location_id=location_id,
        sample_event=sample_event,
        qc=qc,
        collection_start=collection_start_dt,
        collection_end=collection_end_dt,
        recovered_start=recovered_start_dt,
        recovered_end=recovered_end_dt,
        min_flow=min_flow,
        max_flow=max_flow,
        received_by=received_by,
        received_start_date=received_start_date,
        received_end_date=received_end_date,
        min_ph=min_ph,
        max_ph=max_ph,
        skip=skip,
        limit=limit
    )
    
# Get a single sample by sample_id
@router.get("/{sample_id}", response_model=SamplesSchema)
def read_sample(sample_id: str, db: Session = Depends(get_db)):
    """
    Returns a single sample by `sample_id`.
    
    Raises a 404 error if the sample does not exist.
    """
    sample = get_sample_by_id(db=db, sample_id=sample_id)
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    return sample
    