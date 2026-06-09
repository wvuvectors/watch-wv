# app/crud/samples.py
# Defines helper functions to be used throughout app

from sqlalchemy.orm import Session  # SQLAlchemy ORM Session 
from sqlalchemy import select, func, and_
from datetime import date, datetime

from app.models.samples import Samples  # SQLAlchemy Samples table model

# Function to retrieve a single sample by its sample_id
def get_sample_by_id(db: Session, sample_id: str):
    """
    Query the database to return the first Sample record matching the given sample_id.
    Returns None if no matching sample is found.
    """
    return (
        db.query(Samples)                 # Start a query on the Samples table
        .filter(Samples.sample_id == sample_id)  # Filter by sample_id
        .first()                         # Return only the first result
    )

# Function to list multiple samples, with an optional limit (default 100)
def list_samples(db: Session, skip: int = 0, limit: int = 10000):
    """
    Query the database to return a list of Sample records.
    Supports optional skip and limit for pagination.
    """
    return (
        db.query(Samples)        # Query the Samples table
        .offset(skip)            # Skip a number of rows (useful for pagination)
        .limit(limit)            # Limit the number of rows returned
        .all()                   # Return as a list of Sample objects
    )
# Flexible conditional query across multiple columns
def query_samples(
    db: Session,
    *,
    sample_id: str | None = None,
    status: str | None = None,
    location_id: str | None = None,
    sample_event: str | None = None,
    qc: str | None = None,
    collection_start: datetime | None = None,
    collection_end: datetime | None = None,
    recovered_start: datetime | None = None,
    recovered_end: datetime | None = None,
    min_flow: float | None = None,
    max_flow: float | None = None,
    received_by: str | None = None,
    received_start_date: date | None = None,
    received_end_date: date | None = None,
    min_ph: float | None = None,
    max_ph: float | None = None,
    skip: int = 0,
    limit: int = 10000,
):
    filters = []
    
    """
    Dynamically query the Samples table using optional filters.

    Only filters that are not None are applied.
    Designed to scale as the frontend adds more query controls.
    """

    # ---- String / categorical filters ----
    if sample_id is not None:
        filters.append(Samples.sample_id.ilike(f"%{sample_id}%"))
    
    if status is not None:
        filters.append(Samples.sample_status == status)
    
    if location_id is not None:
        filters.append(Samples.location_id == location_id)
    
    if sample_event is not None:
        filters.append(Samples.sample_event == sample_event)
    
    if qc is not None:
        filters.append(Samples.sample_qc == qc)
    
    if received_by is not None:
        filters.append(Samples.sample_received_by == received_by)


    # ---- Numeric filters ----
    if min_flow is not None:
        filters.append(Samples.sample_flow >= min_flow)
    
    if max_flow is not None:
        filters.append(Samples.sample_flow <= max_flow)
    
    if min_ph is not None:
        filters.append(Samples.sample_ph_lab >= min_ph)
    
    if max_ph is not None:
        filters.append(Samples.sample_ph_lab <= max_ph)

        
    # ---- Date filters ----
    if received_start_date is not None:
        filters.append(Samples.sample_received_date >= received_start_date)
    
    if received_end_date is not None:
        filters.append(Samples.sample_received_date <= received_end_date)

    # ---- Datetime filters ----
    if collection_start is not None:
        filters.append(Samples.sample_collection_start_datetime >= collection_start)
    
    if collection_end is not None:
        filters.append(Samples.sample_collection_end_datetime <= collection_end)
    
    if recovered_start is not None:
        filters.append(Samples.sample_recovered_datetime >= recovered_start)
    
    if recovered_end is not None:
        filters.append(Samples.sample_recovered_datetime <= recovered_end)

    stmt = select(Samples).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result