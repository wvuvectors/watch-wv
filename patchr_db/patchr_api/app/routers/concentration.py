# app/routers/concentration.py
# Handles all concentration endpoints

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.schemas.concentration import ConcentrationSchema
from app.crud.concentration import get_concentration_by_id, list_concentrations, query_concentrations

# Create router object
router = APIRouter(
    prefix="/concentration",
    tags=["concentration"]
)

# List concentrations
@router.get("/", response_model=list[ConcentrationSchema])
def read_concentrations(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    concentrations = list_concentrations(db=db, skip=skip, limit=limit)
    return concentrations

# Dynamic querying
@router.get("/query", response_model=list[ConcentrationSchema])
def query_concentrations_endpoint(
    concentration_id: str | None = Query(None),
    sample_id: str | None = Query(None),
    concentration_batch_id: str | None = Query(None),
    concentration_location_in_batch: str | None = Query(None),
    concentration_comment: str | None = Query(None),
    skip: int = Query(0),
    limit: int = Query(100), 
    db: Session = Depends(get_db)
):
    # Query concentration table with optional filters
    
    return query_concentrations(
        db=db,
        concentration_id=concentration_id,
        sample_id=sample_id,
        concentration_batch_id=concentration_batch_id,
        concentration_location_in_batch=concentration_location_in_batch,
        concentration_comment=concentration_comment,
        skip=skip,
        limit=limit
    )
    
# Get single concentration
@router.get("/{concentration_id}", response_model=ConcentrationSchema)
def read_concentration(concentration_id: str, db: Session = Depends(get_db)):
    concentration = get_concentration_by_id(db=db, concentration_id=concentration_id)
    if not concentration:
        raise HTTPException(status_code=404, detail="Concentration record not found")
    return concentration