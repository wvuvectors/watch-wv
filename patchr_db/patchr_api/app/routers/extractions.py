# app/routers/extractions.py
# Handles all extractions endpoints

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.schemas.extractions import ExtractionsSchema
from app.crud.extractions import get_extraction_by_id, list_extractions, query_extractions

# Create router object
router = APIRouter(
    prefix="/extractions",
    tags=["extractions"]
)

# List extractions
@router.get("/", response_model=list[ExtractionsSchema])
def read_extractions(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    extractions = list_extractions(db=db, skip=skip, limit=limit)
    return extractions

# Dynamic querying
@router.get("/query", response_model=list[ExtractionsSchema])
def query_extractions_endpoint(
    extraction_id: str | None = Query(None),
    concentration_id: str | None = Query(None),
    extraction_batch_id: str | None = Query(None),
    extraction_location_in_batch: str | None = Query(None),
    extraction_location_in_storage: str | None = Query(None),
    extraction_comment: str | None = Query(None),
    skip: int = Query(0),
    limit: int = Query(100),
    db: Session = Depends(get_db)
):
    
    # Query extractions table with optional filters
    
    return query_extractions(
        db=db,
        extraction_id=extraction_id,
        concentration_id=concentration_id,
        extraction_batch_id=extraction_batch_id,
        extraction_location_in_batch=extraction_location_in_batch,
        extraction_location_in_storage=extraction_location_in_storage,
        extraction_comment=extraction_comment,
        skip=skip,
        limit=limit
    )
    
# Get single extraction_batch_id
@router.get("/{extraction_id}", response_model=ExtractionsSchema)
def read_extractions(extraction_id: str, db: Session = Depends(get_db)):
    extraction = get_extraction_by_id(db=db, extraction_id=extraction_id)
    if not extraction:
        raise HTTPException(status_code=404, detail="Extraction record not found")
    return extraction