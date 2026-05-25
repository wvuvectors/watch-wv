# app/routers/assay.py
# Handles all assay endpoints 

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.schemas.assay import AssaySchema
from app.crud.assay import get_assay_by_id, list_assays, query_assays

# Create router object
router = APIRouter(
    prefix="/assay",
    tags=["assay"]
)

# List assays
@router.get("/", response_model=list[AssaySchema])
def read_assays(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    assays = list_assays(db=db, skip=skip, limit=limit)
    return assays
    
# Dynamic querying
@router.get("/query", response_model=list[AssaySchema])
def query_assays_endpoint(
    assay_id: str | None = Query(None),
    extraction_id: str | None = Query(None),
    sample_id: str | None = Query(None),
    assay_batch_id: str | None = Query(None),
    assay_location_in_batch: str | None = Query(None), 
    assay_input_ul: float | None = Query(None), 
    assay_class: str | None = Query(None),
    assay_type: str | None = Query(None), 
    assay_target: str | None = Query(None), 
    assay_target_genetic_locus: str | None = Query(None), 
    assay_template: str | None = Query(None), 
    assay_target_macromolecule: str | None = Query(None), 
    assay_target_flourophore: str | None = Query(None), 
    assay_accepted_droplet: str | None = Query(None),
    assay_target_predicted_copies_per_ul_reaction: float | None = Query(None), 
    assay_target_copies_per_ul_reaction: float | None = Query(None), 
    assay_comment: str | None = Query(None), 
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query assays table with optional filters
    
    return query_assays(
        db=db,
        assay_id=assay_id,
        extraction_id=extraction_id,
        sample_id=sample_id,
        assay_batch_id=assay_batch_id,
        assay_location_in_batch=assay_location_in_batch,
        assay_input_ul=assay_input_ul,
        assay_class=assay_class,
        assay_type=assay_type,
        assay_target=assay_target,
        assay_target_genetic_locus=assay_target_genetic_locus,
        assay_template=assay_template,
        assay_target_macromolecule=assay_target_macromolecule,
        assay_target_flourophore=assay_target_flourophore,
        assay_accepted_droplet=assay_accepted_droplet,
        assay_target_predicted_copies_per_ul_reaction=assay_target_predicted_copies_per_ul_reaction,
        assay_target_copies_per_ul_reaction=assay_target_copies_per_ul_reaction,
        assay_comment=assay_comment,
        skip=skip,
        limit=limit
    )

# Get single assay_id
@router.get("/{assay_id}", response_model=AssaySchema)
def read_assays(assay_id: str, db: Session = Depends(get_db)):
    assay = get_assay_by_id(db=db, assay_id=assay_id)
    if not assay:
        raise HTTPException(status_code=404, detail="Assay record not found")
    return assay