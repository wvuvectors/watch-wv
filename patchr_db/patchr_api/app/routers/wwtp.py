# app/routers/wwtp.py
# Handles all WWTP endpoints 

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.schemas.wwtp import WWTPSchema
from app.crud.wwtp import get_wwtp_by_id, list_wwtps, query_wwtps

# Create router object
router = APIRouter(
    prefix="/wwtp",
    tags=["wwtp"]
)

# List WWTPs
@router.get("/", response_model=list[WWTPSchema])
def read_wwtps(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    wwtps = list_wwtps(db=db, skip=skip, limit=limit)
    return wwtps


# Dynamic querying
@router.get("/query", response_model=list[WWTPSchema])
def query_wwtps_endpoint(
    wwtp_id: str | None = Query(None),
    wwtp_site_id: str | None = Query(None),
    wwtp_common_name: str | None = Query(None),
    wwtp_authority_name: str | None = Query(None),
    wwtp_counties_served: str | None = Query(None),
    wwtp_epaid_id: str | None = Query(None),
    wwtp_cwns_id: str | None = Query(None),
    wwtp_capacity_mgd: float | None = Query(None),
    wwtp_population_served: float | None = Query(None),
    min_capacity_mgd: float | None = Query(None),
    max_capacity_mgd: float | None = Query(None),
    min_population_served: float | None = Query(None),
    max_population_served: float | None = Query(None),
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query WWTP table with optional filters
    
    return query_wwtps(
        db=db,
        wwtp_id=wwtp_id,
        wwtp_site_id=wwtp_site_id,
        wwtp_common_name=wwtp_common_name,
        wwtp_authority_name=wwtp_authority_name,
        wwtp_counties_served=wwtp_counties_served,
        wwtp_epaid_id=wwtp_epaid_id,
        wwtp_cwns_id=wwtp_cwns_id,
        wwtp_capacity_mgd=wwtp_capacity_mgd,
        wwtp_population_served=wwtp_population_served,
        min_capacity_mgd=min_capacity_mgd,
        max_capacity_mgd=max_capacity_mgd,
        min_population_served=min_population_served,
        max_population_served=max_population_served,
        skip=skip,
        limit=limit
    )


# Get single wwtp_id
@router.get("/{wwtp_id}", response_model=WWTPSchema)
def read_wwtp(wwtp_id: str, db: Session = Depends(get_db)):
    wwtp = get_wwtp_by_id(db=db, wwtp_id=wwtp_id)
    if not wwtp:
        raise HTTPException(status_code=404, detail="WWTP record not found")
    return wwtp