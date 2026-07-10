# app/routers/county.py
# Handles all county endpoints 

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.schemas.county import CountySchema
from app.crud.county import get_county_by_id, list_counties, query_counties

# Create router object
router = APIRouter(
    prefix="/county",
    tags=["county"]
)

# List counties
@router.get("/", response_model=list[CountySchema])
def read_counties(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    counties = list_counties(db=db, skip=skip, limit=limit)
    return counties


# Dynamic querying
@router.get("/query", response_model=list[CountySchema])
def query_counties_endpoint(
    county_id: str | None = Query(None),
    county_labcode: str | None = Query(None),
    county_fips: str | None = Query(None),
    county_name: str | None = Query(None),
    county_population: float | None = Query(None),
    min_population: float | None = Query(None),
    max_population: float | None = Query(None),
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query county table with optional filters
    
    return query_counties(
        db=db,
        county_id=county_id,
        county_labcode=county_labcode,
        county_fips=county_fips,
        county_name=county_name,
        county_population=county_population,
        min_population=min_population,
        max_population=max_population,
        skip=skip,
        limit=limit
    )


# Get single county_id
@router.get("/{county_id}", response_model=CountySchema)
def read_county(county_id: str, db: Session = Depends(get_db)):
    county = get_county_by_id(db=db, county_id=county_id)
    if not county:
        raise HTTPException(status_code=404, detail="County record not found")
    return county