# app/routers/location.py
# Handles all location endpoints 

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.database import SessionLocal, get_db
from app.schemas.location import LocationSchema
from app.crud.location import get_location_by_id, list_locations, query_locations

# Create router object
router = APIRouter(
    prefix="/location",
    tags=["location"]
)

# List locations
@router.get("/", response_model=list[LocationSchema])
def read_locations(skip: int = 0, limit: int = Query(1000, le=100000), db: Session = Depends(get_db)):
    locations = list_locations(db=db, skip=skip, limit=limit)
    return locations


# Dynamic querying
@router.get("/query", response_model=list[LocationSchema])
def query_locations_endpoint(
    location_id: str | None = Query(None),
    sample_code_prefix: str | None = Query(None),
    location_primary_lab: str | None = Query(None),
    location_status: str | None = Query(None),
    location_common_name: str | None = Query(None),
    location_category: str | None = Query(None),
    location_group: str | None = Query(None),
    location_primary_wwtp_id: str | None = Query(None),
    location_counties_served: str | None = Query(None),
    location_population_served: str | None = Query(None),
    location_sampler_type: str | None = Query(None),
    location_collection_basis: str | None = Query(None),
    location_collection_type: str | None = Query(None),
    location_zipcode: str | None = Query(None),
    location_comment: str | None = Query(None),
    location_lng: float | None = Query(None),
    location_lat: float | None = Query(None),
    min_collection_window_hrs: float | None = Query(None),
    max_collection_window_hrs: float | None = Query(None),
    min_collection_pull_ml: float | None = Query(None),
    max_collection_pull_ml: float | None = Query(None),
    min_collection_step_min: float | None = Query(None),
    max_collection_step_min: float | None = Query(None),
    skip: int = Query(0),
    limit: int = Query(1000),
    db: Session = Depends(get_db)
):

    # Query locations table with optional filters
    
    return query_locations(
        db=db,
        location_id=location_id,
        sample_code_prefix=sample_code_prefix,
        location_primary_lab=location_primary_lab,
        location_status=location_status,
        location_common_name=location_common_name,
        location_category=location_category,
        location_group=location_group,
        location_primary_wwtp_id=location_primary_wwtp_id,
        location_counties_served=location_counties_served,
        location_population_served=location_population_served,
        location_sampler_type=location_sampler_type,
        location_collection_basis=location_collection_basis,
        location_collection_type=location_collection_type,
        location_zipcode=location_zipcode,
        location_comment=location_comment,
        location_lng=location_lng,
        location_lat=location_lat,
        min_collection_window_hrs=min_collection_window_hrs,
        max_collection_window_hrs=max_collection_window_hrs,
        min_collection_pull_ml=min_collection_pull_ml,
        max_collection_pull_ml=max_collection_pull_ml,
        min_collection_step_min=min_collection_step_min,
        max_collection_step_min=max_collection_step_min,
        skip=skip,
        limit=limit
    )


# Get single location_id
@router.get("/{location_id}", response_model=LocationSchema)
def read_location(location_id: str, db: Session = Depends(get_db)):
    location = get_location_by_id(db=db, location_id=location_id)
    if not location:
        raise HTTPException(status_code=404, detail="Location record not found")
    return location