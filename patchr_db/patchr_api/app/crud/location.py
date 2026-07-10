# app/crud/location.py
# Defines helper functions to be used throughout app 

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_

from app.models.location import Location

# Function to retrieve a single location by its location_id
def get_location_by_id(db: Session, location_id: str):
    return (
        db.query(Location)
        .filter(Location.location_id == location_id)
        .first()
    )

# List locations
def list_locations(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(Location)
        .offset(skip)
        .limit(limit)
        .all()
    )

# Dynamic Query
def query_locations(
    db: Session,
    *,
    location_id: str | None = None,
    sample_code_prefix: str | None = None,
    location_primary_lab: str | None = None,
    location_status: str | None = None,
    location_common_name: str | None = None,
    location_category: str | None = None,
    location_group: str | None = None,
    location_primary_wwtp_id: str | None = None,
    location_counties_served: str | None = None,
    location_population_served: str | None = None,
    location_sampler_type: str | None = None,
    location_collection_basis: str | None = None,
    location_collection_type: str | None = None,
    location_zipcode: str | None = None,
    location_comment: str | None = None,
    location_lng: float | None = None,
    location_lat: float | None = None,
    min_collection_window_hrs: float | None = None,
    max_collection_window_hrs: float | None = None,
    min_collection_pull_ml: float | None = None,
    max_collection_pull_ml: float | None = None,
    min_collection_step_min: float | None = None,
    max_collection_step_min: float | None = None,
    skip: int = 0,
    limit: int = 10000,
):
    filters = []
    
    # ---- String / categorical filters ----
    if location_id is not None:
        filters.append(func.lower(func.trim(Location.location_id)) == location_id.strip().lower())

    if sample_code_prefix is not None:
        filters.append(func.lower(func.trim(Location.sample_code_prefix)) == sample_code_prefix.strip().lower())

    if location_primary_lab is not None:
        filters.append(func.lower(func.trim(Location.location_primary_lab)) == location_primary_lab.strip().lower())

    if location_status is not None:
        filters.append(func.lower(func.trim(Location.location_status)) == location_status.strip().lower())

    if location_common_name is not None:
        filters.append(func.lower(func.trim(Location.location_common_name)) == location_common_name.strip().lower())

    if location_category is not None:
        filters.append(func.lower(func.trim(Location.location_category)) == location_category.strip().lower())

    if location_group is not None:
        filters.append(func.lower(func.trim(Location.location_group)) == location_group.strip().lower())

    if location_primary_wwtp_id is not None:
        filters.append(func.lower(func.trim(Location.location_primary_wwtp_id)) == location_primary_wwtp_id.strip().lower())

    if location_counties_served is not None:
        filters.append(func.lower(func.trim(Location.location_counties_served)) == location_counties_served.strip().lower())

    if location_population_served is not None:
        filters.append(func.lower(func.trim(Location.location_population_served)) == location_population_served.strip().lower())

    if location_sampler_type is not None:
        filters.append(func.lower(func.trim(Location.location_sampler_type)) == location_sampler_type.strip().lower())

    if location_collection_basis is not None:
        filters.append(func.lower(func.trim(Location.location_collection_basis)) == location_collection_basis.strip().lower())

    if location_collection_type is not None:
        filters.append(func.lower(func.trim(Location.location_collection_type)) == location_collection_type.strip().lower())

    if location_zipcode is not None:
        filters.append(func.lower(func.trim(Location.location_zipcode)) == location_zipcode.strip().lower())

    if location_comment is not None:
        filters.append(func.lower(func.trim(Location.location_comment)) == location_comment.strip().lower())

    # ---- Numeric filters ----
    if location_lng is not None:
        filters.append(Location.location_lng == location_lng)

    if location_lat is not None:
        filters.append(Location.location_lat == location_lat)

    if min_collection_window_hrs is not None:
        filters.append(Location.location_collection_window_hrs >= min_collection_window_hrs)

    if max_collection_window_hrs is not None:
        filters.append(Location.location_collection_window_hrs <= max_collection_window_hrs)

    if min_collection_pull_ml is not None:
        filters.append(Location.location_collection_pull_ml >= min_collection_pull_ml)

    if max_collection_pull_ml is not None:
        filters.append(Location.location_collection_pull_ml <= max_collection_pull_ml)

    if min_collection_step_min is not None:
        filters.append(Location.location_collection_step_min >= min_collection_step_min)

    if max_collection_step_min is not None:
        filters.append(Location.location_collection_step_min <= max_collection_step_min)

    # Build statement
    stmt = select(Location).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result