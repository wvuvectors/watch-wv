# app/crud/county.py
# Defines helper functions to be used throughout app 

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_

from app.models.county import County

# Function to retrieve a single county by its county_id
def get_county_by_id(db: Session, county_id: str):
    return (
        db.query(County)
        .filter(County.county_id == county_id)
        .first()
    )

# List counties
def list_counties(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(County)
        .offset(skip)
        .limit(limit)
        .all()
    )

# Dynamic Query
def query_counties(
    db: Session,
    *,
    county_id: str | None = None,
    county_labcode: str | None = None,
    county_fips: str | None = None,
    county_name: str | None = None,
    county_population: float | None = None,
    min_population: float | None = None,
    max_population: float | None = None,
    skip: int = 0,
    limit: int = 10000,
):
    filters = []
    
    # ---- String / categorical filters ----
    if county_id is not None:
        filters.append(func.lower(func.trim(County.county_id)) == county_id.strip().lower())

    if county_labcode is not None:
        filters.append(func.lower(func.trim(County.county_labcode)) == county_labcode.strip().lower())

    if county_fips is not None:
        filters.append(func.lower(func.trim(County.county_fips)) == county_fips.strip().lower())

    if county_name is not None:
        filters.append(func.lower(func.trim(County.county_name)) == county_name.strip().lower())

    # ---- Numeric filters ----
    if county_population is not None:
        filters.append(County.county_population == county_population)

    if min_population is not None:
        filters.append(County.county_population >= min_population)

    if max_population is not None:
        filters.append(County.county_population <= max_population)

    # Build statement
    stmt = select(County).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result