# app/crud/wwtp.py
# Defines helper functions to be used throughout app 

from sqlalchemy.orm import Session
from sqlalchemy import select, func, and_

from app.models.wwtp import WWTP

# Function to retrieve a single WWTP by its wwtp_id
def get_wwtp_by_id(db: Session, wwtp_id: str):
    return (
        db.query(WWTP)
        .filter(WWTP.wwtp_id == wwtp_id)
        .first()
    )

# List WWTPs
def list_wwtps(db: Session, skip: int = 0, limit: int = 100000):
    return (
        db.query(WWTP)
        .offset(skip)
        .limit(limit)
        .all()
    )

# Dynamic Query
def query_wwtps(
    db: Session,
    *,
    wwtp_id: str | None = None,
    wwtp_site_id: str | None = None,
    wwtp_common_name: str | None = None,
    wwtp_authority_name: str | None = None,
    wwtp_counties_served: str | None = None,
    wwtp_epaid_id: str | None = None,
    wwtp_cwns_id: str | None = None,
    wwtp_capacity_mgd: float | None = None,
    wwtp_population_served: float | None = None,
    min_capacity_mgd: float | None = None,
    max_capacity_mgd: float | None = None,
    min_population_served: float | None = None,
    max_population_served: float | None = None,
    skip: int = 0,
    limit: int = 10000,
):
    filters = []
    
    # ---- String / categorical filters ----
    if wwtp_id is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_id)) == wwtp_id.strip().lower())

    if wwtp_site_id is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_site_id)) == wwtp_site_id.strip().lower())

    if wwtp_common_name is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_common_name)) == wwtp_common_name.strip().lower())

    if wwtp_authority_name is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_authority_name)) == wwtp_authority_name.strip().lower())

    if wwtp_counties_served is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_counties_served)) == wwtp_counties_served.strip().lower())

    if wwtp_epaid_id is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_epaid_id)) == wwtp_epaid_id.strip().lower())

    if wwtp_cwns_id is not None:
        filters.append(func.lower(func.trim(WWTP.wwtp_cwns_id)) == wwtp_cwns_id.strip().lower())

    # ---- Numeric filters ----
    if wwtp_capacity_mgd is not None:
        filters.append(WWTP.wwtp_capacity_mgd == wwtp_capacity_mgd)

    if wwtp_population_served is not None:
        filters.append(WWTP.wwtp_population_served == wwtp_population_served)

    if min_capacity_mgd is not None:
        filters.append(WWTP.wwtp_capacity_mgd >= min_capacity_mgd)

    if max_capacity_mgd is not None:
        filters.append(WWTP.wwtp_capacity_mgd <= max_capacity_mgd)

    if min_population_served is not None:
        filters.append(WWTP.wwtp_population_served >= min_population_served)

    if max_population_served is not None:
        filters.append(WWTP.wwtp_population_served <= max_population_served)

    # Build statement
    stmt = select(WWTP).where(and_(*filters)).offset(skip).limit(limit)
    result = db.execute(stmt).scalars().all()
    return result