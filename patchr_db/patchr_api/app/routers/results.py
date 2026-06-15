# app/routers/results.py

from datetime import date, datetime

from fastapi import APIRouter, Depends, Query
from sqlalchemy.orm import Session

from app.database import get_db
from app.schemas.results import ResultsSchema
from app.crud.results import query_results

# Create a router object
router = APIRouter(
    prefix="/results",
    tags=["results"]
)

# Results endpoint
@router.get(
    "/query",
    response_model=list[ResultsSchema]
)
def query_results_endpoint(
    location_id: str | None = Query(None),

    recovered_start: date | None = Query(None),
    recovered_end: date | None = Query(None),
    
    assay_target: str | None = Query(None),

    skip: int = Query(0),
    limit: int = Query(100),

    db: Session = Depends(get_db)
):

    # Convert dates to datetimes 
    recovered_start_dt = (
        datetime.combine(
            recovered_start,
            datetime.min.time()
        )
        if recovered_start
        else None
    )

    recovered_end_dt = (
        datetime.combine(
            recovered_end,
            datetime.max.time()
        )
        if recovered_end
        else None
    )

    return query_results(
        db=db,
        location_id=location_id,
        recovered_start=recovered_start_dt,
        recovered_end=recovered_end_dt,
        assay_target=assay_target,
        skip=skip,
        limit=limit
    )