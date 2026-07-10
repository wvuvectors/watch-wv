# app/schemas/location.py

from pydantic import BaseModel
from typing import Optional

class LocationSchema(BaseModel):
    location_id: str
    sample_code_prefix: Optional[str]
    location_primary_lab: Optional[str]
    location_status: Optional[str]
    location_common_name: Optional[str]
    location_category: Optional[str]
    location_group: Optional[str]
    location_lng: Optional[float]
    location_lat: Optional[float]
    location_primary_wwtp_id: Optional[str]
    location_counties_served: Optional[str]
    location_population_served: Optional[str]
    location_sampler_type: Optional[str]
    location_collection_window_hrs: Optional[float]
    location_collection_pull_ml: Optional[float]
    location_collection_step_min: Optional[float]
    location_collection_basis: Optional[str]
    location_collection_type: Optional[str]
    location_zipcode: Optional[str]
    location_comment: Optional[str]

    class Config:
        from_attributes = True