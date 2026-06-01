# app/schemas/samples.py

from pydantic import BaseModel
from datetime import date, datetime
from typing import Optional

class SamplesSchema(BaseModel):
    sample_id: str
    sample_status: Optional[str]
    location_id: Optional[str]
    sample_event: Optional[str]
    sample_qc: Optional[str]
    sample_collection_start_datetime: Optional[datetime]
    sample_collection_end_datetime: Optional[datetime]
    sample_recovered_datetime: Optional[datetime]
    sample_flow: Optional[float]
    sample_received_by: Optional[str]
    sample_received_date: Optional[date]
    sample_ph_lab: Optional[float]
    sample_comment: Optional[str]
    
    class Config:
        from_attributes = True