# app/schemas/concentration.py

from pydantic import BaseModel
from typing import Optional

class ConcentrationSchema(BaseModel):
    concentration_id: str
    sample_id: str | None = None
    concentration_batch_id: str | None = None
    concentration_location_in_batch: str | None = None
    concentration_comment: Optional[str]
    
    class Config:
        from_attributes = True