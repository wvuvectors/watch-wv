# app/schemas/concentration.py

from pydantic import BaseModel
from typing import Optional

class ConcentrationSchema(BaseModel):
    concentration_id: str
    sample_id: str
    concentration_batch_id: Optional[str]
    concentration_location_in_batch: Optional[str]
    concentration_comment: Optional[str]
    
    class Config:
        from_attributes = True