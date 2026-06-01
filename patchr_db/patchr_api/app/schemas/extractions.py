# app/schemas/extractions.py

from pydantic import BaseModel
from typing import Optional

class ExtractionsSchema(BaseModel):
    extraction_id: str
    concentration_id: str
    extraction_batch_id: Optional[str]
    extraction_location_in_batch: Optional[str]
    extraction_location_in_storage: Optional[str]
    extraction_comment: Optional[str]
    
    class Config:
        from_attributes = True