# app/schemas/extractions.py

from pydantic import BaseModel
from typing import Optional

class ExtractionsSchema(BaseModel):
    extraction_id: str
    concentration_id: str | None = None
    extraction_batch_id: str | None = None
    extraction_location_in_batch: str | None = None
    extraction_location_in_storage: str | None = None
    extraction_comment: Optional[str]
    
    class Config:
        from_attributes = True