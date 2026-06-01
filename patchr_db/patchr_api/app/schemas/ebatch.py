# app/schemas/ebatch.py

from pydantic import BaseModel
from datetime import date, datetime
from typing import Optional

class eBatchSchema(BaseModel):
    extraction_batch_id: str
    extraction_date: Optional[date]
    extraction_input_ul: Optional[float]
    extraction_eluant: Optional[str]
    extraction_machine: Optional[str]
    extraction_method: Optional[str]
    extraction_method_lot_id: Optional[str]
    extraction_output_ul: Optional[float]
    extraction_batch_record_version: Optional[str]
    extraction_run_by: Optional[str]
    extraction_batch_comment: Optional[str]
    
    class Config: 
        from_attributes = True