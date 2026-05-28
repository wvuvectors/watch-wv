# app/schemas/cbatch.py

from pydantic import BaseModel
from datetime import date, datetime
from typing import Optional

class cBatchSchema(BaseModel):
    concentration_batch_id: str
    concentration_date: Optional[date]
    concentration_input_ml: Optional[float]
    concentration_machine: Optional[str]
    concentration_method: Optional[str]
    concentration_method_lot_id: Optional[str]
    concentration_output_ml: Optional[float]
    concentration_run_by: Optional[str]
    concentration_batch_record_version: Optional[str]
    concentration_batch_comment: Optional[str]
    
    class Config: 
        from_attributes = True