# app/schemas/abatch.py

from pydantic import BaseModel
from datetime import date, datetime
from typing import Optional

class aBatchSchema(BaseModel):
    assay_batch_id: str
    assay_date: Optional[date]
    assay_reaction_ul: Optional[float]
    assay_machine: Optional[str]
    assay_amplification_method: Optional[str]
    assay_amplification_method_lot_id: Optional[str]
    assay_quantification_method: Optional[str]
    assay_quantification_type: Optional[str]
    assay_batch_record_version: Optional[str]
    assay_method: Optional[str]
    assay_method_lot_id: Optional[str]
    assay_qx_manager_version: Optional[str]
    assay_run_by: Optional[str]
    assay_batch_comment: Optional[str]
    
    class Config: 
        from_attributes = True