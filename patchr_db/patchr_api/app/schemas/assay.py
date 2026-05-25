# app/schemas/assay.py

from pydantic import BaseModel
from typing import Optional

class AssaySchema(BaseModel):
    assay_id: str
    extraction_id: str
    sample_id: str
    assay_batch_id: str | None = None
    assay_location_in_batch: str | None = None
    assay_input_ul: Optional[float]
    assay_class: str | None = None
    assay_type: str | None = None
    assay_target: str | None = None
    assay_target_genetic_locus: str | None = None
    assay_template: str | None = None
    assay_target_marcomolecule: str | None = None
    assay_target_flourophore: str | None = None
    assay_accepted_droplets: Optional[float]
    assay_target_predicted_copies_per_ul_reaction: Optional[float]
    assay_target_copies_per_ul_reaction: Optional[float]
    assay_comment: Optional[str]
    
    class Config:
        from_attributes = True