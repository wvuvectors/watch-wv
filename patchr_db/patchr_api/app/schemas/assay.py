# app/schemas/assay.py

from pydantic import BaseModel
from typing import Optional

class AssaySchema(BaseModel):
    assay_id: str
    extraction_id: str
    sample_id: str
    assay_batch_id: Optional[str]
    assay_location_in_batch: Optional[str]
    assay_input_ul: Optional[float]
    assay_class: Optional[str]
    assay_type: Optional[str]
    assay_target: Optional[str]
    assay_target_genetic_locus: Optional[str]
    assay_template: Optional[str]
    assay_target_macromolecule: Optional[str]
    assay_target_fluorophore: Optional[str]
    assay_accepted_droplets: Optional[float]
    assay_target_predicted_copies_per_ul_reaction: Optional[float]
    assay_target_copies_per_ul_reaction: Optional[float]
    assay_comment: Optional[str]
    
    class Config:
        from_attributes = True