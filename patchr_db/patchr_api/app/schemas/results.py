# app/schemas/results.py

from pydantic import BaseModel
from datetime import datetime

class ResultsSchema(BaseModel):
    sample_id: str
    location_id: str | None = None
    sample_collection_datetime: datetime | None = None
    assay_target: str | None = None
    assay_target_genetic_locus: str | None = None
    assay_target_copies_per_ul_reaction: float | None = None
    concentration_input_ml: float | None = None
    concentration_output_ml: float | None = None
    extraction_input_ul: float | None = None
    extraction_output_ul: float | None = None
    assay_input_ul: float | None = None
    assay_reaction_ul: float | None = None
    copies_per_l_wastewater: float | None = None