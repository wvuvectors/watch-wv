# app/crud/results.py

# app/crud/results.py

from datetime import datetime
from sqlalchemy.orm import Session

from app.models.samples import Samples
from app.models.concentration import Concentration
from app.models.cbatch import cBatch
from app.models.extractions import Extractions
from app.models.ebatch import eBatch
from app.models.assay import Assay
from app.models.abatch import aBatch

def query_results(
    db: Session,
    *,
    location_id: str | None = None,
    recovered_start: datetime | None = None,
    recovered_end: datetime | None = None,
    skip: int = 0,
    limit: int = 100,
):

    query = (
        db.query(
            Samples,
            Concentration,
            cBatch,
            Extractions,
            eBatch,
            Assay,
            aBatch
        )
        .join(
            Concentration,
            Samples.sample_id == Concentration.sample_id
        )
        .join(
            cBatch,
            Concentration.concentration_batch_id == cBatch.concentration_batch_id
        )
        .join(
            Extractions,
            Concentration.concentration_id == Extractions.concentration_id
        )
        .join(
            eBatch,
            Extractions.extraction_batch_id == eBatch.extraction_batch_id
        )
        .join(
            Assay,
            Extractions.extraction_id == Assay.extraction_id
        )
        .join(
            aBatch,
            Assay.assay_batch_id == aBatch.assay_batch_id
        )
    )
    
    filters = []

    if location_id:
        filters.append(
            Samples.location_id == location_id
        )

    if recovered_start:
        filters.append(
            Samples.sample_recovered_datetime >= recovered_start
        )

    if recovered_end:
        filters.append(
            Samples.sample_recovered_datetime <= recovered_end
        )
        
    if assay_target:
        filters.append(
            Assay.assay_target == assay_target

    if filters:
        query = query.filter(*filters)
        
    rows = (
        query
        .offset(skip)
        .limit(limit)
        .all()
    )
    
    results = []

    for (
        sample,
        concentration,
        cbatch,
        extraction,
        ebatch,
        assay,
        abatch
    ) in rows:
    
        required_values = [
            assay.assay_target_copies_per_ul_reaction,
            assay.assay_input_ul,
            abatch.assay_reaction_ul,
            ebatch.extraction_input_ul,
            ebatch.extraction_output_ul,
            cbatch.concentration_input_ml,
            cbatch.concentration_output_ml,
        ]
        
        if any(v is None for v in required_values):
            copies_per_l = None

        else:
            copies_per_l = (
                1_000_000
                * assay.assay_target_copies_per_ul_reaction
                * abatch.assay_reaction_ul
                * ebatch.extraction_output_ul
                * cbatch.concentration_output_ml
            ) / (
                assay.assay_input_ul
                * ebatch.extraction_input_ul
                * cbatch.concentration_input_ml
            )
            
        results.append(
            {
                "sample_id": sample.sample_id,
                "location_id": sample.location_id,
                "sample_recovered_datetime": sample.sample_recovered_datetime,
                "assay_target": assay.assay_target,

                "assay_target_copies_per_ul_reaction":
                    assay.assay_target_copies_per_ul_reaction,

                "concentration_input_ml":
                    cbatch.concentration_input_ml,

                "concentration_output_ml":
                    cbatch.concentration_output_ml,

                "extraction_input_ul":
                    ebatch.extraction_input_ul,

                "extraction_output_ul":
                    ebatch.extraction_output_ul,

                "assay_input_ul":
                    assay.assay_input_ul,

                "assay_reaction_ul":
                    abatch.assay_reaction_ul,

                "copies_per_l_wastewater":
                    copies_per_l,
            }
        )

    return results