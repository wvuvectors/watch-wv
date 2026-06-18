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

def get_location_ids(db: Session):
    return (
        db.query(Samples.location_id)
        .distinct()
        .order_by(Samples.location_id)
        .all()
    )
    
def get_assay_targets(db: Session):
    return (
        db.query(Assay.assay_target)
        .distinct()
        .order_by(Assay.assay_target)
        .all()
    )
    
def get_genetic_loci(db: Session):
    return (
        db.query(Assay.assay_target_genetic_locus)
        .distinct()
        .order_by(Assay.assay_target_genetic_locus)
        .all()
    )
    
def list_results(db: Session, skip: int = 0, limit: int = 100):
    return (
        db.query(Samples)
        .order_by(Samples.sample_recovered_datetime.desc())
        .offset(skip)
        .limit(limit)
        .all()
    )

def query_results(
    db: Session,
    *,
    location_id: str | None = None,
    recovered_start: datetime | None = None,
    recovered_end: datetime | None = None,
    assay_target: str | None = None,
    assay_target_genetic_locus: str | None = None,
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
        filters.append(Samples.location_id == location_id)

    if recovered_start:
        filters.append(Samples.sample_recovered_datetime >= recovered_start)

    if recovered_end:
        filters.append(Samples.sample_recovered_datetime <= recovered_end)
        
    if assay_target:
        filters.append(Assay.assay_target == assay_target)
        
    if assay_target_genetic_locus:
        filters.append(Assay.assay_target_genetic_locus == assay_target_genetic_locus)

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
                "assay_target_genetic_locus": assay.assay_target_genetic_locus,

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