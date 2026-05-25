# R/assay.R

get_assays <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/assay/",
    params = list(limit = limit)
  )
}

query_assays_api <- function(
    limit = as.integer(100000),
    assay_id = NULL,
    extraction_id = NULL,
    sample_id = NULL,
    assay_batch_id = NULL,
    assay_location_in_batch = NULL,
    assay_input_ul = NULL,
    assay_class = NULL,
    assay_type = NULL,
    assay_target = NULL,
    assay_target_genetic_locus = NULL,
    assay_template = NULL,
    assay_target_macromolecule = NULL,
    assay_target_fluorophore = NULL,
    assay_accepted_droplets = NULL,
    assay_target_predicted_copies_per_ul_reaction = NULL,
    assay_target_copies_per_ul_reaction = NULL,
    assay_comment = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(assay_id) && assay_id != "")
    params$assay_id <- assay_id
  
  if (!is.null(extraction_id) && extraction_id != "")
    params$extraction_id <- extraction_id
  
  if (!is.null(sample_id) && sample_id != "")
    params$sample_id <- sample_id
  
  if (!is.null(assay_batch_id) && assay_batch_id != "")
    params$assay_batch_id <- assay_batch_id
  
  if (!is.null(assay_location_in_batch) && assay_location_in_batch != "")
    params$assay_location_in_batch <- assay_location_in_batch
  
  if (!is.null(assay_input_ul))
    params$assay_input_ul <- as.numeric(assay_input_ul)
  
  if (!is.null(assay_class) && assay_class != "")
    params$assay_class <- assay_class
  
  if (!is.null(assay_type) && assay_type != "")
    params$assay_type <- assay_type
  
  if (!is.null(assay_target) && assay_target != "")
    params$assay_target <- assay_target
  
  if (!is.null(assay_target_genetic_locus) && assay_target_genetic_locus != "")
    params$assay_target_genetic_locus <- assay_target_genetic_locus
  
  if (!is.null(assay_template) && assay_template != "")
    params$assay_template <- assay_template
  
  if (!is.null(assay_target_macromolecule) && assay_target_macromolecule != "")
    params$assay_target_macromolecule <- assay_target_macromolecule
  
  if (!is.null(assay_target_fluorophore) && assay_target_fluorophore != "")
    params$assay_target_fluorophore <- assay_target_fluorophore
  
  if (!is.null(assay_accepted_droplets))
    params$assay_accepted_droplets <- as.numeric(assay_accepted_droplets)
  
  if (!is.null(assay_target_predicted_copies_per_ul_reaction))
    params$assay_target_predicted_copies_per_ul_reaction <- as.numeric(assay_target_predicted_copies_per_ul_reaction)
  
  if (!is.null(assay_target_copies_per_ul_reaction))
    params$assay_target_copies_per_ul_reaction <- as.numeric(assay_target_copies_per_ul_reaction)
  
  if (!is.null(assay_comment) && assay_comment != "")
    params$assay_comment <- assay_comment
  
  fetch_api(
    "/assay/query",
    params = params
  )
}