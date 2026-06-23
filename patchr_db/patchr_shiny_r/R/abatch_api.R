# R/abatch_api.R

get_abatch <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/abatch/",
    params = list(limit = limit)
  )
}

query_abatch_api <- function(
    limit = as.integer(100000),
    assay_batch_id = NULL,
    assay_start_date = NULL,
    assay_end_date = NULL,
    assay_amplification_method = NULL,
    assay_quantification_method = NULL,
    assay_method = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(assay_batch_id) && assay_batch_id != "")
    params$assay_batch_id <- assay_batch_id
  
  if (!is.null(assay_start_date) && !is.na(assay_start_date))
    params$assay_start_date <- as.character(assay_start_date)
  
  if (!is.null(assay_end_date) && !is.na(assay_end_date))
    params$assay_end_date <- as.character(assay_end_date)
  
  if (!is.null(assay_amplification_method) && assay_amplification_method != "")
    params$assay_amplification_method <- assay_amplification_method
  
  if (!is.null(assay_quantification_method) && assay_quantification_method != "")
    params$assay_quantification_method <- assay_quantification_method
  
  if (!is.null(assay_method) && assay_method != "")
    params$assay_method <- assay_method
  
  fetch_api(
    "/abatch/query",
    params = params
  )
}