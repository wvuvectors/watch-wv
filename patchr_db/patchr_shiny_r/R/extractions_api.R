# R/extractions_api.R

get_extractions <- function(limit = 10000) {
  
  fetch_api(
    "/extractions/",
    params = list(limit = limit)
  )
}

query_extractions_api <- function(
    extraction_id = NULL,
    concentration_id = NULL,
    extraction_batch_id = NULL,
    extraction_location_in_batch = NULL,
    extraction_location_in_storage = NULL
) {
  
  params <- list()
  
  if (!is.null(extraction_id) && extraction_id != "")
    params$extraction_id <- extraction_id
  
  if (!is.null(concentration_id) && concentration_id != "")
    params$concentration_id <- concentration_id
  
  if (!is.null(extraction_batch_id) && extraction_batch_id != "")
    params$extraction_batch_id <- extraction_batch_id
  
  if (!is.null(extraction_location_in_batch) && extraction_location_in_batch != "")
    params$extraction_location_in_batch <- extraction_location_in_batch
  
  if (!is.null(extraction_location_in_storage) && extraction_location_in_storage != "")
    params$extraction_location_in_storage && extraction_location_in_storage
  
  fetch_api(
    "/extractions/query",
    params = params
  )
}