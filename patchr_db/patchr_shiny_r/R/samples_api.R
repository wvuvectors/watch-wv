# R/samples_api.R

get_samples <- function(limit = as.integer(100000)) {

  fetch_api(
    "/samples/",
    params = list(limit = limit)
  )
}

query_samples_api <- function(
  limit = as.integer(100000),
  sample_id = NULL,
  sample_status = NULL,
  location_id = NULL,
  sample_qc = NULL,
  #recovered_start = NULL,
  #recovered_end = NULL,
  collection_start = NULL,
  collection_end = NULL
) {

  params <- list(limit = limit)

  if (!is.null(sample_id) && sample_id != "")
    params$sample_id <- sample_id
  
  if (!is.null(sample_status) && sample_status != "")
    params$sample_status <- sample_status
  
  if (!is.null(location_id) && location_id != "")
    params$location_id <- location_id
  
  if (!is.null(sample_qc) && sample_qc != "")
    params$sample_qc <- sample_qc
  
  #if (!is.null(recovered_start) && !is.na(recovered_start))
  #  params$recovered_start <- as.character(recovered_start)
  
  #if (!is.null(recovered_end) && !is.na(recovered_end))
  #  params$recovered_end <- as.character(recovered_end)
  
  if (!is.null(collection_start) && !is.na(collection_start))
    params$collection_start <- as.character(collection_start)
  
  if (!is.null(collection_end) && !is.na(collection_end))
    params$collection_end <- as.character(collection_end)

  fetch_api(
    "/samples/query",
    params = params
  )
}