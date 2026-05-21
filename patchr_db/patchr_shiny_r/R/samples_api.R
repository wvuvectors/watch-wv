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
  recovered_start = NULL,
  recovered_end = NULL
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
  
  if (!is.null(recovered_start) && recovered_start != "")
    params$recovered_start <- recovered_start
  
  if (!is.null(recovered_end) && recovered_end != "")
    params$recovered_end <- recovered_end

  fetch_api(
    "/samples/query",
    params = params
  )
}