# R/concentration_api.R

get_concentration <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/concentration/",
    params = list(limit = limit)
  )
}

query_concentration_api <- function(
    limit = as.integer(100000),
    concentration_id = NULL,
    concentration_batch_id = NULL,
    sample_id = NULL,
    concentration_location_in_batch = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(concentration_id) && concentration_id != "")
    params$concentration_id <- concentration_id
  
  if (!is.null(concentration_batch_id) && concentration_batch_id != "")
    params$concentration_batch_id <- concentration_batch_id
  
  if (!is.null(sample_id) && sample_id != "")
    params$sample_id <- sample_id
  
  if (!is.null(concentration_location_in_batch) && concentration_location_in_batch != "")
    params$concentration_location_in_batch <- concentration_location_in_batch
  
  fetch_api(
    "/concentration/query",
    params = params
  )
}