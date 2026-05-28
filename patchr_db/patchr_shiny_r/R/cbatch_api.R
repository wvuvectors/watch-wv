# R/cbatch_api.R

get_cbatch <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/cbatch/",
    params = list(limit = limit)
  )
}

query_cbatch_api <- function(
    limit = as.integer(100000),
    concentration_batch_id = NULL,
    concentration_date = NULL,
    concentration_method = NULL
) {
  if (!is.null(concentration_batch_id) && concentration_batch_id != "")
    params$concentration_batch_id <- concentration_batch_id
  
  if (!is.null(concentration_date) && !is.na(concentration_date))
    params$concentration_date <- as.character(concentration_date)
  
  if (!is.null(concentration_method) && concentration_method)
    params$concentration_method <- concentration_method
  
  fetch_api(
    "/cbatch/query", 
    params = params
  )
}