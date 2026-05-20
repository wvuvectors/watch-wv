# R/concentration_api.R

get_concentration <- function(limit = 10000) {
  
  fetch_api(
    "/concentration/",
    params = list(limit = limit)
  )
}

query_concentration_api <- function(
    concentration_id = NULL,
    concentration_batch_id = NULL,
    sample_id = NULL
) {
  
  params <- list()
  
  if (!is.null(concentration_id) && concentration_id != "")
    params$concentration_id <- concentration_id
  
  if (!is.null(concentration_batch_id) &&
      concentration_batch_id != "")
    params$concentration_batch_id <- concentration_batch_id
  
  if (!is.null(sample_id) && sample_id != "")
    params$sample_id <- sample_id
  
  fetch_api(
    "/concentration/query",
    params = params
  )
}