# R/ebatch_api.R

get_ebatch <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/ebatch/",
    params = list(limit = limit)
  )
}

query_ebatch_api <- function(
    limit = as.integer(100000),
    extraction_batch_id = NULL,
    extraction_date = NULL,
    extraction_method
) {
  
  params <- list(limit = limit)
  
  if (!is.null(extraction_batch_id) && extraction_batch_id != "")
    params$extraction_batch_id <- extraction_batch_id
  
  if (!is.null(extraction_date) && !is.na(extraction_date))
    params$extraction_date <- as.character(extraction_date)
  
  if (!is.null(extraction_method) && extraction_method != "")
    params$extraction_method <- extraction_method
  
  fetch_api(
    "/ebatch/query",
    params = params
  )
}