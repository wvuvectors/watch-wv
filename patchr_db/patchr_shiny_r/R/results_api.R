# R/results_api.R

get_results <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/results/",
    params = list(limit = limit)
  )
}

query_results_api <- function(
    limit = as.integer(100000),
    location_id = NULL,
    recovered_start = NULL,
    recovered_end = NULL,
    assay_target = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(location_id) && location_id != "") 
    params$location_id <- location_id
  
  if (!is.null(recovered_start) && !is.na(recovered_start))
    params$recovered_start <- as.character(recovered_start)
  
  if (!is.null(recovered_end) && !is.na(recovered_end))
    params$recovered_end <- as.character(recovered_end)
  
  if (!is.null(assay_target) && assay_target != "")
    params$assay_target = assay_target
  
  fetch_api(
    "/results/query",
    params = params
  )
}