# R/samples_api.R

get_samples <- function(limit = 10000) {

  fetch_api(
    "/samples/",
    params = list(limit = limit)
  )
}

query_samples_api <- function(
  sample_id = NULL
) {

  params <- list()

  if (!is.null(sample_id) && sample_id != "")
    params$sample_id <- sample_id

  fetch_api(
    "/samples/query",
    params = params
  )
}