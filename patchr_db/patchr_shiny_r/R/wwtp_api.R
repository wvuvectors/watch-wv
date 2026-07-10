# R/wwtp_api.R

get_wwtps <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/wwtp/",
    params = list(limit = limit)
  )
}


query_wwtps_api <- function(
    limit = as.integer(100000),
    wwtp_id = NULL,
    wwtp_site_id = NULL,
    wwtp_common_name = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(wwtp_id) && wwtp_id != "")
    params$wwtp_id <- wwtp_id
  
  if (!is.null(wwtp_site_id) && wwtp_site_id != "")
    params$wwtp_site_id <- wwtp_site_id
  
  if (!is.null(wwtp_common_name) && wwtp_common_name != "")
    params$wwtp_common_name <- wwtp_common_name
  
  fetch_api(
    "/wwtp/query",
    params = params
  )
}