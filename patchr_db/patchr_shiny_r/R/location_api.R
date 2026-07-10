# R/location_api.R

get_locations <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/location/",
    params = list(limit = limit)
  )
}


query_locations_api <- function(
    limit = as.integer(100000),
    location_id = NULL,
    location_common_name = NULL,
    location_primary_wwtp_id = NULL,
    location_zipcode = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(location_id) && location_id != "")
    params$location_id <- location_id
  
  if (!is.null(location_common_name) && location_common_name != "")
    params$location_common_name <- location_common_name
  
  if (!is.null(location_primary_wwtp_id) && location_primary_wwtp_id != "")
    params$location_primary_wwtp_id <- location_primary_wwtp_id
  
  if (!is.null(location_zipcode) && location_zipcode != "")
    params$location_zipcode <- location_zipcode
  
  fetch_api(
    "/location/query",
    params = params
  )
}