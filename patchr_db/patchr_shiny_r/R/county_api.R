# R/county_api.R

get_counties <- function(limit = as.integer(100000)) {
  
  fetch_api(
    "/county/",
    params = list(limit = limit)
  )
}


query_counties_api <- function(
    limit = as.integer(100000),
    county_id = NULL,
    county_labcode = NULL,
    county_name = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(county_id) && county_id != "")
    params$county_id <- county_id
  
  if (!is.null(county_labcode) && county_labcode != "")
    params$county_labcode <- county_labcode
  
  if (!is.null(county_name) && county_name != "")
    params$county_name <- county_name
  
  fetch_api(
    "/county/query",
    params = params
  )
}