# R/api_helpers.R

library(httr)
library(jsonlite)

BASE_API_URL <- "http://127.0.0.1:8000"

# ----------------------------
# Generic API Fetch Helper
# ----------------------------
fetch_api <- function(endpoint, params = list()) {
  
  url <- paste0(BASE_API_URL, endpoint)
  
  response <- GET(url, query = params)
  
  if (status_code(response) != 200) {
    warning(
      paste(
        "API request failed:",
        status_code(response)
      )
    )
    
    return(data.frame())
  }
  
  data <- fromJSON(
    content(response, "text", encoding = "UTF-8")
  )
  
  return(as.data.frame(data))
}