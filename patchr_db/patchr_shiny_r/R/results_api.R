# R/results_api.R

get_results <- function(limit = as.integer(100)) {
  
  fetch_api(
    "/results/",
    params = list(limit = limit)
  )
}

get_results_location_ids <- function() {
  
  data <- get_samples(limit = as.integer(100000))
  
  sort(unique(data$location_id))
}

get_results_assay_targets <- function() {
  
  data <- get_assays(limit = as.integer(100000))
  
  sort(unique(data$assay_target))
}

get_results_genetic_loci <- function() {
  
  data <- get_assays(limit = as.integer(100000))
  
  sort(unique(data$assay_target_genetic_locus))
}

query_results_api <- function(
    limit = as.integer(100000),
    location_id = NULL,
    collection_start = NULL,
    collection_end = NULL,
    assay_target = NULL,
    assay_target_genetic_locus = NULL
) {
  
  params <- list(limit = limit)
  
  if (!is.null(location_id) && location_id != "")
    params$location_id <- location_id
  
  if (!is.null(collection_start) && !is.na(collection_start))
    params$collection_start <- as.character(collection_start)
  
  if (!is.null(collection_end) && !is.na(collection_end))
    params$collection_end <- as.character(collection_end)
  
  if (!is.null(assay_target) && assay_target != "")
    params$assay_target <- assay_target
  
  if (!is.null(assay_target_genetic_locus) &&
      assay_target_genetic_locus != "")
    params$assay_target_genetic_locus <- assay_target_genetic_locus
  
  results <- fetch_api(
    "/results/query",
    params = params
  )
  
  # Move copies_per_l_wastewater to second column
  if (
    nrow(results) > 0 &&
    "copies_per_l_wastewater" %in% names(results)
  ) {
    
    desired_order <- c(
      names(results)[1],
      "copies_per_l_wastewater",
      setdiff(
        names(results)[-1],
        "copies_per_l_wastewater"
      )
    )
    
    results <- results[, desired_order, drop = FALSE]
  }
  
  return(results)
}