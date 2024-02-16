GEOLEVELS <- c("Facility", "County")
GEOLEVELS_DEFAULT <- "County"

TARGETS_RS <- c("Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "SARS-CoV-2", "Respiratory Syncitial Virus, Human (RSV)")
GENLOCI_RS <- c("M", "NEP/NS1", "N2:SARS", "G")
DISEASE_RS <- c("FLUA", "FLUB", "COVID", "RSV A/B")


VIEW_RANGES <- c(1, 3, 6, 12, 24)
VIEW_RANGE_PRIMARY <- 6

MAP_CENTER <- list2env(list(lat = 39.35, lng = -79.4, zoom = 7))

STALE_DATA_THRESHOLDS <- c(10, 14, 21)
STALE_DATA_COLORS <- c("#000000", "#555555", "#999999", "#aaaaaa")

ALERT_LEVEL_COLORS <- c("#85BEFA", "#FBF4A2", "#E47E47", "#DB4F40")
