GEOLEVELS <- c("Facility", "County")
GEOLEVELS_DEFAULT <- "County"

#TARGETS_RS <- c("Influenza Virus A (FluA)", "SARS-CoV-2", "SARS-CoV-2", "Respiratory Syncitial Virus, Human (RSV)")
#GENLOCI_RS <- c("M", "N2", "SC2", "G")
#DISEASE_RS <- c("FLUA", "COVID", "COVID", "RSV")

TARGETS_RS <- c("Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "SARS-CoV-2", "Respiratory Syncitial Virus, Human (RSV)")
DISEASE_RS <- c("FLUA", "FLUB", "COVID", "RSV")

# Do we need this?
GENLOCI_RS <- c("M", "NEP/NS1", "N2", "SC2", "G")

# TARGETS_RS <- c("Influenza Virus A (FluA)", "Human Norovirus GII (HuNoV-GII)", "SARS-CoV-2", "Respiratory Syncitial Virus, Human (RSV)")
# GENLOCI_RS <- c("M", "ORF1_2", "N2:SARS", "G")
# DISEASE_RS <- c("FLUA", "NoV", "COVID", "RSV")

TRENDL_MO_COLOR <- "#00B140"	# MU Green
TRENDL_YR_COLOR <- "#EAAA00"	# WVU Gold

#002855	# WVU Blue

VIEW_RANGES <- c(1, 3, 6, 12, 24)
DATE_BREAKS <- c("5 days", "15 days", "1 month", "2 months", "4 months")
DATE_LABELS <- c("%d-%b", "%d-%b", "%b '%y", "%b '%y", "%b '%y")

VIEW_RANGE_PRIMARY <- 6

MAP_CENTER <- list2env(list(lat = 38.95, lng = -80.2, zoom = 7))

STALE_DATA_THRESHOLDS <- c(14, 21, 30)
STALE_DATA_COLORS <- c("#000000", "#444444", "#808080", "#d0d0d0")

ALERT_LEVEL_THRESHOLDS <- c(50, 100, 150)
ALERT_LEVEL_COLORS <- c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE")


#Brewer Spectral Palette
# [1] "#9E0142" "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#FFFFBF" "#E6F598" "#ABDDA4"
# [9] "#66C2A5" "#3288BD" "#5E4FA2"
