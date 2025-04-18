TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)", "Human Norovirus GII (HuNoV-GII)")
#DISEASES <- c("COVID", "FLUA", "FLUB", "RSV", "Norovirus GII", "Norovirus GI")
DISEASES <- c("COVID", "FLUA", "FLUB", "RSV", "NoVII")

# Do we need this?
GENLOCI <- c("SC2", "M", "NEP/NS1", "N2", "G")

TRENDL_03_COLOR <- "#00B140"	# MU Green
TRENDL_12_COLOR <- "#EAAA00"	# WVU Gold
#TRENDL_12_COLOR <- "#002855"	# WVU Blue

#002855	# WVU Blue

VIEW_RANGES <- c(1, 3, 6, 12, 24)
DATE_BREAKS <- c("5 days", "15 days", "1 month", "2 months", "4 months")
DATE_LABELS <- c("%d-%b", "%d-%b", "%b '%y", "%b '%y", "%b '%y")

VIEW_RANGE_PRIMARY <- 6

MAP_CENTER <- list2env(list(lat = 38.95, lng = -80.2, zoom = 7))

STALE_DATA_THRESHOLDS <- c(14, 21, 30)
STALE_DATA_COLORS <- c("#000000", "#444444", "#808080", "#d0d0d0")

ALERT_LEVEL_THRESHOLDS <- c(0, 100, 150)
ALERT_LEVEL_COLORS <- c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE", "#FF85FF", "#96F786")
ALERT_LEVEL_STRINGS <- c("LOW", "MODERATE", "HIGH", "VERY HIGH", "UNKNOWN")

TREND_STRINGS <- c("DECREASING", "STABLE", "VARIABLE", "INCREASING", "INDETERMINATE", "SPIKING", "DESPIKING")

#Brewer Spectral Palette
# [1] "#9E0142" "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#FFFFBF" "#E6F598" "#ABDDA4"
# [9] "#66C2A5" "#3288BD" "#5E4FA2"

MAP_COLORS <- c(c("COVID"), c("FLUA", "FLUB"), c("RSV"))
MAP_COLORS_DEFAULT <- "COVID"

abundance_level_colors = c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE")
names(abundance_level_colors) <- c("LOW", "MODERATE", "HIGH", "VERY HIGH", "UNKNOWN")

trend_level_colors = c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE", "#FF85FF", "#96F786")
names(trend_level_colors) <- c("DECREASING", "STABLE", "VARIABLE", "INCREASING", "INDETERMINATE", "SPIKING", "DESPIKING")

lab_colors = c("#EAAA00", "#00B140")
names(lab_colors) <- c("ZooWVU", "MUIDSL")

default_colors = c("#ff0000")
