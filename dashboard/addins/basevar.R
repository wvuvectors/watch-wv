TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)", "Human Norovirus GII (HuNoV-GII)", "Human Norovirus GI (HuNoV-GI)")
DISEASES <- c("COVID", "FLUA", "FLUB", "RSV", "Norovirus GII", "Norovirus GI")
#DISEASES <- c("COVID", "FLUA", "FLUB", "RSV", "NoVII", "NoVI")

# Do we need this?
GENLOCI <- c("SC2", "M", "NEP/NS1", "N2", "ORF1_ORF2", "ORF1_ORF2")

TRENDL_03_COLOR <- "#00B140"	# MU Green
TRENDL_12_COLOR <- "#EAAA00"	# WVU Gold
#TRENDL_12_COLOR <- "#002855"	# WVU Blue

#002855	# WVU Blue

VIEW_RANGES <- c(1, 3, 6, 12, 24)
DATE_BREAKS <- c("5 days", "15 days", "1 month", "2 months", "4 months")
DATE_LABELS <- c("%d-%b", "%d-%b", "%b '%y", "%b '%y", "%b '%y")

VIEW_RANGE_PRIMARY <- 12

MAP_CENTER <- list2env(list(lat = 38.95, lng = -80.2, zoom = 7))

STALE_DATA_THRESHOLDS <- c(14, 21, 30)
STALE_DATA_COLORS <- c("#000000", "#444444", "#808080", "#d0d0d0")

ALERT_LEVEL_THRESHOLDS <- c(0, 100, 150)
ALERT_LEVEL_COLORS <- c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE", "#FF85FF", "#96F786")
ALERT_LEVEL_STRINGS <- c("LOW", "MODERATE", "HIGH", "VERY HIGH", "UNKNOWN")
ALERT_LEVEL_DESCRIPTIONS <- c(
"The amount of this disease is less than half of the 3 month average (green dotted line on the graph), and community transmission is estimated to be low.",
"The amount of this disease is between 50% and 100% of the 3 month average (green dotted line on the graph), and community transmission is estimated to be moderate.",
"The amount of this disease is between 100% and 150% of the 3 month average (green dotted line on the graph), and community transmission is estimated to be high.",
"The amount of this disease is greater than 150% of the 3 month average (green dotted line on the graph), and community transmission is estimated to be very high.",
"The data for this disease agent is either missing or too old to make an accurate determination."
)

TREND_STRINGS <- c("DECREASING", "STABLE", "VARIABLE", "INCREASING", "INDETERMINATE", "SPIKING", "DESPIKING")
TREND_DESCRIPTIONS <- c(
"This disease has generally decreased in abundance over the last several weeks.",
"This disease has been relatively stable over the last several weeks.",
"The abundance of this disease has varied over the last several weeks.",
"Over the last several weeks, this disease has risen more than it has declined.",
"There is not enough data for this disease to make an accurate determination of how it may be trending.",
"The abundance of this disease is more than 500% greater than the previous value.",
"The abundance of this disease is less than 500% of the previous value."
)

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
