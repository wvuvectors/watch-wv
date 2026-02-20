TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)")
DISEASES <- c("COVID", "FluA", "FluB", "RSV")
DISEASE_LABELS <- c("COVID", "Influenza A", "Influenza B", "RSV")

# Do we need this?
GENLOCI <- c("SC2", "M", "NEP/NS1", "N2")

TRENDL_03_COLOR <- "#00B140"	# MU Green
TRENDL_12_COLOR <- "#EAAA00"	# WVU Gold

# "#002855"	is WVU Blue
# "#CAE9FD" is a nice light blue

VIEW_RANGES <- c(1, 3, 6, 12, 24)
VIEW_RANGE_PRIMARY <- 12

DATE_BREAKS <- c("5 days", "2 weeks", "1 month", "1 month", "3 months")
DATE_BREAKS_MINOR <- c("1 day", "2 days", "1 week", "1 week", "2 weeks")
DATE_LABELS <- c("%d-%b", "%d-%b", "%b '%y", "%b '%y", "%b '%y")

MAP_CENTER <- list2env(list(lat = 38.95, lng = -80.2, zoom = 7))

# STALE_THRESHOLD_DAYS is used in global.R AND the isStale function in basefun.R.
# Make sure if you change the logic, you change it in BOTH places!
STALE_THRESHOLD_DAYS <- 21


ALEVEL_COLORS <- c("#DDDDDD", "#9289D6", "#3288BD", "#DEE998", "#FDAE61", "#D53E4F")
ALEVEL_STRINGS <- c("UNKNOWN", "VERY LOW", "LOW", "MODERATE", "HIGH", "VERY HIGH")
ALEVEL_DESCRIPTIONS <- c(
"is too old to make an accurate determination.", 
"EXTREMELY LOW. It is approaching the limit of detection for this region.",
"LOW. It is now substantially lower the 3 month average for this region.",
"MODERATE. It is about the same as the 3 month average for this region.",
"HIGH. It is measurably higher than the 3 month average for this region.",
"VERY HIGH. It is now substantially higher than the 3 month average for this region."
)

TLEVEL_COLORS <- c("#BBBBBB", "#FF85FF", "#9289D6", "#3288BD", "#E6F598", "#FDAE61", "#D53E4F")
TLEVEL_STRINGS <- c("UNKNOWN", "SURGING", "DECR. RAPIDLY", "DECREASING", "STABLE", "INCREASING", "INCR. RAPIDLY")
TLEVEL_DESCRIPTIONS <- c(
"Unfortunately, there is not enough recent data to determine whether this disease is trending up or down.",
"It has INCREASED ABRUPTLY within the last week; however, this trend is not consistent.",
"It has been RAPIDLY DECREASING over the last several weeks.",
"It has generally DECREASED over the last several weeks.",
"It has been relatively STABLE over the last several weeks.",
"It has generally INCREASED over the last several weeks.",
"It has been RAPIDLY INCREASING over the last several weeks."
)


MAP_FOOTNOTES = paste0(
	"Click any county on the map above to see specific data for that region. ", 
	"Click anywhere else on the map to return to statewide results. Each map ",
	"county is colored according to the amount (abundance) of the disease present, ", 
	"with the county border colored according to the trend of the disease. ", 
	"IMPORTANT: Data is subject to change as additional sites report results. ",
	"The dashboard updates weekly except during major holidays.",
	sep=""
)

#Brewer Spectral Palette
# [1] "#9E0142" "#D53E4F" "#F46D43" "#FDAE61" "#FEE08B" "#FFFFBF" "#E6F598" "#ABDDA4"
# [9] "#66C2A5" "#3288BD" "#5E4FA2"

abundance_level_colors = c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE")
names(abundance_level_colors) <- c("LOW", "MODERATE", "HIGH", "VERY HIGH", "UNKNOWN")

trend_level_colors = c("#3288BD", "#E6F598", "#FDAE61", "#D53E4F", "#EEEEEE", "#FF85FF", "#96F786")
names(trend_level_colors) <- c("DECREASING", "STABLE", "VARIABLE", "INCREASING", "INDETERMINATE", "SPIKING", "DESPIKING")

lab_colors = c("#EAAA00", "#00B140")
names(lab_colors) <- c("ZooWVU", "MUIDSL")

default_colors = c("#ff0000")
