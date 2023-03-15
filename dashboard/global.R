source("version.R")

# libraries
library(tidyverse)
library(dplyr)
library(data.table)
library(zoo)

library(ggplot2)
library(ggthemes)
library(viridis)
library(RColorBrewer)

library(scales)
library(lubridate)

library(plotly)
library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(leaflet)
library(fontawesome)
library(DT)
#library(leaflet.extras)
library(rsconnect)

library(sf)
library(tigris)
options(tigris_use_cache = TRUE)

#rsconnect::deployApp('path/to/your/app')

# color palette is ggthemes$calc
# 1 Chart 1  #004586	dark blue
# 2 Chart 2  #ff420e	red orange
# 3 Chart 3  #ffd320	yellow
# 4 Chart 4  #579d1c	green
# 5 Chart 5  #7e0021	brown
# 6 Chart 6  #83caff	light blue
# 7 Chart 7  #314004	dark green
# 8 Chart 8  #aecf00	light green
# 9 Chart 9  #4b1f6f	purple
#10 Chart 10 #ff950e	orange
#11 Chart 11 #c5000b	dark red
#12 Chart 12 #0084d1	Carolina blue


GEOLEVELS <- c("Facility", "County", "Region")
GEOLEVELS_DEFAULT <- "County"

TARGETS <- c("SARS-CoV-2", "Influenza A", "Influenza B", "RSV")
TARGETS_DEFAULT <- "SARS-CoV-2"

LOCI <- c("N1", "N2")
LOCI_DEFAULT <- "N1"


MAP_CENTERS <- data.frame("layer" = c("WWTP", "Sewer Network"),
													"lat" = c(38.951883, 39.642414),
													"lng" = c(-80.0534217, -79.9792327),
													"zoom" = c(7, 13))


Sys.setenv(TZ="America/New_York")
today <- Sys.Date()

plot_theme <- function () { 
	theme(axis.text = element_text(size = 8),
				axis.title = element_text(size = 9, color="#333333"),
				axis.line.x = element_line(color="#bbbbbb", size=1),
				axis.line.y = element_line(color="#bbbbbb", size=1),
				axis.ticks.length.y = unit(-0.5, "cm"), 
				strip.text = element_text(size = 8),
				panel.grid.major = element_line(color="#eeeeee", size=1), 
				panel.grid.minor.x = element_line(color="#eeeeee", size=1),
				panel.background = element_rect(fill="transparent"), 
				panel.border = element_rect(fill=NA, color="#bbbbbb", size=1), 
				legend.position = "none",
				#legend.justification = c("right", "top"),
				#legend.box.just = "right",
				#legend.margin = margin(6, 6, 6, 6),
				legend.title = element_text(size = 10, color = "#888888"),
				legend.background = element_rect(fill="transparent"), 
				legend.text = element_text(size = 8, color = "#333333"),
				plot.background = element_rect(fill="transparent"), 
				plot.title = element_text(size = 10, color="#045a8d", hjust=0.5)
)}

ci90 <- function(x) {
	0.5 * qt(0.80, length(x) - 1) * (sd(x) / sqrt(length(x)))
#	m <- mean(x)
#	se <- sd(x) / sqrt(x.length)
#	ci = qt(1 - (0.05 / 2), x.length - 1) * se
#	lower.ci = m - ci
#	upper.ci = m + ci
}

ci95 <- function(x) {
	0.5 * qt(0.95, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

ci98 <- function(x) {
	0.5 * qt(0.98, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

ci99 <- function(x) {
	0.5 * qt(0.99, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

se <- function(x) {
	0.5 * sd(x) / sqrt(length(x))
}

format_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
		true = paste(day, month, years, sep = "\n"), 
		false = if_else(is.na(lag(month)) | lag(month) != month, true = paste(day, month, sep = "\n"), false = day)
	)
}



# Load the data files

county_sf <- read_sf("shapefiles/wv_counties/County_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
state_sf <- read_sf("shapefiles/wv_state/State_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")

df_county <- as.data.frame(read.table("data/county.txt", sep="\t", header=TRUE, check.names=FALSE))
df_lab <- as.data.frame(read.table("data/lab.txt", sep="\t", header=TRUE, check.names=FALSE))
df_location <- as.data.frame(read.table("data/location.txt", sep="\t", header=TRUE, check.names=FALSE))
df_result <- as.data.frame(read.table("data/result.txt", sep="\t", header=TRUE, check.names=FALSE))
#df_sample <- as.data.frame(read.table("data/sample.txt", sep="\t", header=TRUE, check.names=FALSE))
df_wwtp <- as.data.frame(read.table("data/wwtp.txt", sep="\t", header=TRUE, check.names=FALSE))


# Convert date strings into Date objects
df_result$collection_start_datetime <- mdy_hm(df_result$collection_start_datetime)
df_result$collection_end_datetime <- mdy_hm(df_result$collection_end_datetime)
#df_sample$sample_collection_start_datetime <- mdyhm(df_sample$sample_collection_start_datetime)
#df_sample$sample_collection_end_datetime <- mdyhm(df_sample$sample_collection_end_datetime)
#df_sample$sample_recovered_datetime <- mdyhm(df_sample$sample_recovered_datetime)
#df_sample$sample_received_date <- mdy(df_sample$sample_received_date)

# Create a plot date
df_result$date_to_plot <- as_date(df_result$collection_end_datetime)

#df_sample$week_starting <- floor_date(df_sample$day, "week", week_start = 1)
#df_sample$week_ending <- ceiling_date(df_sample$day, "week", week_start = 0)

#week_starting <- floor_date(df_result$collection_end_datetime, "week", week_start = 1)

