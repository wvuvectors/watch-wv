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
library(rgdal)

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


GEOLEVELS <- c("Facility", "County")
GEOLEVELS_DEFAULT <- "County"

TARGETS <- c("SARS-CoV-2", "Flu A", "Flu B", "RSV", "Norovirus")
TARGET_PRIMARY <- "SARS-CoV-2"

LOCI <- c("N1", "N2")
LOCUS_PRIMARY <- "N2"

VIEW_RANGES <- c(1, 3, 6, 12, 24)
VIEW_RANGE_PRIMARY <- 12

MAP_CENTER <- list2env(list(lat = 38.951883, lng = -80.0534217, zoom = 7))

Sys.setenv(TZ="America/New_York")
today <- Sys.Date()
#today <- as.Date("2022-07-12")

plot_theme <- function () { 
	theme(axis.text = element_text(size = 8),
				axis.title = element_text(size = 9, color="#333333"),
				axis.line.x = element_line(color="#bbbbbb", linewidth=1),
				axis.line.y = element_line(color="#bbbbbb", linewidth=1),
				axis.ticks.length.y = unit(-0.5, "cm"), 
				strip.text = element_text(size = 8),
				panel.grid.major = element_line(color="#eeeeee", linewidth=1), 
				panel.grid.minor.x = element_line(color="#eeeeee", linewidth=1),
				panel.background = element_rect(fill="transparent"), 
				panel.border = element_rect(fill=NA, color="#bbbbbb", linewidth=1), 
				legend.position = "none",
				legend.justification = c("left", "top"),
				#legend.direction = "horizontal",
				legend.box.just = "center",
				#legend.margin = margin(6, 6, 6, 6),
				legend.title = element_blank(),
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




# Load shape files
county_spdf <- readOGR( 
  dsn= paste0("shapefiles/wv_counties/") , 
  layer="WV_Counties",
  verbose=FALSE
)
#state_spdf <- read_sf("shapefiles/wv_state/State_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
mypalette <- colorBin(palette="YlOrBr", domain=county_spdf@data$POPCH_PCT, na.color="transparent")



# Load data files
df_county <- as.data.frame(read.table("data/county.txt", sep="\t", header=TRUE, check.names=FALSE))
df_lab <- as.data.frame(read.table("data/lab.txt", sep="\t", header=TRUE, check.names=FALSE))
df_location <- as.data.frame(read.table("data/location.txt", sep="\t", header=TRUE, check.names=FALSE))
df_target <- as.data.frame(read.table("data/target.txt", sep="\t", header=TRUE, check.names=FALSE))
df_result <- as.data.frame(read.table("data/result.txt", sep="\t", header=TRUE, check.names=FALSE))
df_sample <- as.data.frame(read.table("data/sample.txt", sep="\t", header=TRUE, check.names=FALSE))
df_wwtp <- as.data.frame(read.table("data/wwtp.txt", sep="\t", header=TRUE, check.names=FALSE))

TARGET_CLASS = unique(df_target %>% filter(target_id == TARGET_PRIMARY))$target_class
DISEASE_PRIMARY = unique(df_target %>% filter(target_id == TARGET_PRIMARY))$target_disease


df_hosp1 <- as.data.frame(read.table("data/hospitalizations.csv", sep=",", header=TRUE, check.names=TRUE))
colnames(df_hosp1)[1] <- "i"
df_hosp1$time_value <- ymd(df_hosp1$time_value)
df_hosp2 <- df_hosp1 %>% group_by(time_value, location_name) %>% summarize(rolling_weekly_mean = sum(value)*7)
df_hosp2$mmr_year <- lubridate::year(df_hosp2$time_value)
df_hosp2$mmr_week <- lubridate::week(df_hosp2$time_value)
df_hospital <- df_hosp2 %>% group_by(mmr_year, mmr_week, location_name) %>% summarize(weekly_sum = round(sum(rolling_weekly_mean), digits = 0))


# Filter out locations that are not active
df_location <- df_location %>% filter(tolower(location_status) == "active")

# Filter out data entries that did not pass QC
df_result <- df_result %>% filter(target_result_validated == 1)
df_sample <- df_sample %>% filter(tolower(sample_qc) == "pass")


# Convert date strings into Date objects
df_result$collection_start_datetime <- mdy_hm(df_result$collection_start_datetime)
df_result$collection_end_datetime <- mdy_hm(df_result$collection_end_datetime)
df_sample$sample_collection_start_datetime <- mdy_hm(df_sample$sample_collection_start_datetime)
df_sample$sample_collection_end_datetime <- mdy_hm(df_sample$sample_collection_end_datetime)
df_sample$sample_recovered_datetime <- mdy_hm(df_sample$sample_recovered_datetime)
df_sample$sample_received_date <- mdy(df_sample$sample_received_date)


# Create a date to plot
df_result$date_to_plot <- as.Date(df_result$collection_end_datetime)

df_result$mmr_year <- lubridate::year(df_result$date_to_plot)
df_result$mmr_week <- lubridate::week(df_result$date_to_plot)

df_result$week_starting <- floor_date(df_result$date_to_plot, "week", week_start = 1)
df_result$week_ending <- ceiling_date(df_result$date_to_plot, "week", week_start = 1)

