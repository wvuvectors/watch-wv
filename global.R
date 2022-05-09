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
#library(leaflet.extras)
library(rsconnect)

library(sf)
library(tigris)
options(tigris_use_cache = TRUE)

#rsconnect::deployApp('path/to/your/app')


MAP_CENTERS <- data.frame("layer" = c("WWTP", "Sewer Network"),
													"lat" = c(38.951883, 39.642414),
													"lng" = c(-81.0534217, -79.9792327),
													"zoom" = c(8, 13))

INFECTIONS <- c("SARS-CoV-2")
INFECTIONS_DEFAULT <- "SARS-CoV-2"

TARGETS <- c("Mean N1 & N2", "N1", "N2")
TARGET_VALUES <- c("n1n2", "n1", "n2")
TARGETS_DEFAULT <- "n1n2"

TARGETS_DF <- data.frame("infection" = c("SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2"),
												 "target_name" = c("Mean N1 & N2", "N1", "N2"),
												 "target_value" = c("n1n2", "n1", "n2"),
												 "target_color" = c("blue", "dark orange", "dark green")
												)

SMOOTHER_OPTIONS <- c(1, 3, 5, 7, 10, 14)
SMOOTHER_DEFAULT <- 3

L_per_gal <- 3.78541

Sys.setenv(TZ="America/New_York")
today <- Sys.Date()
first_day <- as_date("2021-07-01")
last_day <- today

my_theme <- function () { 
	theme_bw() + theme(axis.text = element_text(size = 8),
										 axis.title = element_text(size = 9, color="#333333"),
										 axis.line.x = element_line(color="#bbbbbb", size=1),
										 axis.line.y = element_line(color="#bbbbbb", size=1),
										 axis.ticks.length.y = unit(-0.5, "cm"), 
										 strip.text = element_text(size = 8),
										 panel.grid.major = element_line(color="#eeeeee", size=1), 
										 #panel.grid.minor.y = element_blank(),
										 panel.grid.minor.x = element_line(color="#eeeeee", size=1),
										 panel.background = element_rect(fill="transparent"), 
										 panel.border = element_rect(fill=NA, color="#bbbbbb", size=1), 
#										 panel.border = element_blank(), 
										 #strip.background = element_rect(fill = 'white', color = 'white'),
										 legend.position = "none",
#										 legend.position = c(.95, .95),
										 legend.justification = c("right", "top"),
										 legend.box.just = "right",
										 legend.margin = margin(6, 6, 6, 6),
										 legend.title = element_text(size = 10, color = "#888888"),
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

#county_sf <- read_sf("data/WV_Counties/County_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
state_sf <- read_sf("data/WV_State/State_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")

watch_file = "data/watch_dashboard.LATEST.txt"
df_watch <- as.data.frame(read.table(watch_file, sep="\t", header=TRUE, check.names=FALSE))
df_watch <- df_watch %>% filter(status == "active" | status == "new")


# Convert date strings into Date objects

df_watch$"Sample Composite Start" <- mdy_hm(df_watch$"Sample Composite Start")
df_watch$"Sample Composite End" <- mdy_hm(df_watch$"Sample Composite End")

df_watch$day <- as_date(df_watch$"Sample Composite End")
df_watch$week_starting <- floor_date(df_watch$day, "week", week_start = 1)
df_watch$week_ending <- ceiling_date(df_watch$day, "week", week_start = 0)

#df_watch$week_num <- week(df_watch$week_ending)
#df_watch$week_alt <- 1 + (df_watch$week_num %% 2)


# Set date constraints on input data
# May want to change this in the future?
df_watch <- df_watch %>% filter(day >= first_day & day <= last_day)


# Make some simpler aliases for common numeric columns
df_watch$daily_flow = df_watch$"Sample Flow (MGD)"

# Clean up some NA entries
df_watch <- df_watch %>% mutate(n1 = replace_na(n1, 0))
df_watch <- df_watch %>% mutate(n2 = replace_na(n2, 0))
df_watch <- df_watch %>% mutate(daily_flow = replace_na(daily_flow, 0))
df_watch <- df_watch %>% mutate(n1n2 = replace_na(n1n2, 0))

df_watch$n1n2.day1.mean <- df_watch$n1n2
df_watch$n1n2.day1.ci <- 0
df_watch$n1n2.load.day1.mean <- df_watch$n1n2.load
df_watch$n1n2.load.day1.ci <- 0

baselines <- c("n1n2.day1.mean", "n1n2.day3.mean", "n1n2.day5.mean", "n1n2.day7.mean", "n1n2.day10.mean", "n1n2.day14.mean", "n1n2.load.day1.mean", "n1n2.load.day3.mean", "n1n2.load.day5.mean", "n1n2.load.day7.mean", "n1n2.load.day10.mean", "n1n2.load.day14.mean")

df_baseline <- data.frame(location_common_name = unique(df_watch$location_common_name))
for (baseline in baselines) {
#	base_mean <- paste0(baseline, ".mean", sep="")
#	base_ci <- paste0(baseline, ".ci", sep="")

	base_day_new <- paste0(baseline, ".day.baseline", sep="")
	base_new <- paste0(baseline, ".baseline", sep="")
#	base_mean_new <- paste0(baseline, ".mean.baseline", sep="")
#	base_ci_new <- paste0(baseline, ".ci.baseline", sep="")

	df_this <- df_watch %>% group_by(location_common_name) %>% 
						 slice_min(.data[[baseline]], n=1, with_ties = FALSE) %>% 
						 select(location_common_name, !!base_day_new := day, !!base_new := baseline)
	df_baseline <- left_join(df_baseline, df_this, by="location_common_name")
}

df_watch <- left_join(df_watch, df_baseline, by="location_common_name")


#df_wwtp <- df_watch %>% filter(daily_flow > 0) %>% mutate(signal_level_1 = 100 * (n1n2.load.day1.mean - n1n2.load.day1.mean.baseline)/n1n2.load.day1.mean.baseline,
#																													signal_level_3 = 100 * (n1n2.load.day3.mean - n1n2.load.day3.mean.baseline)/n1n2.load.day3.mean.baseline,
#																													signal_level_5 = 100 * (n1n2.load.day5.mean - n1n2.load.day5.mean.baseline)/n1n2.load.day5.mean.baseline,
#																													signal_level_7 = 100 * (n1n2.load.day7.mean - n1n2.load.day7.mean.baseline)/n1n2.load.day7.mean.baseline,
#																													signal_level_10 = 100 * (n1n2.load.day10.mean - n1n2.load.day10.mean.baseline)/n1n2.load.day10.mean.baseline,
#																													signal_level_14 = 100 * (n1n2.load.day14.mean - n1n2.load.day14.mean.baseline)/n1n2.load.day14.mean.baseline)
#df_swr <- df_watch %>% filter(daily_flow == 0) %>% mutate(signal_level_1 = 100 * (n1n2.day1.mean - n1n2.day1.mean.baseline)/n1n2.day1.mean.baseline,
#																													signal_level_3 = 100 * (n1n2.day3.mean - n1n2.day3.mean.baseline)/n1n2.day3.mean.baseline,
#																													signal_level_5 = 100 * (n1n2.day5.mean - n1n2.day5.mean.baseline)/n1n2.day5.mean.baseline,
#																													signal_level_7 = 100 * (n1n2.day7.mean - n1n2.day7.mean.baseline)/n1n2.day7.mean.baseline,
#																													signal_level_10 = 100 * (n1n2.day10.mean - n1n2.day10.mean.baseline)/n1n2.day10.mean.baseline,
#																													signal_level_14 = 100 * (n1n2.day14.mean - n1n2.day14.mean.baseline)/n1n2.day14.mean.baseline)

df_wwtp <- df_watch %>% filter(daily_flow > 0) %>% mutate(signal_level_5 = (n1n2.load.day5.mean - n1n2.load.day5.mean.baseline)/n1n2.load.day5.mean.baseline)
df_swr <- df_watch %>% filter(daily_flow == 0) %>% mutate(signal_level_5 = (n1n2.day5.mean - n1n2.day5.mean.baseline)/n1n2.day5.mean.baseline)

df_watch <- rbind(df_wwtp, df_swr)

#alertPal <- colorFactor(c("green", "blue", "#000000", "yellow", "orange", "red"), as.factor(df_watch$alertLevel))

#scale_colour_manual(values = alertPal)


TREND_TXT <- "No trend yet"

alertPal <- colorBin(
	palette = c("green", "black", "yellow", "orange", "red", "dark red"),
	bins = c(-100, 0, 1, 11, 101, 1001, 100001),
	na.color = "#808080",
#	reverse=TRUE,
	domain = df_watch$signal_level_5)
  
