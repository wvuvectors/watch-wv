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


MAP_CENTERS <- data.frame("layer" = c("WWTP", "Sewer Network"),
													"lat" = c(38.991883, 39.642414),
													"lng" = c(-81.6534217, -79.9792327),
													"zoom" = c(8, 13))

INFECTIONS <- c("SARS-CoV-2")
INFECTIONS_DEFAULT <- "SARS-CoV-2"

TARGETS <- c("Mean N1 & N2", "N1", "N2")
TARGET_VALUES <- c("n1n2", "n1", "n2")
TARGETS_DEFAULT <- "n1n2"

TARGET_COLORS <- c("#4b1f6f", "#0084d1", "#aecf00")
TARGET_FILLS <- c("#4b1f6f", "#0084d1", "#aecf00")

TARGETS_DF <- data.frame("infection" = c("SARS-CoV-2", "SARS-CoV-2", "SARS-CoV-2"),
												 "target_name" = c("Mean N1 & N2", "N1", "N2"),
												 "target_value" = c("n1n2", "n1", "n2")
												)

SIGNAL_TREND_WINDOW <- 5
STRENGTH_WEIGHT <- 1
TREND_WEIGHT <- 1.5

SMOOTHER_OPTIONS <- c(1, 3, 5, 7, 10, 14)
SMOOTHER_DEFAULT <- 5

L_per_gal <- 3.78541

Sys.setenv(TZ="America/New_York")
today <- Sys.Date()

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

#county_sf <- read_sf("data/WV_Counties/County_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
state_sf <- read_sf("data/WV_State/State_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")

watch_file = "data/watch_dashboard.LATEST.txt"
df_watch_pre <- as.data.frame(read.table(watch_file, sep="\t", header=TRUE, check.names=FALSE))
df_watch_pre <- df_watch_pre %>% filter(status == "active" | status == "new")


# Convert date strings into Date objects

df_watch_pre$"Sample Composite Start" <- mdy(df_watch_pre$"Sample Composite Start")
df_watch_pre$"Sample Composite End" <- mdy(df_watch_pre$"Sample Composite End")
df_watch_pre$"Sample Received Date" <- mdy(df_watch_pre$"Sample Received Date")

df_watch_pre$day <- as_date(df_watch_pre$"Sample Composite End")
df_watch_pre$week_starting <- floor_date(df_watch_pre$day, "week", week_start = 1)
df_watch_pre$week_ending <- ceiling_date(df_watch_pre$day, "week", week_start = 7)

df_watch_pre$day_received <- as_date(df_watch_pre$"Sample Received Date")
#df_watch_pre$week_num <- week(df_watch_pre$week_ending)
#df_watch_pre$week_alt <- 1 + (df_watch_pre$week_num %% 2)


# Set date constraints on input data
# May want to change this in the future?
last_day <- max(df_watch_pre$day)
first_day <- last_day - years(1)
df_watch_pre <- df_watch_pre %>% filter(day >= first_day & day <= last_day)

# TEMP KLUDGE!!!
#df_watch_pre <- df_watch_pre %>% filter(location_common_name != "Huntington")

# Set the primary plot column names for each group
#df_watch_pre <- df_watch_pre %>% mutate(plot_view_col = ifelse(group == "WWTP", yes = "n1n2.loadcap.day5.mean", no = "n1n2.day5.mean"),
#																				plot_view_colbase = ifelse(group == "WWTP", yes = "n1n2.loadcap.day5.mean.baseline", no = "n1n2.day5.mean.baseline"))

# Make some simpler aliases for common numeric columns
#df_watch_pre$daily_flow = df_watch_pre$"Sample Flow (MGD)"

df_watch_pre$n1n2.day1.mean <- df_watch_pre$n1n2
df_watch_pre$n1n2.day1.ci <- 0

df_watch_pre$n1n2.load.day1.mean <- df_watch_pre$n1n2.load
df_watch_pre$n1n2.load.day1.ci <- 0

df_watch_pre$n1n2.loadcap.day1.mean <- df_watch_pre$n1n2.loadcap
df_watch_pre$n1n2.loadcap.day1.ci <- 0

baselines <- c("n1n2", "n1n2.load", "n1n2.loadcap", "n1n2.day5.mean", "n1n2.load.day5.mean", "n1n2.loadcap.day5.mean")

df_baseline <- data.frame(location_common_name = unique(df_watch_pre$location_common_name))
for (baseline in baselines) {
#	base_mean <- paste0(baseline, ".mean", sep="")
#	base_ci <- paste0(baseline, ".ci", sep="")

	base_day_new <- paste0(baseline, ".day.baseline", sep="")
	base_new <- paste0(baseline, ".baseline", sep="")
#	base_mean_new <- paste0(baseline, ".mean.baseline", sep="")
#	base_ci_new <- paste0(baseline, ".ci.baseline", sep="")

	df_this <- df_watch_pre %>% group_by(location_common_name) %>% 
						 slice_min(.data[[baseline]], n=1, with_ties = FALSE) %>% 
						 select(location_common_name, !!base_day_new := day, !!base_new := baseline)
	df_baseline <- left_join(df_baseline, df_this, by="location_common_name")
}

df_watch_pre <- left_join(df_watch_pre, df_baseline, by="location_common_name")

# calculate signal strength, fold change, and percent change from baseline for each day
df_wwtp <- df_watch_pre %>% filter(group == "WWTP" & !is.na(daily_flow) & !is.na(n1) & !is.na(n2) & !is.na(n1n2.day5.mean))
df_wwtp <- df_wwtp %>% mutate(fold_change_smoothed = (n1n2.loadcap.day5.mean - n1n2.loadcap.day5.mean.baseline)/n1n2.loadcap.day5.mean.baseline,
															percent_change_smoothed = 100*(n1n2.loadcap.day5.mean - n1n2.loadcap.day5.mean.baseline)/n1n2.loadcap.day5.mean.baseline,
															signal_strength_smoothed = n1n2.loadcap.day5.mean,
															fold_change = (n1n2.loadcap - n1n2.loadcap.day5.mean.baseline)/n1n2.loadcap.day5.mean.baseline,
															signal_strength = n1n2.loadcap)
df_wwtp <- df_wwtp[!is.na(df_wwtp$n1n2.loadcap.day5.mean), ]												

df_swr <- df_watch_pre %>% filter(group == "Sewer Network" & !is.na(n1) & !is.na(n2) & !is.na(n1n2.day5.mean))
df_swr <- df_swr %>% mutate(fold_change_smoothed = (n1n2.day5.mean - n1n2.day5.mean.baseline)/n1n2.day5.mean.baseline, 
														percent_change_smoothed = 100*(n1n2.day5.mean - n1n2.day5.mean.baseline)/n1n2.day5.mean.baseline,
														signal_strength_smoothed = n1n2.day5.mean,
														fold_change = (n1n2 - n1n2.day5.mean.baseline)/n1n2.day5.mean.baseline,
														signal_strength = n1n2)
df_swr <- df_swr[!is.na(df_swr$n1n2.day5.mean), ]												

df_watch_pre <- rbind(df_wwtp, df_swr)
df_watch <- df_watch_pre


# Calculate signal trend (trajectory) for most recent window
# Uses the raw signal strength for the n most recent samples (copies/L for sewers, or per capita load for treatment plants)
df_trend <- df_watch %>% 
						 group_by(location_common_name, day) %>% 
						 arrange(day) %>%
						 summarize(mean_signal_strength = mean(signal_strength, na.rm = TRUE)) %>% 
						 slice_tail(n=SIGNAL_TREND_WINDOW) %>%
						 summarize(current_scaled_signal_trend = coef(lm(formula = ((mean_signal_strength-min(mean_signal_strength, na.rm = TRUE))/(max(mean_signal_strength, na.rm = TRUE)-min(mean_signal_strength, na.rm = TRUE))) ~ day))[2],
						 					 current_signal_trend = coef(lm(formula = mean_signal_strength ~ day))[2])


# Calculate strength and fold change for most recent day
# Uses the 5-day mean already pre-calculated
df_level <- df_watch %>% group_by(location_common_name) %>% summarize(current_fold_change_smoothed = fold_change_smoothed[day == max(day)])
df_strength <- df_watch %>% group_by(location_common_name) %>% summarize(current_signal_strength_smoothed = signal_strength_smoothed[day == max(day)])


# Merge signal trend, strength, and fold change data frames
df_currents <- left_join(df_level, df_strength, by=c("location_common_name" = "location_common_name"))
df_signal <- left_join(df_trend, df_currents, by=c("location_common_name" = "location_common_name"))

# Add group identifiers
df_loc <- df_watch %>% select(location_common_name, group)
df_signal <- left_join(df_signal, distinct(df_loc), by=c("location_common_name" = "location_common_name"))

df_maxmins <- df_watch %>% 
							group_by(location_common_name) %>% 
							summarize(max_signal_strength = max(signal_strength, na.rm = TRUE),
												min_signal_strength = min(signal_strength, na.rm = TRUE),
												max_signal_strength_smoothed = max(signal_strength_smoothed, na.rm = TRUE),
												min_signal_strength_smoothed = min(signal_strength_smoothed, na.rm = TRUE),
												max_fold_change_smoothed = max(fold_change_smoothed, na.rm = TRUE),
												min_fold_change_smoothed = min(fold_change_smoothed, na.rm = TRUE),
												max_percent_change_smoothed = max(percent_change_smoothed, na.rm = TRUE),
												min_percent_change_smoothed = min(percent_change_smoothed, na.rm = TRUE))
df_signal <- left_join(df_signal, df_maxmins, by=c("location_common_name" = "location_common_name"))


# Scale strength so it falls between 0-1
# Scale trend so it falls between -1 and 1
# Weight accordingly
# Take the sum of the weighted & scaled values
# scaled signal goes from -1 (very low and decreasing rapidly) to 2 (very high and increasing rapidly)
# weighted signal goes from total weight times scaled signal range
#
#df_signal <- df_signal %>% mutate(scaled_signal_strength = (current_signal_strength - min_signal_strength)/(max_signal_strength - min_signal_strength),
#																	scaled_fold_change = (current_fold_change - min_fold_change)/(max_fold_change - min_fold_change),
#																	scaled_signal_trend = current_signal_trend)
#																	#scaled_signal_trend = ifelse(current_signal_trend < 0, yes = ((current_signal_trend - min_signal_trend)/(0 - min_signal_trend))-1, no = current_signal_trend/max_signal_trend))
#
#df_signal <- df_signal %>% mutate(scaled_signal_strength = ifelse(scaled_signal_strength > 1, yes = 1, no = scaled_signal_strength),
#																	scaled_fold_change = ifelse(scaled_fold_change > 1, yes = 1, no = scaled_fold_change))
#																	#scaled_signal_trend = ifelse(scaled_signal_trend > 1, yes = 1, no = scaled_signal_trend))
#
#df_signal <- df_signal %>% mutate(scaled_signal = (STRENGTH_WEIGHT * scaled_signal_strength) + (TREND_WEIGHT * scaled_signal_trend))
#
wt <- STRENGTH_WEIGHT * TREND_WEIGHT
SIGNAL_BINS <- data.frame(fold_change_smoothed = c(0, 26, 101, 251, 501, 5001),
													signal_trend = c(-500000001, -50, -2, 2, 50, 500000001)) 
													#scaled_signal_indicator = c(wt*-1.1, wt*-0.49, wt*0.16, wt*0.34, wt*0.67, wt*1.01, wt*1.51, wt*2.1))

SIGNAL_CODES <- data.frame(change = c("low", "moderate", "high", "very high", "extremely high"), 
													 trend = c("downward", "downward", "stable", "upward", "upward"),
													 color = c("#579d1c", "#ffd320", "#eedc82", "#ff950e", "#c5000b"))

names(TARGET_COLORS) <- levels(factor(c(levels(as.factor(TARGET_VALUES)))))
names(TARGET_FILLS) <- levels(factor(c(levels(as.factor(TARGET_VALUES)))))

ALERT_COLORS = c("#579d1c", "#ffd320", "#eedc82", "#ff950e", "#c5000b")
ALERT_TXT = c("low", "moderate", "high", "very high", "extremely high")

alertPal <- colorBin(
	palette = SIGNAL_CODES$color,
	bins = SIGNAL_BINS$fold_change_smoothed,
	na.color = "#808080",
#	reverse=TRUE,
	domain = df_watch$fold_change_smoothed)
  

targetPal <- colorFactor(
	palette = TARGETS_DF$target_color,
	domain = TARGETS_DF$target_value,
	ordered = TRUE,
	na.color = "#aaaaaa",
	alpha = TRUE
)
