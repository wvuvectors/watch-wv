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
													"lat" = c(38.951883, 39.642414),
													"lng" = c(-81.0534217, -79.9792327),
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
first_day <- as_date("2021-07-01")
last_day <- today

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

#county_sf <- read_sf("data/WV_Counties/County_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")
state_sf <- read_sf("data/WV_State/State_wld.shp") %>% st_transform("+proj=longlat +datum=WGS84 +no_defs")

watch_file = "data/watch_dashboard.LATEST.txt"
df_watch_pre <- as.data.frame(read.table(watch_file, sep="\t", header=TRUE, check.names=FALSE))
df_watch_pre <- df_watch_pre %>% filter(status == "active" | status == "new")


# Convert date strings into Date objects

df_watch_pre$"Sample Composite Start" <- mdy_hm(df_watch_pre$"Sample Composite Start")
df_watch_pre$"Sample Composite End" <- mdy_hm(df_watch_pre$"Sample Composite End")
df_watch_pre$"Sample Received Date" <- mdy(df_watch_pre$"Sample Received Date")

df_watch_pre$day <- as_date(df_watch_pre$"Sample Composite End")
df_watch_pre$week_starting <- floor_date(df_watch_pre$day, "week", week_start = 1)
df_watch_pre$week_ending <- ceiling_date(df_watch_pre$day, "week", week_start = 0)

df_watch_pre$day_received <- as_date(df_watch_pre$"Sample Received Date")

#df_watch_pre$week_num <- week(df_watch_pre$week_ending)
#df_watch_pre$week_alt <- 1 + (df_watch_pre$week_num %% 2)


# Set date constraints on input data
# May want to change this in the future?
df_watch_pre <- df_watch_pre %>% filter(day >= first_day & day <= last_day)


# Make some simpler aliases for common numeric columns
df_watch_pre$daily_flow = df_watch_pre$"Sample Flow (MGD)"

# Clean up some NA entries
df_watch_pre <- df_watch_pre %>% mutate(n1 = replace_na(n1, 0))
df_watch_pre <- df_watch_pre %>% mutate(n2 = replace_na(n2, 0))
df_watch_pre <- df_watch_pre %>% mutate(daily_flow = replace_na(daily_flow, 0))
df_watch_pre <- df_watch_pre %>% mutate(n1n2 = replace_na(n1n2, 0))

df_watch_pre$n1n2.day1.mean <- df_watch_pre$n1n2
df_watch_pre$n1n2.day1.ci <- 0
df_watch_pre$n1n2.load.day1.mean <- df_watch_pre$n1n2.load
df_watch_pre$n1n2.load.day1.ci <- 0

baselines <- c("n1n2.day5.mean", "n1n2.load.day5.mean")

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

# calculate signal level for each day
df_wwtp <- df_watch_pre %>% filter(daily_flow > 0) %>% mutate(signal_level = (n1n2.load.day5.mean - n1n2.load.day5.mean.baseline)/n1n2.load.day5.mean.baseline)
df_swr <- df_watch_pre %>% filter(daily_flow == 0) %>% mutate(signal_level = (n1n2.day5.mean - n1n2.day5.mean.baseline)/n1n2.day5.mean.baseline)
df_watch_pre <- rbind(df_wwtp, df_swr)

df_watch <- df_watch_pre
df_watch <- df_watch[!is.na(df_watch$n1n2.load.day5.mean) | !is.na(df_watch$n1n2.day5.mean), ]												


# calculate signal trend (trajectory) for most recent day
df_signal <- df_watch %>% 
						 group_by(location_common_name, day) %>% 
						 arrange(day) %>%
						 summarize(delta = mean(signal_level, na.rm = TRUE)) %>% 
						 slice_tail(n=SIGNAL_TREND_WINDOW) %>%
						 summarize(trend = coef(lm(formula = delta ~ day))[2])

# merge signal trend and strength
tmp <- df_watch %>% group_by(location_common_name) %>% summarize(strength = signal_level[day == max(day)])
#max_strength = max(signal_level, na.rm = TRUE),
#max_strength = ifelse(max(signal_level, na.rm = TRUE) > 1000, yes = 1000, no = max(signal_level, na.rm = TRUE)),
#min_strength = min(signal_level, na.rm = TRUE))
df_signal <- left_join(df_signal, tmp, by=c("location_common_name" = "location_common_name"))


# Scale strength so it falls between 0-1
# Scale trend so it falls between -1 and 1
# Weight accordingly
# Take the sum of the weighted & scaled values
# scaled signal goes from -1 (very low and decreasing rapidly) to 2 (very high and increasing rapidly)
#
df_signal <- df_signal %>% mutate(max_strength = max(strength, na.rm = TRUE),
																	min_trend = ifelse(trend < -30, yes = -30, no = min(trend, na.rm = TRUE)),
																	max_trend = ifelse(trend > 30, yes = 30, no = max(trend, na.rm = TRUE)),
																	scaled_trend = ifelse(trend < 0, yes = ((trend - min_trend)/(0 - min_trend))-1, no = trend/max_trend))

df_signal <- df_signal %>% mutate(scaled_strength = ifelse(max_strength > 1000, yes = strength/1000, no = strength/max_strength),
																	scaled_signal = (STRENGTH_WEIGHT * scaled_strength) + (TREND_WEIGHT * scaled_trend))

wt <- STRENGTH_WEIGHT * TREND_WEIGHT
SIGNAL_BINS <- data.frame(scaled_strength = c(-1.1, 0, .03, .1, .26, .51, .76, 1.1), 
													scaled_trend = c(-1.1, -0.75, -0.5, -0.25, 0.11, 0.3, 0.67, 1.01), 
													scaled_indicator = c(wt*-1.1, wt*-0.49, wt*0.16, wt*0.34, wt*0.67, wt*1.01, wt*1.51, wt*2.1))

SIGNAL_CODES <- data.frame(strength = c("undetectable", "very low", "low", "moderate", "high", "very high", "extremely high"), 
													 trend = c("decreasing rapidly", "decreasing moderately", "decreasing slightly", "relatively stable", "increasing", "increasing substantially", "increasing rapidly"),
													 color = c("#579d1c", "#579d1c", "#ffd320", "#ff950e", "#ff950e", "#c5000b", "#EE2221"))

names(TARGET_COLORS) <- levels(factor(c(levels(as.factor(TARGET_VALUES)))))
names(TARGET_FILLS) <- levels(factor(c(levels(as.factor(TARGET_VALUES)))))

ALERT_COLORS = c("#579d1c", "#aecf00", "#ffd320", "#ff950e", "#ff950e", "#ff420e", "#c5000b")
ALERT_TXT = c("Low concern", "Low concern", "Watchful", "Concerning", "Concerning", "Alarming", "Critical")

alertPal <- colorBin(
	palette = SIGNAL_CODES$color,
	bins = SIGNAL_BINS$scaled_indicator,
	na.color = "#808080",
#	reverse=TRUE,
	domain = df_signal$scaled_signal)
  

targetPal <- colorFactor(
	palette = TARGETS_DF$target_color,
	domain = TARGETS_DF$target_value,
	ordered = TRUE,
	na.color = "#aaaaaa",
	alpha = TRUE
)
