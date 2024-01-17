library(tidyverse)
library(dplyr)
library(lubridate)

library(ggplot2)
library(ggthemes)

library(scales)


format_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
		true = paste(day, month, years, sep = "\n"), 
		false = if_else(is.na(lag(month)) | lag(month) != month, true = paste(day, month, sep = "\n"), false = day)
	)
}

MAX_LOADCAP <- 1000

pcr_df <- as.data.frame(read.table("data/watch_dashboard.LATEST.txt", header=TRUE, sep="\t", check.names = FALSE))

pcr_df$day_end <- as_date(mdy(pcr_df$"Sample Composite End"))
#pcr_df$week <- floor_date(as_date(mdy(pcr_df$day_end)), unit="week")

recent_pcr_df <- subset(pcr_df, day_end > today() - days(90))

#display_pcr_df <- recent_pcr_df %>% filter(location_common_name == "Morgantown Star City" | location_common_name == "Cheat Lake")
display_pcr_df <- recent_pcr_df %>% filter(location_common_name == "Morgantown Star City")

p_pcr_log <- ggplot(display_pcr_df) + labs(y = "log2(viral particles per person)", x = "") + 
	geom_point(aes(y=log2(n2.loadcap), x=day_end), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = day_end, y = log2(n2.loadcap)), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = day_end, y = log2(n2.loadcap.day5.mean)), color = "#000000", na.rm = TRUE) + 
#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in Wastewater for Morgantown WV (", max(pcr_df$day_end), ")", sep=""), subtitle = "On a log2 scale. Black line is the 5 day rolling trend.") + 
	scale_y_continuous(labels = comma) + 
	scale_x_date(breaks = "1 week", labels = format_dates)


limited_pcr_df <- display_pcr_df %>% mutate(n2.loadcap = case_when(
	n2.loadcap > MAX_LOADCAP ~ MAX_LOADCAP,
	TRUE ~ n2.loadcap
))

limited_pcr_df <- limited_pcr_df %>% mutate(n2.loadcap.day5.mean = case_when(
	n2.loadcap.day5.mean > MAX_LOADCAP ~ MAX_LOADCAP,
	TRUE ~ n2.loadcap.day5.mean
))

p_pcr_lin <- ggplot(limited_pcr_df) + labs(y = "viral particles per person", x = "") + 
	geom_point(aes(y=n2.loadcap, x=day_end), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = day_end, y = n2.loadcap), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = day_end, y = n2.loadcap.day5.mean), color = "#000000", na.rm = TRUE) + 
#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in Wastewater for Morgantown WV (", max(pcr_df$day_end), ")", sep=""), subtitle = "Black line is the 5 day rolling trend.") + 
	scale_y_continuous(labels = comma, limits = c(0, MAX_LOADCAP)) + 
	scale_x_date(breaks = "1 week", labels = format_dates)
