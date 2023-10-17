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


var_df <- as.data.frame(read.table("analysis/updates/sarvardb.LATEST.txt", header=TRUE, sep="\t", check.names = FALSE))

var_df$day_end <- as_date(mdy_hm(var_df$end_datetime))
var_df$day_start <- as_date(mdy_hm(var_df$start_datetime))
var_df$week <- floor_date(as_date(mdy_hm(var_df$end_datetime)), unit="week")
var_df$percent <- as.numeric(var_df$proportion) * 100

recent_var_df <- subset(var_df, day_end > today() - days(90))

recent_weekly_var_df <- recent_var_df %>% filter(percent > 10) %>% 
						 group_by(week, variant) %>% 
						 summarize(MUT = mean(percent, na.rm = TRUE)) %>% 
						 arrange(week, variant)

p_var <- ggplot(recent_weekly_var_df, aes(fill=variant, y=MUT, x=week)) + 
	geom_bar(position="stack", stat="identity") + 
	theme(legend.position = "bottom")


pcr_df <- as.data.frame(read.table("data/watch_dashboard.LATEST.txt", header=TRUE, sep="\t", check.names = FALSE))

pcr_df$day_end <- as_date(mdy(pcr_df$"Sample Composite End"))
#pcr_df$week <- floor_date(as_date(mdy(pcr_df$day_end)), unit="week")

recent_pcr_df <- subset(pcr_df, day_end > today() - days(90))

#display_pcr_df <- recent_pcr_df %>% filter(location_common_name == "Morgantown Star City" | location_common_name == "Cheat Lake")
display_pcr_df <- recent_pcr_df %>% filter(location_common_name == "Morgantown Star City")

p_pcr <- ggplot(display_pcr_df) + labs(y = "log2(viral particles per person)", x = "") + 
	geom_point(aes(y=log2(n2.loadcap), x=day_end), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = day_end, y = log2(n2.loadcap)), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = day_end, y = log2(n2.loadcap.day5.mean)), color = "#000000", na.rm = TRUE) + 
#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in Wastewater for Morgantown WV (", max(pcr_df$day_end), ")", sep=""), subtitle = "On a log2 scale. Black line is the 5 day rolling trend") + 
	scale_y_continuous(labels = comma) + 
	scale_x_date(breaks = "1 week", labels = format_dates)








