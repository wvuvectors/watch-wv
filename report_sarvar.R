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

recent_var_df <- subset(var_df %>% filter(str_detect(facility, "WWTP")), day_end > today() - days(90))

recent_daily_var_df <- recent_var_df %>% filter(percent > 5) %>% 
						 group_by(day_end, lineage) %>% 
						 summarize(MUT = sum(percent, na.rm = TRUE)) %>% 
						 arrange(day_end, lineage)

recent_weekly_var_df <- recent_var_df %>% filter(percent > 5) %>% 
						 group_by(week, lineage) %>% 
						 summarize(MUT = mean(percent, na.rm = TRUE)) %>% 
						 arrange(week, lineage)

recent_daily_starcity_df <- recent_var_df %>% filter(percent > 2 & facility == "StarCityWWTP-01") %>% 
	group_by(day_end, lineage) %>% 
	summarize(MUT = sum(percent, na.rm = TRUE)) %>% 
	arrange(day_end, lineage)


p_daily_var <- ggplot(recent_daily_var_df, aes(fill=lineage, y=MUT, x=day_end)) + 
	geom_bar(position="stack", stat="identity") + 
	scale_x_date(date_breaks = "5 days") + 
	theme(legend.position = "bottom")

p_daily_starcity <- ggplot(recent_daily_starcity_df, aes(fill=lineage, y=MUT, x=day_end)) + 
	geom_bar(position="stack", stat="identity") + 
	scale_x_date(date_breaks = "5 days") + 
	theme(legend.position = "bottom")
	
p_var_stadium <- ggplot(var_df %>% filter(str_detect(facility, "Stadium")), aes(fill=lineage, y=percent, x=facility)) + 
	geom_bar(position="stack", stat="identity") + 
	facet_grid(~day_start) + 
	theme(legend.position = "bottom")







