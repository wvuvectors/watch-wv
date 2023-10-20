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


all_df <- as.data.frame(read.table("analysis/updates/sarvardb.LATEST.txt", header=TRUE, sep="\t", check.names = FALSE))
all_df$day_end <- as_date(mdy_hm(all_df$end_datetime))
all_df$day_start <- as_date(mdy_hm(all_df$start_datetime))
all_df$week <- floor_date(as_date(mdy_hm(all_df$end_datetime)), unit="week")
all_df$percent <- as.numeric(all_df$proportion) * 100

recent_df <- subset(all_df %>% filter(str_detect(facility, "WWTP")), day_end > today() - days(90))

recent_2pct_df <- recent_df %>% filter(percent >= 2)

#recent_weekly_var_df <- recent_var_df %>% filter(percent > 5) %>% 
#						 group_by(week, parental) %>% 
#						 summarize(MUT = mean(percent, na.rm = TRUE)) %>% 
#						 arrange(week, parental)


p_all <- ggplot(recent_2pct_df, aes(fill=parental, y=percent, x=day_end)) + 
	geom_bar(position="stack", stat="identity") + 
	facet_wrap(~facility, nrow=3) + 
	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
	theme(legend.position = "bottom")

p_starcity <- ggplot(recent_2pct_df %>% filter(facility == "StarCityWWTP-01"), aes(fill=parental, y=percent, x=day_end)) + 
	geom_bar(position="stack", stat="identity") + 
	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
	theme(legend.position = "bottom")
	
p_var_stadium <- ggplot(all_df %>% filter(str_detect(facility, "Stadium")), aes(fill=lineage, y=percent, x=facility)) + 
	geom_bar(position="stack", stat="identity") + 
	facet_grid(~day_start) + 
	scale_x_date(date_breaks = "1 day", labels = format_dates) + 
	theme(legend.position = "bottom")







