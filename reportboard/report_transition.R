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

format_dates_wide <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	#day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
					true = paste(month, years, sep = "\n"), 
					false = if_else(is.na(lag(month)) | lag(month) != month, true = paste(month, sep = "\n"), false = month)
	)
}

MAX_LOADCAP <- 1000

#pcr_df <- as.data.frame(read.table("data/merged_nwss.LATEST.csv", header=TRUE, sep=",", check.names = FALSE))
#pcr_df$report_day <- as_date(ymd(pcr_df$sample_collect_date))
#pcr_df$report_week <- floor_date(as_date(ymd(pcr_df$report_day)), unit="week")

# mean N2 over report_week, all facilities:
#sum_n2_df <- pcr_df %>% 
#	filter(pcr_gene_target == "n2") %>% 
#	group_by(report_week) %>% 
#	summarise_at(vars(pcr_target_avg_conc), list(name = mean))


pcr_df <- as.data.frame(read.table("data/watchdb.result.LATEST.txt", header=TRUE, sep="\t", check.names = FALSE))
#pcr_df$report_day <- as_date(parse_date_time(pcr_df$collection_end_date, "%m/%d/%y"))
pcr_df$report_day <- as_date(mdy(pcr_df$collection_end_date))
pcr_df$report_week <- floor_date(as_date(ymd(pcr_df$report_day)), unit="week")

# mean N2 over report_week, all facilities:
sum_n2_df <- pcr_df %>% 
	filter(target_genetic_locus == "N2:SARS" & !is.na(target_copies_fn_per_cap)) %>% 
	group_by(report_week) %>% 
	summarise_at(vars(target_copies_fn_per_cap), list(mean_n2 = mean))


# look at just the last 6 months:
recent_pcr_df <- subset(sum_n2_df, report_week > today() - months(6))

p_recent <- ggplot(recent_pcr_df) + labs(y = "SARS-CoV-2 load per 100 persons", x = "") + 
	geom_point(aes(y=mean_n2, x=report_week), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = report_week, y = mean_n2), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = report_week, y = mean_n2), color = "#000000", na.rm = TRUE) + 
	#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
	#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in Wastewater for WV, as of ", max(sum_n2_df$report_week), sep=""), subtitle = "(weekly mean over all facilities)") + 
	scale_y_continuous() + 
	scale_x_date(breaks = "1 week", labels = format_dates)


p_all <- ggplot(sum_n2_df) + labs(y = "SARS-CoV-2 load per 100 persons", x = "") + 
	geom_point(aes(y=mean_n2, x=report_week), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = report_week, y = mean_n2), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = report_week, y = mean_n2), color = "#000000", na.rm = TRUE) + 
	#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
	#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in Wastewater for WV, as of ", max(sum_n2_df$report_week), sep=""), subtitle = "(weekly mean over all facilities)") + 
	scale_y_continuous() + 
	scale_x_date(breaks = "1 month", labels = format_dates_wide)


