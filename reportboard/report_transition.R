library(tidyverse)
library(dplyr)
library(lubridate)

library(ggplot2)
library(ggthemes)

library(scales)

plot_theme <- function () { 
	theme(axis.text = element_text(size = 8),
				axis.title = element_text(size = 9, color="#333333"),
				axis.line.x = element_line(color="#bbbbbb", linewidth=1),
				axis.line.y = element_line(color="#bbbbbb", linewidth=1),
				axis.ticks.length.y = unit(-0.5, "cm"), 
				strip.text = element_text(size = 8),
				panel.grid.major = element_line(color="#eeeeee", linewidth=0), 
				panel.grid.minor.x = element_line(color="#eeeeee", linewidth=0),
				panel.background = element_rect(fill="#ffffff"), 
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

format_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
		true = paste(month, years, sep = "\n"), 
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


pcr_df <- as.data.frame(read.table("../patchr/data/latest/watchdb.result.txt", header=TRUE, sep="\t", check.names = FALSE))
#pcr_df$report_day <- as_date(parse_date_time(pcr_df$collection_end_date, "%m/%d/%y"))
pcr_df$report_day <- as_date(mdy_hm(pcr_df$collection_end_date))
pcr_df$report_week <- floor_date(as_date(ymd(pcr_df$report_day)), unit="week")

# mean N2 over report_week, by facility:
wwtp_n2_df <- pcr_df %>% 
	filter(target_genetic_locus == "N2" & !is.na(target_copies_fn_per_cap)) %>% 
	group_by(location_id, report_week) %>% 
	summarise_at(vars(target_copies_fn_per_cap), list(mean_n2 = mean))

recent_wwtp_df <- subset(wwtp_n2_df, report_week > today() - months(9))

p_byfac_recent <- ggplot(recent_wwtp_df) + labs(y = "SARS-CoV-2 load per 100 persons", x = "") + 
	geom_point(aes(y=mean_n2, x=report_week), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = report_week, y = mean_n2), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = report_week, y = mean_n2), color = "#000000", na.rm = TRUE) + 
	#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
	facet_wrap(~location_id, nrow=4, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in WV Wastewater, as of ", max(recent_wwtp_df$report_week), sep=""), subtitle = "(weekly mean by facility)") + 
	scale_y_continuous() + 
	scale_x_date(breaks = "1 month", labels = format_dates_wide, limits=c(today() - months(9), today()))

# ggsave(paste("SARS.WV_ALL.", today(), ".png", sep=""), 
# 			 plot=p_byfac_recent, 
# 			 path="out/", 
# 			 width=1600, 
# 			 height=1200, 
# 			 units="px", 
# 			 dpi=72
# )

sentinel_df <- pcr_df %>% 
	filter(target_genetic_locus == "N2" & !is.na(target_copies_fn_per_cap)) %>% 
	filter(location_id == "StarCityWWTP-01" | location_id == "PrincetonWWTP-01" | location_id == "WheelingWWTP-01")

recent_sentinel_df <- subset(sentinel_df, report_week > today() - months(9))


p_sentinel_recent <- ggplot(recent_sentinel_df) + labs(y = "SARS-CoV-2 load per 100 persons", x = "") + 
	geom_point(aes(y=target_copies_fn_per_cap, x=report_week), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = report_week, y = target_copies_fn_per_cap), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = report_week, y = target_copies_fn_per_cap), color = "#000000", na.rm = TRUE) + 
	facet_wrap(~location_id, nrow=1, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in WV Wastewater sentinel sites, as of ", max(recent_sentinel_df$report_week), sep="")) + 
	scale_y_continuous() + 
	scale_x_date(breaks = "1 month", labels = format_dates_wide, limits=c(today() - months(9), today()))

# ggsave(paste("SARS.WV_SENTINELS.", today(), ".png", sep=""), 
# 			 plot=p_sentinel_recent, 
# 			 path="out/", 
# 			 width=1600, 
# 			 height=1200, 
# 			 units="px", 
# 			 dpi=72
# )

# mean N2 over report_week, all facilities:
sum_n2_df <- pcr_df %>% 
	filter(target_genetic_locus == "N2" & !is.na(target_copies_fn_per_cap) & target_copies_fn_per_cap < 500000) %>% 
	group_by(report_week) %>% 
	summarise_at(vars(target_copies_fn_per_cap), list(mean_n2 = mean))

end_date <- as.Date(c("2024-05-01"), "%Y-%m-%d")
start_date <- as.Date(c("2021-01-01"), "%Y-%m-%d")

mobile_lab_df <- subset(sum_n2_df, report_week > start_date)
mobile_lab_df <- subset(mobile_lab_df, report_week < end_date)

# look at just the last 9 months:
recent_pcr_df <- subset(sum_n2_df, report_week > today() - months(9))

p_recent_sum <- ggplot(recent_pcr_df) + labs(y = "SARS-CoV-2 load per 100 persons", x = "") + 
	geom_point(aes(y=mean_n2, x=report_week), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = report_week, y = mean_n2), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	geom_line(aes(x = report_week, y = mean_n2), color = "#000000", na.rm = TRUE) + 
	#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
	#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	theme(legend.position = "none") +
	ggtitle(paste("COVID in Wastewater for WV, as of ", max(sum_n2_df$report_week), sep=""), subtitle = "(weekly mean over all facilities)") + 
	scale_y_continuous() + 
	scale_x_date(breaks = "1 week", labels = format_dates)

p_mob <- ggplot(mobile_lab_df) + labs(y = "", x = "") + 
	geom_point(aes(y=mean_n2, x=report_week), color = "#EAAA00", shape = 16, size = 12, alpha=1.0, na.rm = TRUE) + 
	geom_line(aes(x = report_week, y = mean_n2), color = "#002855", linewidth = 8, na.rm = TRUE) + 
	#geom_col(aes(x = report_week, y = mean_n2), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	#	geom_ribbon(aes(x = day_end, y = n2.loadcap.day5.mean, xmin = min(day_end), xmax = max(day_end), ymin = n2.loadcap.day5.mean-n2.loadcap.day5.ci, ymax = n2.loadcap.day5.mean+n2.loadcap.day5.ci), fill = "#945200", alpha = 0.2) + 
	#	facet_wrap(~location_common_name, nrow=2, scales = "free") +
	#theme(legend.position = "none") +
	plot_theme() + 
	#ggtitle(paste("COVID in Wastewater for WV, as of ", max(sum_n2_df$report_week), sep=""), subtitle = "(weekly mean over all facilities)") + 
	scale_y_continuous() + 
	scale_x_date(breaks = "3 months", labels = format_dates)

# ggsave(paste("SARS.WV_SUM.", today(), ".png", sep=""), 
# 			 plot=p_recent_sum, 
# 			 path="out/", 
# 			 width=1600, 
# 			 height=1200, 
# 			 units="px", 
# 			 dpi=72
# 			)




