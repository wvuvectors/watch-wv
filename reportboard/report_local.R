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


format_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
		true = paste(day, month, years, sep = "\n"), 
		false = if_else(is.na(lag(month)) | lag(month) != month, true = paste(day, month, sep = "\n"), false = day)
	)
}

df_pcr <- as.data.frame(read.table("../patchr/data/latest/watchdb.result.txt", header=TRUE, sep="\t", check.names = FALSE))
df_nov <- df_pcr %>% filter(target_genetic_locus == "ORF1_2")

df_nov$collection_start_datetime <- mdy_hm(df_nov$collection_start_datetime)
df_nov$collection_end_datetime <- mdy_hm(df_nov$collection_end_datetime)
df_nov$date_primary <- as.Date(df_nov$collection_end_datetime)

df_nov$target_copies_fn_per_cap <- df_nov$target_copies_fn_per_cap/df_nov$target_per_capita_basis
df_nov$target_copies_per_ldcap <- df_nov$target_copies_per_ldcap/df_nov$target_per_capita_basis
df_nov$target_per_capita_basis <- 1

df_nov$val <- df_nov$target_copies_fn_per_cap

p1 <- ggplot(df_nov) + labs(y = "viral particles per person", x = "") + 
	geom_point(aes(x=date_primary, y=val), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = date_primary, y = val), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	facet_wrap(~location_id, nrow=4, scales = "free") +
	plot_theme() + 
	theme(legend.position = "none") +
	ggtitle(paste("Norovirus (GII) in WV Wastewater (", max(df_nov$date_primary), ")", sep="")) + 
#	scale_y_continuous(labels = comma, limits = c(0, MAX_LOADCAP)) + 
	scale_x_date(breaks = "1 day", labels = format_dates)

p2 <- ggplot(df_nov %>% filter(location_id == "StarCityWWTP-01" | location_id == "WheelingWWTP-01")) + labs(y = "viral particles per person", x = "") + 
	geom_point(aes(x=date_primary, y=val), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = date_primary, y = val), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	facet_wrap(~location_id, nrow=1, scales = "free") +
	plot_theme() + 
	theme(legend.position = "none") +
	ggtitle(paste("Norovirus (GII) in WV Wastewater (", max(df_nov$date_primary), ")", sep="")) + 
	#	scale_y_continuous(labels = comma, limits = c(0, MAX_LOADCAP)) + 
	scale_x_date(breaks = "1 day", labels = format_dates)

df_mean <- df_nov %>% 
	group_by(date_primary) %>%
	summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE)) %>%
	arrange(date_primary)

	
p3 <- ggplot(df_mean) + labs(y = "viral particles per person", x = "") + 
	geom_point(aes(x=date_primary, y=val), fill = "#945200", shape = 16, size = 2, alpha=0.3, na.rm = TRUE) + 
	geom_col(aes(x = date_primary, y = val), fill = "#945200", alpha = 0.7, na.rm = TRUE) + 
	#facet_wrap(~location_id, nrow=4, scales = "free") +
	plot_theme() +
	#theme(legend.position = "none") +
	ggtitle(paste("Norovirus (GII) in WV Wastewater (", max(df_nov$date_primary), ")", sep="")) + 
	#	scale_y_continuous(labels = comma, limits = c(0, MAX_LOADCAP)) + 
	scale_x_date(breaks = "1 day", labels = format_dates)
