library(tidyverse)
library(dplyr)
library(lubridate)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)

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
				legend.position = "right",
				#legend.justification = c("left", "top"),
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

lineage_pal <- c(
	"0" = "red", 
	"1" = "blue", 
	"2" = "darkgreen", 
	"3" = "orange"
)


all_df <- as.data.frame(read.table("../seqr/data/latest/seqrdb.txt", header=TRUE, sep="\t", check.names = FALSE))
all_df$day_end <- as_date(mdy_hm(all_df$sample_collection_end_datetime))
all_df$day_start <- as_date(mdy_hm(all_df$sample_collection_start_datetime))
all_df$report_week <- floor_date(as_date(mdy_hm(all_df$sample_collection_end_datetime)), unit="week")
all_df$percent <- as.numeric(all_df$variant_proportion) * 100

all_df <- all_df %>% mutate(
	color_group = case_when(
		str_detect(lineage_group, "^B\\.") ~ "B*",
		str_detect(lineage_group, "^BA") ~ "BA*",
		str_detect(lineage_group, "^BA\\.2\\.86") ~ "BA.2.86",
		str_detect(lineage_group, "^XBB") ~ "XBB*",
		str_detect(lineage_group, "^EG") ~ "EG*",
		str_detect(lineage_group, "^HV") ~ "HV*",
		str_detect(lineage_group, "^JN") ~ "JN*"
	)
)
all_df$color_group <- replace_na(all_df$color_group, "Other XBB")


recent_df <- subset(all_df %>% filter(str_detect(location, "WWTP")), day_end > today() - days(120))

recent_2pct_df <- recent_df %>% filter(percent >= 2)

recent_10pct_df <- recent_df %>% filter(percent >= 10)

end_date <- as.Date(c("2024-05-01"), "%Y-%m-%d")
start_date <- as.Date(c("2021-01-01"), "%Y-%m-%d")

mobile_lab_df <- subset(all_df, report_week > start_date)
mobile_lab_df <- subset(mobile_lab_df, report_week < end_date)

mobile_lab_df2 <- mobile_lab_df %>% 
	group_by(report_week, location, color_group) %>% 
	summarise_at(vars(percent), list(weekly_mean_pct = mean))

mobile_lab_df3 <- mobile_lab_df %>% 
	group_by(report_week, color_group) %>% 
	summarise_at(vars(percent), list(weekly_mean_pct = mean))

#recent_weekly_var_df <- recent_var_df %>% filter(percent > 5) %>% 
#						 group_by(week, parental) %>% 
#						 summarize(MUT = mean(percent, na.rm = TRUE)) %>% 
#						 arrange(week, parental)


# p_all <- ggplot(mobile_lab_df, aes(fill=lineage_group, y=percent, x=day_end)) + 
# 	geom_bar(position="stack", stat="identity", aes(fill=factor(color_group))) + 
# 	scale_fill_brewer(type="qual", palette = "Dark2") + 
# 	facet_wrap(~location, nrow=3) + 
# 	labs(x="", y="") + 
# 	#ggtitle("Proportion of SARS-CoV-2 Lineages in WV Wastewater, by Location") +
# 	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
# 	theme(legend.position = "bottom", legend.title=element_blank())
# 
# 
# p_sentinel <- ggplot(recent_df %>% filter(str_detect(location, "StarCity") | str_detect(location, "Princeton") | str_detect(location, "Wheeling")), aes(fill=lineage_group, y=percent, x=day_end)) + 
# 	geom_bar(position="stack", stat="identity", aes(fill=factor(color_group))) + 
# 	scale_fill_brewer(type="qual", palette = "Dark2") + 
# 	facet_wrap(~location, nrow=3) + 
# 	labs(x="", y="") + 
# 	ggtitle("Proportion of SARS-CoV-2 Lineages in WV Wastewater, by Location") +
# 	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
# 	theme(legend.position = "bottom", legend.title=element_blank())
# 
# 
# p_all_2 <- ggplot(recent_2pct_df, aes(fill=lineage_group, y=percent, x=day_end)) + 
# 	geom_bar(position="stack", stat="identity") + 
# 	facet_wrap(~location, nrow=3) + 
# 	labs(x="", y="") + 
# 	ggtitle("Proportion of SARS-CoV-2 Lineages in WV Wastewater, by Location", subtitle="Only variants > 2% abundance are shown.") +
# 	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
# 	theme(legend.position = "bottom", legend.title=element_blank())
# 
# p_all_10 <- ggplot(recent_10pct_df, aes(fill=lineage_group, y=percent, x=day_end)) + 
# 	geom_bar(position="stack", stat="identity") + 
# 	facet_wrap(~location, nrow=3) + 
# 	labs(x="", y="") + 
# 	ggtitle("Proportion of SARS-CoV-2 Lineages in WV Wastewater, by Location", subtitle="Only variants > 10% abundance are shown.") +
# 	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
# 	theme(legend.position = "bottom", legend.title=element_blank())
# 
p_starcity <- ggplot(mobile_lab_df2 %>% filter(location == "StarCityWWTP-01"), aes(fill=color_group, y=weekly_mean_pct, x=report_week)) +
	geom_bar(position="stack", stat="identity") +
	labs(x="", y="") +
	plot_theme() + 
	#ggtitle("Proportion of SARS-CoV-2 Lineages in the Star City Wastewater Facility", subtitle="Only variants > 2% abundance are shown.") +
	scale_x_date(date_breaks = "3 months", labels = format_dates)
	#theme(legend.position = "bottom", legend.title=element_blank())


p_sum <- ggplot(mobile_lab_df3, aes(fill=color_group, y=weekly_mean_pct, x=report_week)) +
	geom_bar(position="stack", stat="identity") +
	labs(x="", y="") +
	plot_theme() + 
	#ggtitle("Proportion of SARS-CoV-2 Lineages in the Star City Wastewater Facility", subtitle="Only variants > 2% abundance are shown.") +
	scale_x_date(date_breaks = "3 months", labels = format_dates)
#theme(legend.position = "bottom", legend.title=element_blank())


# p_starcity_10 <- ggplot(recent_10pct_df %>% filter(location == "StarCityWWTP-01"), aes(fill=lineage_group, y=percent, x=day_end)) + 
# 	geom_bar(position="stack", stat="identity") + 
# 	labs(x="", y="") + 
# 	ggtitle("Proportion of SARS-CoV-2 Lineages in the Star City Wastewater Facility", subtitle="Only variants > 10% abundance are shown.") +
# 	scale_x_date(date_breaks = "5 days", labels = format_dates) + 
# 	theme(legend.position = "bottom", legend.title=element_blank())
# 

p_var_stadium <- ggplot(all_df %>% filter(str_detect(location, "Stadium")), aes(fill=lineage_group, y=percent, x=facility)) + 
	geom_bar(position="stack", stat="identity") + 
	facet_grid(~day_start) + 
	labs(x="", y="") + 
	ggtitle("Proportion of SARS-CoV-2 Lineages during WVU Home Football Games", subtitle="Only variants > 2% abundance are shown.") +
	scale_x_date(date_breaks = "1 day", labels = format_dates) + 
	theme(legend.position = "bottom", legend.title=element_blank())





focus_df <- subset(all_df %>% filter(str_detect(location, "WheelingWWTP") & percent > 1.0 & day_end > today() - months(9) & day_end < today() - months(3)))

#p_focus <- ggplot(focus_df %>% filter(str_detect(lineage_group, "EG.5")), aes(fill=variant, y=percent, x=day_end)) + 
p_focus <- ggplot(focus_df, aes(fill=lineage_group, y=percent, x=day_end)) + 
	geom_bar(position="stack", stat="identity") + 
#	scale_fill_brewer(type="qual", palette = "Dark2") + 
#	scale_color_manual(values = lineage_pal) + 
	labs(x="", y="") + 
	ggtitle("Proportion of SARS-CoV-2 Lineages in Wheeling WWTP") +
	scale_x_date(date_breaks = "7 days", labels = format_dates) + 
	theme(legend.position = "right", legend.title=element_blank())


