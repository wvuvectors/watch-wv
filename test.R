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

#library(plotly)
#library(shiny)
#library(shinythemes)
#library(shinyWidgets)
#library(leaflet)
#library(leaflet.extras)
#library(rsconnect)

#library(sf)
#library(tigris)
#options(tigris_use_cache = TRUE)

#rsconnect::deployApp('path/to/your/app')

L_per_gal=3.78541

Sys.setenv(TZ="America/New_York")
today <- Sys.Date()

first_day <- as_date("2021-10-03")
last_day <- today

blues7 <- c("1" = "#DEEBF7", "2" = "#C6DBEF", "3" = "#9ECAE1", "4" = "#6BAED6", "5" = "#4292C6", "6" = "#2171B5", "7" = "#084594") 
reds7 <- c("1" = "#FCBBA1", "2" = "#FC9272", "3" = "#FB6A4A", "4" = "#EF3B2C", "5" = "#CB181D", "6" = "#A50F15", "7" = "#67000D") 
greens7 <- c("1" = "#C7E9C0", "2" = "#A1D99B", "3" = "#74C476", "4" = "#41AB5D", "5" = "#238B45", "6" = "#006D2C", "7" = "#00441B") 
grays7 <- c("1" = "#D9D9D9", "2" = "#BDBDBD", "3" = "#969696", "4" = "#737373", "5" = "#525252", "6" = "#252525", "7" = "#000000") 

reds3 <- c("1" = "#FC9272", "2" = "#CB181D", "3" = "#67000D") 
blues2 <- c("1" = "#9ECAE1", "2" = "#2171B5") 

my_theme <- function () { 
	theme_bw() + theme(axis.text = element_text(size = 8),
										 axis.title = element_text(size = 9, color="#333333"),
										 axis.line.x = element_line(color="#bbbbbb", size=1),
										 axis.line.y = element_line(color="#bbbbbb", size=1),
										 axis.ticks.length.y = unit(-.5, "cm"), 
										 strip.text = element_text(size = 8),
										 panel.grid.major = element_line(color="#eeeeee", size=1), 
										 panel.grid.minor = element_blank(),
										 panel.background = element_blank(), 
										 legend.position = "bottom",
#										 panel.border = element_rect(fill=NA, color="#bbbbbb", size=1), 
										 panel.border = element_blank(), 
										 strip.background = element_rect(fill = 'white', color = 'white')
)}
  
#county_shapes <- tigris::counties(state = "WV", class="sf") %>% st_transform(4326)
#county_shapes <- read_sf("data/cb_2018_54_cousub_500k/cb_2018_54_cousub_500k.shp") %>% st_transform(4326)
#wvdnr_sf <- read_sf("https://services9.arcgis.com/SQbkdxLkuQJuLGtx/ArcGIS/rest/services/WVDNR_Districts/FeatureServer/0")


watch_file = "WaTCH-WV/data/watch_dashboard.LATEST.txt"
watch_df <- as.data.frame(read.table(watch_file, sep="\t", header=TRUE, check.names=FALSE))
watch_df <- watch_df %>% filter(status == "active" | status == "new")

#locations <- unique(watch_df$location_common_name)

# Handle dates
watch_df$"Sample Composite Start" <- mdy_hm(watch_df$"Sample Composite Start")
#watch_df$"Sample Composite End" <- mdy_hm(watch_df$"Sample Composite End")

watch_df$day <- as_date(watch_df$"Sample Composite Start")
watch_df$week_starting <- floor_date(watch_df$day, "week", week_start = 1)
watch_df$week_ending <- ceiling_date(watch_df$day, "week", week_start = 0)

watch_df$week_num <- week(watch_df$week_ending)
watch_df$week_alt <- 1 + (watch_df$week_num %% 2)


# Handle test result numbers
watch_df$n1 = watch_df$"Assay Target 1 Result (CN/L)"
watch_df$n2 = watch_df$"Assay Target 2 Result (CN/L)"
watch_df$daily_flow = watch_df$"Sample Flow (MGD)"

watch_df %>% mutate(n1 = replace_na(n1, 0))
watch_df %>% mutate(n2 = replace_na(n2, 0))
watch_df %>% mutate(daily_flow = replace_na(daily_flow, 0))

watch_df$n1n2 = rowMeans(watch_df[,c("n1", "n2")], na.rm=TRUE)

watch_df <- watch_df %>% filter(day >= first_day & day <= last_day)

watch_df <- watch_df %>% mutate(daily_flow = replace_na(daily_flow, 0))
watch_df <- watch_df %>% mutate(n1n2 = replace_na(n1n2, 0))

watch_df <- watch_df %>%
	mutate(rnamass = (daily_flow*n1n2*L_per_gal),
				 rnamass_log10 = log10(rnamass))

watch_df[watch_df == -Inf] <- NA

ci95 <- function(x) {
	0.5 * qt((0.95 / 2), length(x) - 1) * (sd(x) / sqrt(length(x)))
#	m <- mean(x)
#	se <- sd(x) / sqrt(x.length)
#	ci = qt(1 - (0.05 / 2), x.length - 1) * se
#	lower.ci = m - ci
#	upper.ci = m + ci
}

ci99 <- function(x) {
	0.5 * qt(0.99, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

se <- function(x) {
	0.5 * sd(x) / sqrt(length(x))
}

rollwin=5

loc_df <- watch_df %>% filter(location_common_name == "Morgantown Star City")
loc_zoo <- zoo(loc_df$rnamass, loc_df$day)
loc_mean_zoo <- rollmean(loc_zoo, rollwin, fill=NA, align="right")
loc_mean_df <- fortify(loc_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.rnamass"))

loc_se_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = se)
loc_se_df <- fortify(loc_se_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.se"))

loc_ci95_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci95)
loc_ci95_df <- fortify(loc_ci95_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci95"))

loc_ci99_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci99)
loc_ci99_df <- fortify(loc_ci99_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci99"))

loc_join1_df <- left_join(loc_mean_df, loc_se_df, by = c("day" = "day"))
loc_join2_df <- left_join(loc_join1_df, loc_ci95_df, by = c("day" = "day"))
loc_join_df <- left_join(loc_join2_df, loc_ci99_df, by = c("day" = "day"))

loc_rolled_df <- left_join(loc_df, loc_join_df, by = c("day" = "day"))
loc_rolled_df <- select(loc_rolled_df, -c(Series.x, Series.y, Series.x.x, Series.y.y))

lims_x_date <- as.Date(strptime(c(first_day, last_day), format = "%Y-%m-%d"))

p1 <- ggplotly(ggplot(loc_rolled_df, aes(x = day, y = rollmean.rnamass)) + 
									#labs(y = "Mean Weekly Mass Load", x = "") + 
									labs(y = paste0(rollwin, "-Day Rolling Mean Mass Load", sep=""), x = "") + 
									geom_point(color="red", shape = 16, size = 1, alpha=0.5) + 
									geom_line(show.legend = FALSE, color = "#000000") + 
#									geom_errorbar(aes(ymin=rollmean.rnamass-rollmean.se, ymax=rollmean.rnamass+rollmean.se), width=.2, position=position_dodge(0.05)) + 
									geom_ribbon(aes(ymin=rollmean.rnamass-rollmean.ci95, ymax=rollmean.rnamass+rollmean.ci95), fill="yellow", alpha=0.2) + 
									geom_ribbon(aes(ymin=rollmean.rnamass-rollmean.ci99, ymax=rollmean.rnamass+rollmean.ci99), fill="orange", alpha=0.2) + 
#									scale_y_continuous(limits=c(0,NA)) + 
									scale_y_continuous() + 
									scale_x_date(breaks = "1 month", date_labels = '%b-%Y', limits = lims_x_date) + 
									my_theme() + 
									#scale_fill_manual(values = reds7) +
									#ggtitle("Weekly Mean COVID, All WV Treatment Facilities") + 
									theme(axis.text.x = element_text(angle = 60, hjust=0.9)))

