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
#library(leaflet.extras)
library(rsconnect)

library(sf)
library(tigris)
options(tigris_use_cache = TRUE)

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
										 axis.ticks.length.y = unit(-0.5, "cm"), 
										 strip.text = element_text(size = 8),
										 panel.grid.major = element_line(color="#eeeeee", size=1), 
										 panel.grid.minor.y = element_blank(),
										 panel.grid.minor.x = element_line(color="#eeeeee", size=1),
										 panel.background = element_blank(), 
										 legend.position = "bottom",
										 panel.border = element_rect(fill=NA, color="#bbbbbb", size=1), 
#										 panel.border = element_blank(), 
										 strip.background = element_rect(fill = 'white', color = 'white')
)}

ci90 <- function(x) {
	0.5 * qt(0.90, length(x) - 1) * (sd(x) / sqrt(length(x)))
#	m <- mean(x)
#	se <- sd(x) / sqrt(x.length)
#	ci = qt(1 - (0.05 / 2), x.length - 1) * se
#	lower.ci = m - ci
#	upper.ci = m + ci
}

ci95 <- function(x) {
	0.5 * qt(0.95, length(x) - 1) * (sd(x) / sqrt(length(x)))
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
df_watch <- as.data.frame(read.table(watch_file, sep="\t", header=TRUE, check.names=FALSE))
df_watch <- df_watch %>% filter(status == "active" | status == "new")


# Convert date strings into Date objects

df_watch$"Sample Composite Start" <- mdy_hm(df_watch$"Sample Composite Start")
#df_watch$"Sample Composite End" <- mdy_hm(df_watch$"Sample Composite End")

df_watch$day <- as_date(df_watch$"Sample Composite Start")
df_watch$week_starting <- floor_date(df_watch$day, "week", week_start = 1)
df_watch$week_ending <- ceiling_date(df_watch$day, "week", week_start = 0)

df_watch$week_num <- week(df_watch$week_ending)
df_watch$week_alt <- 1 + (df_watch$week_num %% 2)


# Set date constraints
# May want to change this in the future?
df_watch <- df_watch %>% filter(day >= first_day & day <= last_day)


# Make some simpler aliases for common numeric columns
df_watch$n1 = df_watch$"Assay Target 1 Result (CN/L)"
df_watch$n2 = df_watch$"Assay Target 2 Result (CN/L)"
df_watch$daily_flow = df_watch$"Sample Flow (MGD)"

# Clean up some N/A entries
df_watch <- df_watch %>% mutate(n1 = replace_na(n1, 0))
df_watch <- df_watch %>% mutate(n2 = replace_na(n2, 0))
df_watch <- df_watch %>% mutate(daily_flow = replace_na(daily_flow, 0))

df_watch$n1n2 = rowMeans(df_watch[,c("n1", "n2")], na.rm=TRUE)
df_watch <- df_watch %>% mutate(n1n2 = replace_na(n1n2, 0))

#df_watch <- df_watch %>% mutate(daily_flow = replace_na(daily_flow, 0))
df_wwtp <- df_watch %>% filter(category == "wwtp" & daily_flow > 0) %>%
					 mutate(rnamass = (daily_flow*n1n2*L_per_gal),
					 				rnamass.n1 = (daily_flow*n1*L_per_gal), 
					 				rnamass.n2 = (daily_flow*n2*L_per_gal)) %>% 
					 arrange(day)
df_wwtp[df_wwtp == -Inf] <- NA


#for (loc in unique(df_wwtp$location_common_name)) {
#	df_loc <- df_wwtp %>% filter(location_common_name)
#	zoo_loc <- zoo(df_loc$rnamass, df_loc$day)
#	zoo_ci90 <- rollapply(zoo_loc, width=rollwin, fill=NA, align="right", FUN = ci90)
#	df_ci90 <- fortify(zoo_ci90, melt=TRUE, names=c(Index="day", Value="rollmean.ci90"))
#}


# set  up some global numbers
FACILITY_TOTAL_WWTP = n_distinct(df_wwtp$location_common_name)
CAP_TOTAL_WWTP = sum(distinct(df_wwtp, location_common_name, capacity_mgd)$capacity_mgd)+1
POP_TOTAL_WWTP = sum(distinct(df_wwtp, location_common_name, population_served)$population_served)
CTY_TOTAL_WWTP = n_distinct(df_wwtp$counties_served)

SAMPLE_TOTAL_WWTP = n_distinct(df_wwtp$"Sample ID")
FIRST_DATE_WWTP = min(df_wwtp$day)
LAST_DATE_WWTP = max(df_wwtp$day)



df_upstream <- df_watch %>% filter(category == "upstream") %>% arrange(day)

  				
# set  up some global numbers
FACILITY_TOTAL_UPSTREAM = n_distinct(df_upstream$location_common_name)

SAMPLE_TOTAL_UPSTREAM = n_distinct(df_upstream$"Sample ID")
FIRST_DATE_UPSTREAM = min(df_upstream$day)
LAST_DATE_UPSTREAM = max(df_upstream$day)



#map_pal <- colorNumeric(
#	palette = c('gold', 'orange', 'dark orange', 'orange red', 'red', 'dark red'),
#	domain = df_wwtp_day$capacity_mgd)
  
