#! /usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(data.table)
library(DT)
library(zoo)
library(rlang)
library(glue)
library(readxl)
library(scales)
library(lubridate)

# county
# target
# epi_year
# epi_week
# abundance_per_person
# total_abundance
# county_fips
# abundance_alert_level
# trend_alert_level

# epi_date
# mean_abundance
# location_id
# target
# cal_year
# cal_month
# cal_week
# epi_year
# epi_week
# rolling_mean
# abundance_fold_change
# abundance_level
# trend_slope
# trend_level
# spike_slope

TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)")
GENLOCI <- c("SC2", "N2", "M", "NEP/NS1", "G")

ALEVEL_STRINGS <- c("UNKNOWN", "VERY LOW", "LOW", "MODERATE", "HIGH", "VERY HIGH")
TLEVEL_STRINGS <- c("UNKNOWN", "SURGING", "DECR. RAPIDLY", "DECREASING", "STABLE", "INCREASING", "INCR. RAPIDLY")


DB_BASE <- "../dashboard/data"
DB_RESULTS <- paste(DB_BASE, "/wvdash.rsstable.txt", sep="")
RES_ALL <- paste(DB_BASE, "/watchdb.all_tables.xlsx", sep="")


Sys.setenv(TZ="America/New_York")
today <- Sys.Date()
#today <- as.Date("2022-07-12")


# Load data files
df_rss <- as.data.frame(read.table(DB_RESULTS, sep="\t", header=TRUE, check.names=FALSE))
df_res_county <- as.data.frame(read_excel(RES_ALL, sheet = "county"))
df_res_county <- df_res_county %>% select(county=county_id, county_fips)

df_agg <- df_rss %>% 
	filter(location_id %in% df_res_county$county) %>% 
	select(
		county=location_id, 
		target, 
		epi_year, 
		epi_week, 
		abundance_per_person=mean_abundance, 
		abundance_alert_level=abundance_level, 
		trend_alert_level=trend_level) %>%
	mutate(
		abundance_alert_level = ALEVEL_STRINGS[as.integer(abundance_alert_level)], 
		trend_alert_level = TLEVEL_STRINGS[as.integer(trend_alert_level)]
	)

df_agg <- left_join(df_agg, df_res_county, by = "county")

write.table(df_agg, file = stdout(), sep = "\t", row.names = FALSE)


