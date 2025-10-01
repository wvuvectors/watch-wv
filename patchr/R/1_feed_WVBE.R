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

#TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)", "Human Norovirus GII (HuNoV-GII)", "Human Norovirus GI (HuNoV-GI)")
TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)")
GENLOCI <- c("SC2", "M", "NEP/NS1", "N2")

DB_BASE <- "../dashboard/data"
DB_RESULTS_WVU <- paste(DB_BASE, "/watchdb.result.txt", sep="")
DB_SAMPLES_WVU <- paste(DB_BASE, "/watchdb.sample.txt", sep="")

DB_RESULTS_MU <- paste(DB_BASE, "/mu.result.txt", sep="")
DB_SAMPLES_MU <- paste(DB_BASE, "/mu.sample.txt", sep="")

DB_ALERTS <- paste(DB_BASE, "/wvdash.alerts.txt", sep="")

RES_BASE <- "resources"
RES_ALL <- paste(RES_BASE, "/watchdb.all_tables.xlsx", sep="")


Sys.setenv(TZ="America/New_York")
today <- Sys.Date()
#today <- as.Date("2022-07-12")


# Load data files
df_pcr_wvu <- as.data.frame(read.table(DB_RESULTS_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr_mu <- as.data.frame(read.table(DB_RESULTS_MU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr <- rbind(df_pcr_wvu, df_pcr_mu)

df_s_wvu <- as.data.frame(read.table(DB_SAMPLES_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_s_mu <- as.data.frame(read.table(DB_SAMPLES_MU, sep="\t", header=TRUE, check.names=FALSE))
df_sample <- rbind(df_s_wvu, df_s_mu)

df_alerts <- as.data.frame(read.table(DB_ALERTS, sep="\t", header=TRUE, check.names=FALSE))


df_res_loc <- as.data.frame(read_excel(RES_ALL, sheet = "location"))
df_res_county <- as.data.frame(read_excel(RES_ALL, sheet = "county"))
df_res_county <- df_res_county %>% select(location_counties_served = county_id, county_fips)
df_res_loc <- left_join(df_res_loc, df_res_county, by = "location_counties_served")

df_active_loc <- df_res_loc %>% filter(tolower(location_status) == "active")


# Restrict results to those that pass sample QC and have NTC below threshold
#
df_pcr <- df_pcr %>% filter(location_id %in% df_active_loc$location_id)
df_pcr <- df_pcr %>% filter(tolower(sample_qc) == "pass")
df_pcr <- df_pcr %>% filter(tolower(target_result_validated) != "ntc above threshold")
df_pcr <- df_pcr %>% filter(target %in% TARGETS)
df_pcr <- df_pcr %>% filter(target_genetic_locus %in% GENLOCI)

#df_pcr <- df_pcr %>% filter(location_id != "CheatLakeWWTP-01")

# Convert date strings into Date objects.
#
df_pcr$collection_start_datetime <- mdy_hm(df_pcr$collection_start_datetime)
df_pcr$collection_end_datetime <- mdy_hm(df_pcr$collection_end_datetime)

df_sample$sample_collection_start_datetime <- mdy_hm(df_sample$sample_collection_start_datetime)
df_sample$sample_collection_end_datetime <- mdy_hm(df_sample$sample_collection_end_datetime)
df_sample$sample_recovered_datetime <- mdy_hm(df_sample$sample_recovered_datetime)
df_sample$sample_received_date <- mdy(df_sample$sample_received_date)


# Create a date to use as the primary reference and make sure all tables use the same date.
#
df_pcr$date_primary <- as.Date(df_pcr$collection_end_datetime)


# Add some date objects, for use in summary stats.
#
df_pcr$year <- lubridate::year(df_pcr$date_primary)
df_pcr$month <- lubridate::month(df_pcr$date_primary)
df_pcr$epi_week <- lubridate::epiweek(df_pcr$date_primary)
df_pcr <- df_pcr %>%
	mutate(epi_year = case_when(
		df_pcr$month == 12 & df_pcr$epi_week == 1 ~ df_pcr$year + 1,
		TRUE ~ df_pcr$year
	))
	
#df_pcr$mmr_year <- lubridate::year(df_pcr$date_primary)
#df_pcr$mmr_week <- lubridate::week(df_pcr$date_primary)
#df_pcr$week_starting <- floor_date(df_pcr$date_primary, "week", week_start = 1)
#df_pcr$week_ending <- ceiling_date(df_pcr$date_primary, "week", week_start = 1)


#################################################

df_pcr$target_copies_fn_per_cap <- df_pcr$target_copies_fn_per_cap/df_pcr$target_per_capita_basis
df_pcr$target_copies_per_ldcap <- df_pcr$target_copies_per_ldcap/df_pcr$target_per_capita_basis
df_pcr$target_per_capita_basis <- 1

df_rs <- df_pcr %>% 
	filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow)) %>% 
	select(collection_end_datetime, assay_id, location_id, target, total_abundance = target_copies_flownorm, abundance_per_person = target_copies_fn_per_cap, epi_year, epi_week)

# average over epi week and facility
df_wkly <- df_rs %>% 
	group_by(location_id, target, epi_year, epi_week) %>% 
	summarise(abundance_per_person = mean(abundance_per_person, na.rm=TRUE),
						total_abundance = mean(total_abundance, na.rm=TRUE))

# Add county names
df_counties <- df_active_loc %>% select(location_id, county = location_counties_served, county_fips)
df_wkly <- left_join(df_wkly, df_counties, by = "location_id")


df_agg <- df_wkly %>% 
	group_by(county, target, epi_year, epi_week) %>% 
	summarise(abundance_per_person = sum(abundance_per_person, na.rm=TRUE),
						total_abundance = sum(total_abundance, na.rm=TRUE))
df_agg <- left_join(df_agg, unique(df_counties %>% select(county, county_fips)), by = "county")


write.table(df_agg, file = stdout(), sep = "\t", row.names = FALSE)


