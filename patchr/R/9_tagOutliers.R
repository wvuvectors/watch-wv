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
library(rstatix)

TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)")
GENLOCI <- c("SC2", "N2", "M", "NEP/NS1", "G")

DB_BASE <- "../dashboard/data"
DB_RESULTS_WVU <- paste(DB_BASE, "/watchdb.result.txt", sep="")

DB_RESULTS_MU <- paste(DB_BASE, "/mu.result.txt", sep="")

RES_ALL <- paste(DB_BASE, "/watchdb.all_tables.xlsx", sep="")


# Load data files
df_pcr_wvu <- as.data.frame(read.table(DB_RESULTS_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr_mu <- as.data.frame(read.table(DB_RESULTS_MU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr <- rbind(df_pcr_wvu, df_pcr_mu)

df_res_loc <- as.data.frame(read_excel(RES_ALL, sheet = "location"))
df_active_loc <- df_res_loc %>% filter(tolower(location_status) == "active")
LOCATIONS <- unname(df_active_loc$location_id)

# Restrict results to those that pass sample QC and have NTC below threshold
#
df_pcr <- df_pcr %>% filter(location_id %in% df_active_loc$location_id)
df_pcr <- df_pcr %>% filter(tolower(sample_qc) == "pass")
df_pcr <- df_pcr %>% filter(tolower(target_result_validated) != "ntc above threshold")
df_pcr <- df_pcr %>% filter(target %in% TARGETS)
df_pcr <- df_pcr %>% filter(target_genetic_locus %in% GENLOCI)

# Convert date strings into Date objects.
#
df_pcr$collection_start_datetime <- mdy_hm(df_pcr$collection_start_datetime)
df_pcr$collection_end_datetime <- mdy_hm(df_pcr$collection_end_datetime)


df_rs <- df_pcr %>% filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow))

df_agg <- data.frame()

for (i in 1:length(TARGETS)) {
	targ <- TARGETS[i]
	for (j in 1:length(LOCATIONS)) {
		loc <- LOCATIONS[j]
		df_this <- df_rs %>% filter(target == targ & location_id == loc)
#		df_added <- df_this %>% identify_outliers(target_copies_fn_per_cap)
#		df_added <- df_added %>% select(assay_id, is.outlier, is.extreme)

		df_this$z <- (df_this$target_copies_fn_per_cap-mean(df_this$target_copies_fn_per_cap))/sd(df_this$target_copies_fn_per_cap)
		#subset data frame where z-score of points value is greater than 3
		df_this$z_out_5 <- (df_this$z>5)
		df_this$z_out_10 <- (df_this$z>10)
		df_this$is.extreme <- (df_this$z>10)
		#outliers <- df[df$z>3, ]

		#low <- median(df_this$target_copies_fn_per_cap) - 3 * mad(df_this$target_copies_fn_per_cap, constant=1)
		#high <- median(df_this$target_copies_fn_per_cap) + 3 * mad(df_this$target_copies_fn_per_cap, constant=1)
		#subset dataframe where points value is outside of low and high bounds
		#df_this$hampel <- (df_this$target_copies_fn_per_cap<low | df_this$target_copies_fn_per_cap>high)

		df_added <- df_this %>% select(assay_id, is.extreme)
		df_agg <- rbind(df_agg, df_added)
	}
}

df_pcr_wvu_out <- left_join(df_pcr_wvu, df_agg, by="assay_id")
df_pcr_mu_out <- left_join(df_pcr_mu, df_agg, by="assay_id")

write.table(df_pcr_wvu_out, file = paste(DB_BASE, "/watchdb.result.tagged.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pcr_mu_out, file = paste(DB_BASE, "/mu.result.tagged.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)


