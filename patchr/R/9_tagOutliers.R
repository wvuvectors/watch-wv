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

f <- file("stdin")
open(f)
while(length(line <- readLines(f, n = 1)) > 0) {
  input_fnames <- strsplit(line, " ")[[1]]
}
close(f)

WVU_RESULTS_F <- input_fnames[1]
MU_RESULTS_F <- input_fnames[2]
ALL_RESOURCE_F <- input_fnames[3]

# WVU_RESULTS_F <- "../patchr/data/latest/watchdb.result.txt"
# MU_RESULTS_F <- "../dashboard/data/mu.result.txt"
# ALL_RESOURCE_F <- "../patchr/resources/watchdb.all_tables.xlsx"


WVD_BASE <- "../dashboard/data"
WVD_THRESHOLDS_F <- paste(WVD_BASE, "/wvdash.thresholds.txt", sep="")

# Load data files.
df_pcr_wvu <- as.data.frame(read.table(WVU_RESULTS_F, sep="\t", header=TRUE, check.names=FALSE))
df_pcr_mu <- as.data.frame(read.table(MU_RESULTS_F, sep="\t", header=TRUE, check.names=FALSE))
df_pcr <- rbind(df_pcr_wvu, df_pcr_mu)

df_res_loc <- as.data.frame(read_excel(ALL_RESOURCE_F, sheet = "location"))
df_active_loc <- df_res_loc %>% filter(tolower(location_status) == "active")

df_thresholds <- as.data.frame(read.table(WVD_THRESHOLDS_F, sep="\t", header=TRUE, check.names=FALSE))

#
# Dates, dates, dates.
#

# Convert date strings into Date objects.
df_pcr$collection_start_datetime <- mdy_hm(df_pcr$collection_start_datetime)
df_pcr$collection_end_datetime <- mdy_hm(df_pcr$collection_end_datetime)

# Create a date to use as the primary reference and make sure all tables use the same date.
df_pcr$date_primary <- as.Date(df_pcr$collection_end_datetime)

# Add epi week info.
df_pcr$year <- lubridate::year(df_pcr$date_primary)
df_pcr$month <- lubridate::month(df_pcr$date_primary)
df_pcr$epi_week <- lubridate::epiweek(df_pcr$date_primary)
df_pcr <- df_pcr %>%
  mutate(epi_year = case_when(
    df_pcr$month == 12 & df_pcr$epi_week == 1 ~ df_pcr$year + 1,
    TRUE ~ df_pcr$year
  ))

df_pcr <- df_pcr %>% mutate(epi_date = as.integer(epi_year) + as.integer(epi_week)/100)

#
# To keep numbers small and more digestable, make sure all samples are per capita
#
df_pcr$target_copies_fn_per_cap <- df_pcr$target_copies_fn_per_cap/df_pcr$target_per_capita_basis
df_pcr$target_copies_per_ldcap <- df_pcr$target_copies_per_ldcap/df_pcr$target_per_capita_basis
df_pcr$target_per_capita_basis <- 1

#
# Data cleaning.
#

# Restrict results to those that pass sample QC and have NTC below threshold.
df_pcr <- df_pcr %>% filter(tolower(sample_qc) == "pass")
df_pcr <- df_pcr %>% filter(tolower(target_result_validated) != "ntc above threshold")

# Restrict results to those included in the primary target and locus lists.
df_pcr <- df_pcr %>% filter(target %in% TARGETS)
df_pcr <- df_pcr %>% filter(target_genetic_locus %in% GENLOCI)

# Restrict results to those from active locations only. Also store these location ids in 
# an unnamed array for later use.
df_pcr <- df_pcr %>% filter(location_id %in% df_active_loc$location_id)
LOCATIONS <- unname(unique(df_pcr$location_id))

df_agg <- data.frame()

for (i in 1:length(TARGETS)) {
	targ <- TARGETS[i]
	for (j in 1:length(LOCATIONS)) {
		loc <- LOCATIONS[j]
		df_this <- df_pcr %>% filter(target == targ & location_id == loc)
#		df_added <- df_this %>% identify_outliers(target_copies_fn_per_cap)
#		df_added <- df_added %>% select(assay_id, is.outlier, is.extreme)

    mean_norm <- mean(df_this$target_copies_fn_per_cap, na.rm=TRUE)
    sd_norm <- sd(df_this$target_copies_fn_per_cap, na.rm=TRUE)
    
    mean_raw <- mean(df_this$target_copies_per_ld, na.rm=TRUE)
    sd_raw <- sd(df_this$target_copies_per_ld, na.rm=TRUE)
    
    df_this <- df_this %>%
		  mutate(z = if_else(
		    tolower(event_type) == "routine surveillance" & !is.na(sample_flow),
		    (target_copies_fn_per_cap - mean_norm) / sd_norm,
		    (target_copies_per_ld - mean_raw) / sd_raw)
		  )

    df_z <- df_thresholds %>% 
      filter(category == "z" & target == targ) %>% 
      select(step)
    z_step <- as.numeric(df_z$step)
    
    df_this$is.outlier <- (df_this$z > (z_step * 1))
		df_this$is.extreme <- (df_this$z > (z_step * 2))

		#low <- median(df_this$target_copies_fn_per_cap) - 3 * mad(df_this$target_copies_fn_per_cap, constant=1)
		#high <- median(df_this$target_copies_fn_per_cap) + 3 * mad(df_this$target_copies_fn_per_cap, constant=1)
		#subset dataframe where points value is outside of low and high bounds
		#df_this$hampel <- (df_this$target_copies_fn_per_cap<low | df_this$target_copies_fn_per_cap>high)

		df_added <- df_this %>% select(assay_id, z, is.outlier, is.extreme)
		df_agg <- rbind(df_agg, df_added)
	}
}

df_pcr_wvu_out <- left_join(df_pcr_wvu, df_agg, by="assay_id")
df_pcr_mu_out <- left_join(df_pcr_mu, df_agg, by="assay_id")

write.table(df_pcr_wvu_out, file = paste(WVD_BASE, "/wvu.result_tagged.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_pcr_mu_out, file = paste(WVD_BASE, "/mu.result_tagged.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)


# df_all <- left_join(df_pcr, df_agg, by="assay_id")
# 
# p1 <- ggplot(df_all %>% filter(target == TARGETS[1]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-5, 5)) + 
#   labs(y="count", title=paste0(TARGETS[1], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[1], ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p2 <- ggplot(df_all %>% filter(target == TARGETS[2]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-1, 2)) + 
#   labs(y="count", title=paste0(TARGETS[2], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[2], ".png", sep=""), plot = p2, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p3 <- ggplot(df_all %>% filter(target == TARGETS[3]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-2, 2)) + 
#   labs(y="count", title=paste0(TARGETS[3], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[3], ".png", sep=""), plot = p3, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p4 <- ggplot(df_all %>% filter(target == TARGETS[4]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-2, 3)) +
#   labs(y="count", title=paste0(TARGETS[4], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[4], ".png", sep=""), plot = p4, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 

