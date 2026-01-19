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

get_LTcoefs <- function(window_data) {
  # window_data is a data frame slice
  model <- lm(mean_abundance ~ epi_date, data = as.data.frame(window_data))
  return(coef(model))
}

ABUND_ALERT_WINDOW <- 12			# Number of consecutive samples to use in determining abundance alert level.
TREND_ALERT_WINDOW <- 4				# Number of consecutive samples to use in trend calculations.
ABUND_RESPONSE_THRESHOLD <- 0.05	# Mean abundance values below this threshold (fraction of the rolling mean) are not used to fire alerts.

TARGETS <- c("SARS-CoV-2", "Influenza Virus A (FluA)", "Influenza Virus B (FluB)", "Respiratory Syncitial Virus, Human (RSV)")
GENLOCI <- c("SC2", "N2", "M", "NEP/NS1", "G")

WVD_BASE <- "../dashboard/data"
WVU_RESULTS_F <- paste(WVD_BASE, "/wvu.result_tagged.txt", sep="")
MU_RESULTS_F <- paste(WVD_BASE, "/mu.result_tagged.txt", sep="")
ALL_RESOURCE_F <- paste(WVD_BASE, "/watchdb.all_tables.xlsx", sep="")
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
df_pcr$date_of_collection <- as.Date(df_pcr$collection_end_datetime)

# Add epi week info.
df_pcr$year <- lubridate::year(df_pcr$date_of_collection)
df_pcr$month <- lubridate::month(df_pcr$date_of_collection)
df_pcr$week <- lubridate::week(df_pcr$date_of_collection)
df_pcr$epi_week <- lubridate::epiweek(df_pcr$date_of_collection)
df_pcr <- df_pcr %>%
  mutate(epi_year = case_when(
    df_pcr$month == 12 & df_pcr$epi_week == 1 ~ df_pcr$year + 1,
    df_pcr$month == 1 & df_pcr$epi_week == 53 ~ df_pcr$year - 1,
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

# skip any extreme values
df_pcr <- df_pcr %>% filter(is.extreme != TRUE)
# skip any outlier values
df_pcr <- df_pcr %>% filter(is.outlier != TRUE)

# Only want to process routine surveillance data. Alerts rely on consecutive samples, and 
# mass-gathering events don't fit this model.
df_rs <- df_pcr %>% filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow))

# Restrict results to those from active locations only. Also store these location ids in 
# an unnamed array for later use.
df_rs <- df_rs %>% filter(location_id %in% df_active_loc$location_id)
LOCATIONS <- unname(unique(df_rs$location_id))

df_counties <- df_active_loc %>% filter(location_id %in% LOCATIONS) %>% select(location_id, location_counties_served)
COUNTIES <- unname(unique(df_counties$location_counties_served))

# We're going to parse out sub-tables by location and target to calculate alerts.
# Then we'll rowbind them all back together into this dataframe, df_agg.
df_agg <- data.frame()

# First calculate by location (facility)
#
for (i in 1:length(LOCATIONS)) {
	loc <- LOCATIONS[i]
	for (j in 1:length(TARGETS)) {
		targ <- TARGETS[j]

		df_this <- df_rs %>% 
			filter(target == targ & location_id == loc) %>% 
			group_by(epi_date) %>% 
			reframe(mean_abundance = mean(target_copies_fn_per_cap, na.rm = TRUE),
								location_id = first(location_id), 
								target = first(target),
								cal_year = first(year), 
								cal_month = first(month), 
								cal_week = first(week), 
								epi_year = first(epi_year), 
								epi_week = first(epi_week))
		
		if (nrow(df_this) < 1) {
		  next
		}
		
		#print(paste0(loc, ": ", targ, ". MAX = ", max_abund, sep=""))

		df_added <- df_this %>% 
			arrange(epi_date) %>% 
			mutate(
				rolling_mean = rollmean(mean_abundance, k = ABUND_ALERT_WINDOW, fill = NA, align = "right")
			)

		df_added <- df_added %>% 
			mutate(
				abundance_fold_change = case_when(
					is.na(rolling_mean) ~ NaN,
					(mean_abundance <= ABUND_RESPONSE_THRESHOLD*rolling_mean) ~ NA,
					.default =  (mean_abundance - rolling_mean) / rolling_mean
				)
			)

		df_added <- df_added %>% 
			mutate(
				abundance_level = case_when(
					is.nan(abundance_fold_change) ~ 1,	# rolling mean is missing
					is.na(abundance_fold_change) ~ 2,		# rolling mean is too small to be useful
					(abundance_fold_change <= -0.5) ~ 3,
					(abundance_fold_change > -0.5 & abundance_fold_change <= 0.5) ~ 4,
					(abundance_fold_change > 0.5 & abundance_fold_change <= 1.5) ~ 5,
					(abundance_fold_change > 1.5) ~ 6
				)
			)

		trend_results <- rollapplyr(
			df_added,
			width = TREND_ALERT_WINDOW,
			FUN = get_LTcoefs,
			by.column = FALSE, # Apply the function to the entire window data frame
			fill = NA            # Fill initial empty values with NA
		)
		# Convert results to a data frame for easier use
		trend_coefs_df <- as.data.frame(trend_results)
		names(trend_coefs_df) <- c("trend_intercept", "trend_slope")
		# Add the results back to the original data frame
		df_added <- cbind(df_added, trend_slope = trend_coefs_df$trend_slope)
    
		trend_step <- df_thresholds %>% 
		  filter(category == "trend" & target == targ) %>% 
		  select(step)
		trend_step <- as.numeric(trend_step$step)
		
		df_added <- df_added %>%
			mutate(
				trend_level = case_when(
					is.na(trend_slope) ~ 1,
					(trend_slope <= -1.5*trend_step) ~ 3,
					(trend_slope > -1.5*trend_step & trend_slope <= -0.5*trend_step) ~ 4,
					(trend_slope > -0.5*trend_step & trend_slope <= 0.5*trend_step) ~ 5,
					(trend_slope > 0.5*trend_step & trend_slope <= 1.5*trend_step) ~ 6,
					(trend_slope > 1.5*trend_step) ~ 7
				)
			)
		
		# A drastic increase in mean abundance from one sample to the next is categorized as 
		# a "spike" and is handled differently in the alert system. To determine this, we 
		# calculate the trend line between each 2 consecutive samples.
		spike_results <- rollapplyr(
			df_added,
			width = 2,
			FUN = get_LTcoefs,
			by.column = FALSE, # Apply the function to the entire window data frame
			fill = NA            # Fill initial empty values with NA
		)
		# Convert results to a data frame for easier use
		spike_coefs_df <- as.data.frame(spike_results)
		names(spike_coefs_df) <- c("spike_intercept", "spike_slope")
		# Add the results back to the original data frame
		df_added <- cbind(df_added, spike_slope = spike_coefs_df$spike_slope)

		spike_threshold <- df_thresholds %>% 
		  filter(category == "spike" & target == targ) %>% 
		  select(step)
		spike_threshold <- as.numeric(spike_threshold$step)
		
		df_added <- df_added %>%
			mutate(
				trend_level = if_else(spike_slope >= spike_threshold, 2, trend_level)
			)

		df_agg <- rbind(df_agg, df_added)
	}
}

# Now roll up locations and calculate by county.
#
for (i in 1:length(COUNTIES)) {
  county <- COUNTIES[i]
  county_map <- df_counties %>% filter(location_counties_served == county)
  
  for (j in 1:length(TARGETS)) {
    targ <- TARGETS[j]
    
    df_this <- df_rs %>% 
      filter(target == targ & location_id %in% county_map$location_id) %>% 
      group_by(epi_date) %>% 
      reframe(mean_abundance = mean(target_copies_fn_per_cap, na.rm = TRUE),
              location_id = county, 
							target = first(target),
							cal_year = first(year), 
							cal_month = first(month), 
							cal_week = first(week), 
							epi_year = first(epi_year), 
							epi_week = first(epi_week))
    
    #print(paste0(county, ": ", targ, ". COUNT = ", nrow(df_this), sep=""))

        if (nrow(df_this) < 1) {
      next
    }
    
    df_added <- df_this %>% 
      arrange(epi_date) %>% 
      mutate(
        rolling_mean = rollmean(mean_abundance, k = ABUND_ALERT_WINDOW, fill = NA, align = "right")
      )
    
    df_added <- df_added %>% 
      mutate(
        abundance_fold_change = case_when(
          is.na(rolling_mean) ~ NaN,
          (mean_abundance <= ABUND_RESPONSE_THRESHOLD*rolling_mean) ~ NA,
          .default =  (mean_abundance - rolling_mean) / rolling_mean
        )
      )
    
    df_added <- df_added %>% 
      mutate(
        abundance_level = case_when(
          is.nan(abundance_fold_change) ~ 1,	# rolling mean is missing
          is.na(abundance_fold_change) ~ 2,		# rolling mean is too small to be useful
          (abundance_fold_change <= -0.5) ~ 3,
          (abundance_fold_change > -0.5 & abundance_fold_change <= 0.5) ~ 4,
          (abundance_fold_change > 0.5 & abundance_fold_change <= 1.5) ~ 5,
          (abundance_fold_change > 1.5) ~ 6
        )
      )
    
    trend_results <- rollapplyr(
      df_added,
      width = TREND_ALERT_WINDOW,
      FUN = get_LTcoefs,
      by.column = FALSE, # Apply the function to the entire window data frame
      fill = NA            # Fill initial empty values with NA
    )
    # Convert results to a data frame for easier use
    trend_coefs_df <- as.data.frame(trend_results)
    names(trend_coefs_df) <- c("trend_intercept", "trend_slope")
    # Add the results back to the original data frame
    df_added <- cbind(df_added, trend_slope = trend_coefs_df$trend_slope)
    
    trend_step <- df_thresholds %>% 
      filter(category == "trend" & target == targ) %>% 
      select(step)
    trend_step <- as.numeric(trend_step$step)
    
    df_added <- df_added %>%
      mutate(
        trend_level = case_when(
          is.na(trend_slope) ~ 1,
          (trend_slope <= -1.5*trend_step) ~ 3,
          (trend_slope > -1.5*trend_step & trend_slope <= -0.5*trend_step) ~ 4,
          (trend_slope > -0.5*trend_step & trend_slope <= 0.5*trend_step) ~ 5,
          (trend_slope > 0.5*trend_step & trend_slope <= 1.5*trend_step) ~ 6,
          (trend_slope > 1.5*trend_step) ~ 7
        )
      )
    
    # A drastic increase in mean abundance from one sample to the next is categorized as 
    # a "spike" and is handled differently in the alert system. To determine this, we 
    # calculate the trend line between each 2 consecutive samples.
    spike_results <- rollapplyr(
      df_added,
      width = 2,
      FUN = get_LTcoefs,
      by.column = FALSE, # Apply the function to the entire window data frame
      fill = NA            # Fill initial empty values with NA
    )
    # Convert results to a data frame for easier use
    spike_coefs_df <- as.data.frame(spike_results)
    names(spike_coefs_df) <- c("spike_intercept", "spike_slope")
    # Add the results back to the original data frame
    df_added <- cbind(df_added, spike_slope = spike_coefs_df$spike_slope)
    
    spike_threshold <- df_thresholds %>% 
      filter(category == "spike" & target == targ) %>% 
      select(step)
    spike_threshold <- as.numeric(spike_threshold$step)
    
    df_added <- df_added %>%
      mutate(
        trend_level = if_else(spike_slope >= spike_threshold, 2, trend_level)
      )
    
    df_agg <- rbind(df_agg, df_added)
  }
}

# Finally, do the entire state.
#
for (j in 1:length(TARGETS)) {
  targ <- TARGETS[j]
  
  df_this <- df_rs %>% 
    filter(target == targ) %>% 
    group_by(epi_date) %>% 
    reframe(mean_abundance = mean(target_copies_fn_per_cap, na.rm = TRUE),
            location_id = "WV", 
						target = first(target),
						cal_year = first(year), 
						cal_month = first(month), 
						cal_week = first(week), 
						epi_year = first(epi_year), 
						epi_week = first(epi_week))
  
  #print(paste0("WV: ", targ, ". COUNT = ", nrow(df_this), sep=""))
  
  if (nrow(df_this) < 1) {
    next
  }
  
  df_added <- df_this %>% 
    arrange(epi_date) %>% 
    mutate(
      rolling_mean = rollmean(mean_abundance, k = ABUND_ALERT_WINDOW, fill = NA, align = "right")
    )
  
  df_added <- df_added %>% 
    mutate(
      abundance_fold_change = case_when(
        is.na(rolling_mean) ~ NaN,
        (mean_abundance <= ABUND_RESPONSE_THRESHOLD*rolling_mean) ~ NA,
        .default =  (mean_abundance - rolling_mean) / rolling_mean
      )
    )
  
  df_added <- df_added %>% 
    mutate(
      abundance_level = case_when(
        is.nan(abundance_fold_change) ~ 1,	# rolling mean is missing
        is.na(abundance_fold_change) ~ 2,		# rolling mean is too small to be useful
        (abundance_fold_change <= -0.5) ~ 3,
        (abundance_fold_change > -0.5 & abundance_fold_change <= 0.5) ~ 4,
        (abundance_fold_change > 0.5 & abundance_fold_change <= 1.5) ~ 5,
        (abundance_fold_change > 1.5) ~ 6
      )
    )
  
  trend_results <- rollapplyr(
    df_added,
    width = TREND_ALERT_WINDOW,
    FUN = get_LTcoefs,
    by.column = FALSE, # Apply the function to the entire window data frame
    fill = NA            # Fill initial empty values with NA
  )
  # Convert results to a data frame for easier use
  trend_coefs_df <- as.data.frame(trend_results)
  names(trend_coefs_df) <- c("trend_intercept", "trend_slope")
  # Add the results back to the original data frame
  df_added <- cbind(df_added, trend_slope = trend_coefs_df$trend_slope)
  
  trend_step <- df_thresholds %>% 
    filter(category == "trend" & target == targ) %>% 
    select(step)
  trend_step <- as.numeric(trend_step$step)
  
  df_added <- df_added %>%
    mutate(
      trend_level = case_when(
        is.na(trend_slope) ~ 1,
        (trend_slope <= -1.5*trend_step) ~ 3,
        (trend_slope > -1.5*trend_step & trend_slope <= -0.5*trend_step) ~ 4,
        (trend_slope > -0.5*trend_step & trend_slope <= 0.5*trend_step) ~ 5,
        (trend_slope > 0.5*trend_step & trend_slope <= 1.5*trend_step) ~ 6,
        (trend_slope > 1.5*trend_step) ~ 7
      )
    )
  
  # A drastic increase in mean abundance from one sample to the next is categorized as 
  # a "spike" and is handled differently in the alert system. To determine this, we 
  # calculate the trend line between each 2 consecutive samples.
  spike_results <- rollapplyr(
    df_added,
    width = 2,
    FUN = get_LTcoefs,
    by.column = FALSE, # Apply the function to the entire window data frame
    fill = NA            # Fill initial empty values with NA
  )
  # Convert results to a data frame for easier use
  spike_coefs_df <- as.data.frame(spike_results)
  names(spike_coefs_df) <- c("spike_intercept", "spike_slope")
  # Add the results back to the original data frame
  df_added <- cbind(df_added, spike_slope = spike_coefs_df$spike_slope)
  
  spike_threshold <- df_thresholds %>% 
    filter(category == "spike" & target == targ) %>% 
    select(step)
  spike_threshold <- as.numeric(spike_threshold$step)
  
  df_added <- df_added %>%
    mutate(
      trend_level = if_else(spike_slope >= spike_threshold, 2, trend_level)
    )
  
  df_agg <- rbind(df_agg, df_added)
}


# Write the sum table to file.
# 
write.table(df_agg, file = paste(WVD_BASE, "/wvdash.rsstable.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)


# p_afc <- ggplot(df_agg, aes()) + 
#   geom_col(aes(x=epi_date, y=abundance_fold_change), fill="purple") + 
#   facet_wrap(~target, nrow=2, scales = "free")
# ggsave(paste0("abundance_fold_change.png", sep=""), plot = p_afc, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p_al <- ggplot(df_agg, aes()) + 
#   geom_histogram(aes(x=abundance_level, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), fill="yellow", color="darkblue") + 
#   labs(y="count") + 
#   facet_wrap(~target, nrow=2, scales = "free_y")
# ggsave(paste0("abundance_level.png", sep=""), plot = p_al, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p_tl <- ggplot(df_agg, aes()) + 
#   geom_histogram(aes(x=trend_level, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), fill="orange", color="black") + 
#   labs(y="count") + 
#   facet_wrap(~target, nrow=2, scales = "free_y")
# ggsave(paste0("trend_level.png", sep=""), plot = p_tl, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# 
# p_ssall <- ggplot(df_agg, aes()) + 
#   geom_histogram(aes(x=spike_slope, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), fill="yellow", color="darkblue") + 
#   labs(y="count") + 
#   facet_wrap(~target, nrow=2, scales = "free")
# ggsave(paste0("spike_slope_all.png", sep=""), plot = p_ssall, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# 
# p1 <- ggplot(df_agg %>% filter(target == TARGETS[1]), aes()) + 
#   geom_histogram(aes(x=trend_slope, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=10) + 
#   labs(y="count", title=paste0(TARGETS[1], " at binwidth 10", sep="")) + 
#   scale_x_continuous(limits=c(-250, 250))
# ggsave(paste0("trend_slope_", TARGETS[1], ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p2 <- ggplot(df_agg %>% filter(target == TARGETS[2]), aes()) + 
#   geom_histogram(aes(x=trend_slope, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=2) + 
#   labs(y="count", title=paste0(TARGETS[2], " at binwidth 2", sep="")) + 
#   scale_x_continuous(limits=c(-30, 30))
# ggsave(paste0("trend_slope_", TARGETS[2], ".png", sep=""), plot = p2, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p3 <- ggplot(df_agg %>% filter(target == TARGETS[3]), aes()) + 
#   geom_histogram(aes(x=trend_slope, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=1) + 
#   labs(y="count", title=paste0(TARGETS[3], " at binwidth 1", sep="")) + 
#   scale_x_continuous(limits=c(-20, 20)) 
# ggsave(paste0("trend_slope_", TARGETS[3], ".png", sep=""), plot = p3, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# p4 <- ggplot(df_agg %>% filter(target == TARGETS[4]), aes()) + 
#   geom_histogram(aes(x=trend_slope, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=1) + 
#   labs(y="count", title=paste0(TARGETS[4], " at binwidth 1", sep="")) + 
#   scale_x_continuous(limits=c(-20, 20))
# ggsave(paste0("trend_slope_", TARGETS[4], ".png", sep=""), plot = p4, width = 3306, height = 1892, units = c("px"), path = "analyses/")
# 
# 
