#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# We use a ratio of the most recent X samples compared to the most recent Y samples. 
# X is set by the ABUND_ALERT_WINDOW variable and Y by the ABUND_ALERT_WINDOW_BASIS variable.
#
ABUND_ALERT_WINDOW_BASIS <- 12		# Number of consecutive samples to use as the basis for determining abundance alert level.
ABUND_ALERT_WINDOW <- 4						# Number of consecutive samples to use as the numerator in determining abundance alert level.
ABUND_RESPONSE_THRESHOLD <- 0.05	# Mean abundance values below this threshold (fraction of the rolling mean) are not used to calculate alerts.

# We use the slope over Z samples to determine the trend. Z is set by the TREND_ALERT_WINDOW variable.
#
TREND_ALERT_WINDOW <- 4						# Number of consecutive samples to use in trend calculations.


# Accepts 2 arguments on STDIN:
#		OUTDIR is the directory that will contain the output of this script.
#		DTSTAMP is a unique id to tag the archive file name (usually a datetime stamp).
#
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) > 0) {
#	OUTDIR <- args[1]
# DTSTAMP <- args[2]
#}
OUTDIR <- "../data/xtend/watch"
WFILE <- file.path(OUTDIR, "watch_join.LATEST.txt")

if (!exists("DTSTAMP")) {
	DTSTAMP <- format(Sys.time(), "%Y%m%d%H%M")
}

# Our output file needs to have the following columns:
#
# epi_date
# mean_abundance
# location_id
# location_category
# target
# target_disease
# watch_dashboard_category
# cal_year
# cal_month
# cal_week
# epi_year
# epi_week
# rolling_mean
# rolling_mean_basis
# abundance_fold_change
# abundance_level
# trend_slope
# trend_level
# spike_slope
# is_normalized


# We need info from the resources Excel file.
#
resource_tables <- c("target")
platef_in <- excel2df(file.path("resources", "watchdb.all_tables.xlsx"))
resource_list <- list()
for (tname in resource_tables) {
	r_df <- as.data.frame(platef_in[[tname]])
	colnames(r_df) <- as.character(r_df[1, ])
	r_df <- r_df[-1, ]
	r_df <- r_df %>% mutate(across(everything(), as.character))
	resource_list[[tname]] <- r_df
}
target_df <- resource_list[["target"]]


# Load in the watch table.
#
watch_df <- as.data.frame(
	read.table(
		WFILE, 
		quote="", 
		sep="\t", 
		header=TRUE, 
		fill=TRUE,
		comment.char="", 
		check.names=FALSE)
)
watch_df <- watch_df %>% mutate(across(everything(), as.character))
watch_df <- watch_df %>% mutate(across(starts_with("epi_"), as.numeric))

watch_df <- watch_df %>% 
	mutate(across(c(
		target_copies_per_ml_per_person_per_day, 
		target_copies_per_ml, 
		year, 
		month, 
		week), 
	as.numeric)
)

watch_df <- watch_df %>% filter(!is.na(epi_date)) %>% filter(!is.na(target_copies_per_ml))

watch_df$target_copies_per_L <- 1000 * watch_df$target_copies_per_ml
watch_df$target_copies_per_L_per_person_per_day <- 1000 * watch_df$target_copies_per_ml_per_person_per_day

# Convert date strings into Date objects.
#
watch_df$collection_start_datetime <- mdy_hm(watch_df$collection_start_datetime)
watch_df$collection_end_datetime <- mdy_hm(watch_df$collection_end_datetime)
watch_df$date_of_collection <- mdy(watch_df$date_of_collection)

# Add the dashboard category and filter for a valid dashboard assignment.
#
watch_df <- watch_df %>% left_join(
  target_df %>% select(target = target_id, watch_dashboard_category),
  by = "target"
)
watch_df <- watch_df %>% filter(watch_dashboard_category != "None")

# For convenience, we pull out the locations and targets into unnamed lists.
#
LOCATIONS <- unname(unique(watch_df$location_id))
TARGETS <- unname(unique(watch_df$target))

# Skip any extreme values.
watch_df <- watch_df %>% filter(is.extreme != TRUE)
# Skip any outlier values.
watch_df <- watch_df %>% filter(is.outlier != TRUE)

# We're going to parse out sub-tables by location and target to calculate alerts.
# Then we'll rowbind them all back together into this dataframe, agg_df.
agg_df <- data.frame()


# First rollup and calculate trend and level by location (facility).
#
for (i in 1:length(LOCATIONS)) {
	loc <- LOCATIONS[i]
	for (j in 1:length(TARGETS)) {
		targ <- TARGETS[j]
		#loc <- "StarCityWWTP-01"
		#targ <- "SARS-CoV-2 (N1)"
		this_df <- watch_df %>% filter(target == targ & location_id == loc)

		if (nrow(this_df) < 2) {
			#print(paste0("Skipping ", loc, ": ", targ, ". Only ", nrow(this_df), " rows.", sep=""))
		  next
		}

		no_flow_fraction <- as.numeric(this_df %>% summarise(na_count = sum(is.na(sample_flow)))) / nrow(this_df)

    if (no_flow_fraction <= 0.2) {
			this_df <- this_df %>% 
				group_by(epi_date) %>% 
				reframe(mean_abundance = mean(target_copies_per_L_per_person_per_day, na.rm = TRUE),
									location_id = first(location_id), 
									target = first(target),
									cal_year = first(year), 
									cal_month = first(month), 
									cal_week = first(week), 
									epi_year = first(epi_year), 
									epi_week = first(epi_week),
									watch_dashboard_category = first(watch_dashboard_category),
									is_normalized = TRUE)
    } else {
      this_df <- this_df %>% 
        group_by(epi_date) %>% 
        reframe(mean_abundance = mean(target_copies_per_L, na.rm = TRUE),
                location_id = first(location_id), 
                target = first(target),
                cal_year = first(year), 
                cal_month = first(month), 
                cal_week = first(week), 
                epi_year = first(epi_year), 
                epi_week = first(epi_week),
                watch_dashboard_category = first(watch_dashboard_category),
                is_normalized = FALSE)
    }
				
		if (nrow(this_df) < 5) {
			#print(paste0("Skipping ", loc, ": ", targ, ". Only ", nrow(this_df), " rows.", sep=""))
		  next
		}

		added_df <- this_df %>% 
			arrange(epi_date) %>% 
			mutate(
				rolling_mean = rollmean(mean_abundance, k = ABUND_ALERT_WINDOW, fill = NA, align = "right"),
				rolling_mean_basis = rollmean(mean_abundance, k = ABUND_ALERT_WINDOW_BASIS, fill = NA, align = "right")
			)

		added_df <- added_df %>% 
			mutate(
				abundance_fold_change = case_when(
					is.na(rolling_mean) | is.na(rolling_mean_basis) ~ NaN,
					(mean_abundance <= ABUND_RESPONSE_THRESHOLD*rolling_mean) ~ NA,
					.default =  (rolling_mean - rolling_mean_basis) / rolling_mean_basis
				)
			)

		added_df <- added_df %>% 
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
			added_df,
			width = TREND_ALERT_WINDOW,
			FUN = get_LTcoefs,
			by.column = FALSE, # Apply the function to the entire window data frame
			fill = NA            # Fill initial empty values with NA
		)
		# Convert results to a data frame for easier use
		trend_coefs_df <- as.data.frame(trend_results)
		names(trend_coefs_df) <- c("trend_intercept", "trend_slope")
		# Add the results back to the original data frame
		added_df <- cbind(added_df, trend_slope = trend_coefs_df$trend_slope)
    
		trend_step <- target_df %>% 
		  filter(target_id == targ) %>% 
		  select(watch_stepsize_trend)
		trend_step <- as.numeric(trend_step$watch_stepsize_trend)
		
		added_df <- added_df %>%
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
			added_df,
			width = 2,
			FUN = get_LTcoefs,
			by.column = FALSE, # Apply the function to the entire window data frame
			fill = NA            # Fill initial empty values with NA
		)
		# Convert results to a data frame for easier use
		spike_coefs_df <- as.data.frame(spike_results)
		names(spike_coefs_df) <- c("spike_intercept", "spike_slope")
		# Add the results back to the original data frame
		added_df <- cbind(added_df, spike_slope = spike_coefs_df$spike_slope)

		spike_threshold <- target_df %>% 
		  filter(target_id == targ) %>% 
		  select(watch_stepsize_spike)
		spike_threshold <- as.numeric(spike_threshold$watch_stepsize_spike)
		
		added_df <- added_df %>%
			mutate(
				trend_level = if_else(spike_slope >= spike_threshold, 2, trend_level)
			)

		agg_df <- rbind(agg_df, added_df)
	}
}

cfn <- file.path(OUTDIR, "watch_rollbyfac.LATEST.txt")
write.table(agg_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

afn <- file.path(OUTDIR, "archives", paste0("watch_rollbyfac.", DTSTAMP, ".txt", sep=""))
write.table(agg_df, file = afn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)


# p1 <- ggplot(df_all %>% filter(target == TARGETS[1]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-5, 5)) + 
#   labs(y="count", title=paste0(TARGETS[1], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[1], ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = "analyses/")



ret_val <- nrow(agg_df)
cat(ret_val)

