#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 4 argument on STDIN:
#   INDIR is the directory that contains the latest watch data tables.
#		OUTDIR is the directory that will contain the output of this script.
#		MUF is the full path to the Marshall U data file, a csv in NWSS format.
#		DTSTAMP is a unique id to tag the archive file name (usually a datetime stamp).
#
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) > 0) {
#	INDIR <- args[1]
#	OUTDIR <- args[2]
#	MUF <- args[3]
# DTSTAMP <- args[4]
#}
INDIR <- "../data/latest"
OUTDIR <- "../data/xtend/watch"
MUF <- "/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/mu_nwss.LATEST.csv"

if (!exists("DTSTAMP")) {
	DTSTAMP <- format(Sys.time(), "%Y%m%d%H%M")
}

# We need info from the resources Excel file.
#
resource_tables <- c("target", "location")
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
location_df <- resource_list[["location"]]



# Load in the sewer result table.
#
sewer_df <- as.data.frame(
	read.table(
		file.path(INDIR, "watchdb.result.txt"), 
		quote="", 
		sep="\t", 
		header=TRUE, 
		fill=TRUE,
		comment.char="", 
		check.names=FALSE)
)
sewer_df <- sewer_df %>% mutate(across(everything(), as.character))

# Convert date strings into Date objects.
#
sewer_df$collection_start_datetime <- mdy_hm(sewer_df$collection_start_datetime)
sewer_df$collection_end_datetime <- mdy_hm(sewer_df$collection_end_datetime)

# Restrict our results to just those targets that are found in our target dataframe.
# This gets rid of BRSV, human RP, etc.
# Also restrict to active locations only.
#
sewer_df <- sewer_df %>% 
  filter(target %in% target_df$target_id) %>% 
  filter(location_id %in% location_df$location_id) %>% 
  filter(tolower(sample_qc) == "pass")

# Trim to a minimum watch format.
sewer_df <- sewer_df %>% 
	select(!c(
		run_id, 
		sample_qc, 
		assay_block_handle_reps, 
		aggregate_count, 
		target_copies_per_day, 
		target_copies_per_person_per_day, 
		target_copies_per_ml_per_day
		)
	)

# Add the target_disease.
#
sewer_df <- sewer_df %>% 
	left_join(
		target_df %>% select(
			target = target_id, 
			target_disease
	), by = "target")



# Load in the MU data. This should be a csv file in NWSS format.
#
mu_nwss_df <- fread(MUF, quote="\"")
mu_nwss_df <- mu_nwss_df %>% mutate(across(everything(), as.character))
mu_nwss_df <- mu_nwss_df %>% filter(tolower(ntc_amplify) != "yes")

# Convert MU data from NWSS to minimal WaTCH format.
#
# result_id
# sample_id
# lab_id
# event_type
# collection_start_datetime
# collection_end_datetime
# sample_flow
# target_copies_per_ml
### target
### location_id
### location_category
## collection_window_hours
## target_copies_per_ml_per_person_per_day
## location_population_served
## location_counties_served
# run_id	OMIT
# sample_qc	OMIT
# assay_block_handle_reps	OMIT
# aggregate_count	OMIT
# target_copies_per_day	OMIT
# target_copies_per_person_per_day	OMIT
# target_copies_per_ml_per_day	OMIT

mu_df <- data.frame(
  sample_id = mu_nwss_df$sample_id,
	lab_id = mu_nwss_df$lab_id, 
	event_type = "Routine Surveillance",
	collection_start_datetime = mdy_hm(paste0(mu_nwss_df$sample_collect_date, " ", mu_nwss_df$sample_collect_time)), # 3/12/2022,0:01
	collection_end_datetime = mdy_hm(paste0(mu_nwss_df$sample_collect_date, " ", mu_nwss_df$sample_collect_time)), # 3/12/2022,0:01
	sample_flow = mu_nwss_df$flow_rate,
	target_copies_per_ml = (as.numeric(mu_nwss_df$pcr_target_avg_conc))/1000,
	nwss_target = tolower(mu_nwss_df$pcr_target),
	nwss_gene_target = tolower(mu_nwss_df$pcr_gene_target),
	nwss_site_id = mu_nwss_df$site_id
)


target_df <- target_df %>% mutate(
  nwss_target = tolower(target_df$nwss_target),
  nwss_gene_target = tolower(target_df$nwss_gene_target)	
)

# Add the target and target_disease, using the NWSS data as combined key for target_df. 
# Also add a suffix to the sample ID to make a unique result id.
#
mu_df <- mu_df %>% 
	left_join(
		target_df %>% select(
			nwss_gene_target, 
			nwss_target, 
			target=target_id, 
			target_disease
	), by = join_by(nwss_target, nwss_gene_target))

# Create a unique result ID for every row.
#
mu_df$TMP_indx <- seq_len(nrow(mu_df))

mu_df <- mu_df %>% mutate(
  result_id = paste0(sample_id, "-", TMP_indx, sep="")
) %>% select(-TMP_indx)
 
    
# Fix the location_id and category.
#
mu_df <- mu_df %>% 
	left_join(
		location_df %>% select(
			nwss_site_id, 
			location_id, 
			location_category,
			location_counties_served,
			location_population_served,
			collection_window_hours = location_collection_window_hrs
	), by = "nwss_site_id")

# mu_df <- mu_df %>% mutate(
#   collection_window_hours = if_else(
#     as.numeric(difftime(collection_end_datetime, collection_start_datetime, units="hours")) > 3, 
#     ceiling(as.numeric(difftime(collection_end_datetime, collection_start_datetime, units="hours"))),
#     -1
#   )
# )

# Calculate normalized abundance from target copies, pop served, and collection window.
# 
# Some samples with no flow recorded have a dash instead; convert these to NAs.
#
mu_df <- mu_df %>% mutate(sample_flow = na_if(sample_flow, "-"))
mu_df <- mu_df %>% mutate(
	target_copies_per_ml_per_person_per_day = if_else(
		!is.na(mu_df$sample_flow) & !is.na(mu_df$location_population_served),
		as.numeric(mu_df$target_copies_per_ml) * 
		as.numeric(mu_df$sample_flow) * 3785.41 * 10^6 * 
		as.numeric(mu_df$collection_window_hours) / 
		as.numeric(mu_df$location_population_served),
		NA
	)
)

mu_df <- mu_df %>% select(!c(nwss_gene_target, nwss_target, nwss_site_id))


# rbind SEWER and MU tables. This will be used to build the dashboard table.
#
watch_df <- rbind(sewer_df, mu_df)


#
# Dates, dates, dates.
#

# Create a date to use as the primary reference and make sure all tables use the same date.
watch_df$date_of_collection <- as.Date(watch_df$collection_end_datetime)

# Add epi week info.
watch_df$year <- lubridate::year(watch_df$date_of_collection)
watch_df$month <- lubridate::month(watch_df$date_of_collection)
watch_df$week <- lubridate::week(watch_df$date_of_collection)
watch_df$epi_week <- lubridate::epiweek(watch_df$date_of_collection)
watch_df <- watch_df %>%
  mutate(epi_year = case_when(
    watch_df$month == 12 & watch_df$epi_week == 1 ~ watch_df$year + 1,
    watch_df$month == 1 & watch_df$epi_week == 53 ~ watch_df$year - 1,
    TRUE ~ watch_df$year
  ))

watch_df <- watch_df %>% 
	mutate(epi_date = as.integer(epi_year) + as.integer(epi_week)/100)


# For convenience, we pull out the locations and targets into unnamed lists.
#
LOCATIONS <- unname(unique(watch_df$location_id))
TARGETS <- unname(unique(watch_df$target))

watch_df <- watch_df %>% mutate(
  target_copies_per_ml = case_when(
    is.nan(target_copies_per_ml) ~ NA,
    is.infinite(target_copies_per_ml) ~ NA,
    .default = as.numeric(target_copies_per_ml)
  ),
  target_copies_per_ml_per_person_per_day = case_when(
    is.nan(target_copies_per_ml_per_person_per_day) ~ NA,
    is.infinite(target_copies_per_ml_per_person_per_day) ~ NA,
    .default = as.numeric(target_copies_per_ml_per_person_per_day)
  ),
  sample_flow = case_when(
    is.nan(sample_flow) ~ NA,
    is.infinite(sample_flow) ~ NA,
    .default = as.numeric(sample_flow)
  )
)

  
# Loop over every target for every location and calculate a z-score.
# We use this to identify outliers and extreme values; the latter are excluded 
# from the calculation of trend and level.
#
agg_df <- data.frame()
for (i in 1:length(TARGETS)) {
  targ <- TARGETS[i]
  for (j in 1:length(LOCATIONS)) {
    loc <- LOCATIONS[j]
    #targ <- "SARS-CoV-2 (N2)"
    #loc <- "StarCityWWTP-01"
    this_df <- watch_df %>% filter(target == targ & location_id == loc)
    #		df_added <- this_df %>% identify_outliers(target_copies_fn_per_cap)
    #		df_added <- df_added %>% select(assay_id, is.outlier, is.extreme)
    
    if (nrow(this_df) < 1) {
      next
    }

    mean_norm <- mean(this_df$target_copies_per_ml_per_person_per_day, na.rm=TRUE)
    sd_norm <- sd(this_df$target_copies_per_ml_per_person_per_day, na.rm=TRUE)
    
    mean_raw <- mean(as.numeric(this_df$target_copies_per_ml), na.rm=TRUE)
    sd_raw <- sd(as.numeric(this_df$target_copies_per_ml), na.rm=TRUE)
    
    this_df <- this_df %>%
      mutate(z = if_else(
        tolower(event_type) == "routine surveillance" & !is.na(sample_flow),
        (as.numeric(target_copies_per_ml_per_person_per_day) - mean_norm) / sd_norm,
        (as.numeric(target_copies_per_ml) - mean_raw) / sd_raw)
      )
    
    zscore_df <- target_df %>% 
      filter(target_id == targ) %>% 
      select(watch_stepsize_z)
    z_step <- as.numeric(zscore_df$watch_stepsize_z)
    
    this_df$is.outlier <- (this_df$z > (z_step * 1))
    this_df$is.extreme <- (this_df$z > (z_step * 2))
    
    #low <- median(this_df$target_copies_fn_per_cap) - 3 * mad(this_df$target_copies_fn_per_cap, constant=1)
    #high <- median(this_df$target_copies_fn_per_cap) + 3 * mad(this_df$target_copies_fn_per_cap, constant=1)
    #subset dataframe where points value is outside of low and high bounds
    #this_df$hampel <- (this_df$target_copies_fn_per_cap<low | this_df$target_copies_fn_per_cap>high)
    
    added_df <- this_df %>% 
      select(result_id, z, is.outlier, is.extreme) %>% 
      drop_na(z, is.outlier, is.extreme)
    agg_df <- rbind(agg_df, added_df)
  }
}

scored_df <- distinct(watch_df %>% left_join(agg_df, by="result_id", relationship = "many-to-many"))
scored_df <- scored_df %>% mutate(
  collection_end_datetime = format(collection_end_datetime, "%m/%d/%Y %H:%M"),
  collection_start_datetime = format(collection_start_datetime, "%m/%d/%Y %H:%M"),
  date_of_collection = format(date_of_collection, "%m/%d/%Y")
)

# Write output table to files.
#
cfn <- file.path(OUTDIR, "/watch_join.LATEST.txt")
write.table(scored_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

afn <- file.path(OUTDIR, "archives", paste0("watch_join.", DTSTAMP, ".txt", sep=""))
write.table(scored_df, file = afn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)


# p1 <- ggplot(df_all %>% filter(target == TARGETS[1]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-5, 5)) + 
#   labs(y="count", title=paste0(TARGETS[1], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[1], ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = "analyses/")


ret_val <- nrow(scored_df)
cat(ret_val)

#nrow(scored_df %>% filter(collection_window_hours < 0))
#nrow(scored_df %>% filter(target_copies_per_ml_per_person_per_day < 0))

