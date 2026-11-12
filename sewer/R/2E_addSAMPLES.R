#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 argument on STDIN:
#   UPDIR is the directory that contains the update files for run being processed. UPDIR must 
# 	exist, and it must contain a folder "run_files" that contains "fixed" copies of the 
#		original run files, including the samples csv file.
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	UPDIR <- args[1]
}
#UPDIR <- "../data/updates/2026-10-15"
RUNID <- basename(UPDIR)


# Get all the sample IDs that are found across this update (conc, extract, or assay).
conc_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.concentration.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
		fill=TRUE,
		comment.char="",
    check.names=FALSE)
)
extract_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.extraction.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
		fill=TRUE,
		comment.char="",
    check.names=FALSE)
)
assay_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.qassay.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
		fill=TRUE,
		comment.char="",
    check.names=FALSE)
)

samples_in_update_df <- data.frame(
  sample_id = unique(unlist(list(conc_df$sample_id, extract_df$sample_id, assay_df$sample_id)))
)

update_df <- data.frame(
	sample_id = character(), 
	location_id = character(), 
	sample_collection_method= character(), 
	sample_collection_by = character(), 
	sample_event = character(),
	sample_qc = character(),
	sample_collection_start_datetime = character(),
	sample_collection_end_datetime = character(),
	sample_recovered_datetime = character(),
	sample_flow = character(),
	sample_received_by = character(),
	sample_received_date = character(),
	sample_comment = character(),
	location_type = character(), 
	county_id = character(),
	sample_validation_key = character()
)

runfile_list <- list.files(
	path = paste0(UPDIR, "/run_files", sep=""),
	pattern = "\\.xlsx|m$",
	full.names = TRUE
)

for (fpath in runfile_list) {
  #fpath <- "../data/updates/2026-10-15/run_files/collections.xlsx"

  # Read the batch file into a named list, with each sheet as a df.
	platef_in <- excel2df(fpath)
	
	# Do we need this AND the column name check below?
	if (length(names(platef_in)) > 1) {
		next
	}
	
	# Read the first sheet.
	allsample_df <- as.data.frame(platef_in[1])
	colnames(allsample_df) <- as.character(allsample_df[1, ])
	allsample_df <- allsample_df[-1, ]

	if (!"Asset Tag ID" %in% colnames(allsample_df)) {
		# This is not a Samples table and is ignored.
		next
	}
	
	#allsample_df <- allsample_df %>% mutate(across(where(is.character), ~ na_if(.,"")))
	
	# Cleanup the keys so they are consistent.
	batch_keys <- colnames(allsample_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- str_replace_all(batch_keys, "\\.", "")
	batch_keys <- str_replace_all(batch_keys, "\\)", "")
	batch_keys <- str_replace_all(batch_keys, "\\(", "")
	batch_keys <- str_replace_all(batch_keys, "\\[", "")
	batch_keys <- str_replace_all(batch_keys, "\\]", "")
	batch_keys <- str_replace_all(batch_keys, "µ", "u")
	batch_keys <- str_replace_all(batch_keys, "/", "")
	batch_keys <- tolower(batch_keys)
	colnames(allsample_df) <- batch_keys
	
	# This function is part of the openxlsx package. It should take an Excel DateTime (which is numeric)
	# and convert it to a POSIXct object in the format YYYY-MM-DD. It appears the AssetTiger export 
	# uses strings instead, though, so no need.
	#
	# allsample_df <- allsample_df %>% mutate(
	# 	sample_collection_start_datetime = if_else(
	# 		is.na(sample_collection_start) | sample_collection_start == "", 
	# 		NA, 
	# 		convertToDateTime(as.numeric(sample_collection_start)),
	# 		missing = NA
	# 	),
	# 	sample_collection_end_datetime = if_else(
	# 		is.na(sample_collection_end) | sample_collection_end == "", 
	# 		NA, 
	# 		convertToDateTime(as.numeric(sample_collection_end)),
	# 		missing = NA
	# 	),
	# 	sample_recovered_datetime = if_else(
	# 		is.na(sample_retrieved_datetime) | sample_retrieved_datetime == "", 
	# 		NA, 
	# 		convertToDateTime(as.numeric(sample_retrieved_datetime)),
	# 		missing = NA
	# 	)
	# )

	# Restrict allsamples to just those assets in the update files.
	#
	this_df <- allsample_df %>% 
		filter(asset_tag_id %in% samples_in_update_df$sample_id) %>% 
		select(
			sample_id = asset_tag_id, 
			location_id = location, 
			sample_collection_method, 
			sample_collection_by, 
			sample_event = event_type,
			sample_qc = sample_qc_check,
			sample_collection_start_datetime = sample_collection_start,
			sample_collection_end_datetime = sample_collection_end,
			sample_recovered_datetime = sample_retrieved_datetime,
			sample_flow = sample_flow_mgd,
			sample_received_by,
			sample_received_date,
			sample_comment = comments,
			location_type = department, 
			county_id = site
		)
		
	# this_df <- this_df %>% mutate(
	# 	sample_collection_start_datetime = format(sample_collection_start_datetime, "%m/%d/%Y %H:%M"),
	# 	sample_collection_end_datetime = format(sample_collection_end_datetime, "%m/%d/%Y %H:%M"),
	# 	sample_recovered_datetime = format(sample_recovered_datetime, "%m/%d/%Y %H:%M")
	# )
	
	update_df <- rbind(update_df, this_df)
}

#
#
# If sample id in the assay, conc, or extraction tables is NOT in the AT database, remove 
# and note. May need to re-write these tables.
# This requires filtering samples_in_update_vec for ids NOT in allsample_df.
#
external_samples_df <- samples_in_update_df %>% filter(!sample_id %in% allsample_df$asset_tag_id)

if (nrow(external_samples_df) > 0) {
	conc_ext_df <- conc_df %>% filter(sample_id %in% external_samples_df$sample_id)
	conc_filt_df <- conc_df %>% filter(!sample_id %in% external_samples_df$sample_id)
	if (nrow(conc_ext_df) > 0 & nrow(conc_ext_df) < nrow(conc_df)) {
		write.table(conc_ext_df, file = paste0(UPDIR, "/external.concentration.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
		write.table(conc_filt_df, file = paste0(UPDIR, "/update.concentration.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	}

	extract_ext_df <- extract_df %>% filter(sample_id %in% external_samples_df$sample_id)
	extract_filt_df <- extract_df %>% filter(!sample_id %in% external_samples_df$sample_id)
	if (nrow(extract_filt_df) > 0 & nrow(extract_filt_df) < nrow(extract_df)) {
		write.table(extract_ext_df, file = paste0(UPDIR, "/external.extraction.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
		write.table(extract_filt_df, file = paste0(UPDIR, "/update.extraction.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	}

	assay_ext_df <- assay_df %>% filter(sample_id %in% external_samples_df$sample_id)
	assay_filt_df <- assay_df %>% filter(!sample_id %in% external_samples_df$sample_id)
	if (nrow(assay_filt_df) > 0 & nrow(assay_filt_df) < nrow(assay_df)) {
		write.table(assay_ext_df, file = paste0(UPDIR, "/external.qassay.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
		write.table(assay_filt_df, file = paste0(UPDIR, "/update.qassay.txt", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	}
}


# Only validate and print if there is update data.
if (nrow(update_df) > 0) {

	# Validate the input data against required columns.
	# Add a column to each output df to hold a validation key.
	# 	1 = all data present and properly formed.
	# 	0 = all required data present and properly formed, but at least one optional data field is missing or malformed.
	# 	-1 = at least one required data value missing or malformed.
	# If any rows are missing required data, alert the control script.
	# That is done by making the output value negative.
		
	update_df <- update_df %>% mutate(
		sample_validation_key = case_when(
			!is.na(sample_id) & sample_id != "" ~ "PASS", 
			!is.na(location_id) & location_id != "" ~ "PASS", 
			#!is.na(sample_event) & sample_event != "" ~ "PASS", 
			!is.na(sample_collection_start_datetime) & sample_collection_start_datetime != "" ~ "PASS", 
			!is.na(sample_collection_end_datetime) & sample_collection_end_datetime != "" ~ "PASS", 
			.default = "FAIL"
		)
	)
	
	# Print the update df.
	rfn <- paste0(UPDIR, "/update.sample.txt", sep="")
	write.table(update_df, file = rfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

ret_val <- nrow(update_df)
if (ret_val > 0) {
	if (nrow(update_df %>% filter(sample_validation_key == "FAIL")) > 0) {
		ret_val <- -1 * ret_val
	}
}

cat(ret_val)

