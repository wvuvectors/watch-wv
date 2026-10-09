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
#UPDIR <- "test_example/updates/2026-42-42"
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
samples_vec <- unique(unlist(
	c(
		conc_df$sample_id, 
		extract_df$sample_id, 
		assay_df$sample_id
	), use.names = FALSE))



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
	county_id = character()
)

runfile_list <- list.files(
	path = paste0(UPDIR, "/run_files", sep=""),
	pattern = "\\.xlsx|m$",
	full.names = TRUE
)

for (fpath in runfile_list) {
	# Read the batch file into a named list, with each sheet as a df.
	platef_in <- excel2df(fpath)
	
	if (!("Export" %in% names(platef_in))) {
		next
	}
	
	# Read the Export sheet.
	allsample_df <- as.data.frame(platef_in[["Export"]])
	colnames(allsample_df) <- as.character(allsample_df[1, ])
	allsample_df <- allsample_df[-1, ]

	if (colnames(allsample_df)[1] != "Asset Tag ID") {
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
	
	update_df <- allsample_df %>% 
		filter(asset_tag_id %in% samples_vec) %>% 
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
		validation_key = case_when(
			is.na(sample_id) | sample_id == "" ~ -1, 
			is.na(location_id) | location_id == "" ~ -1, 
			#is.na(sample_event) | sample_event == "" ~ -1, 
			is.na(sample_collection_start_datetime) | sample_collection_start_datetime == "" ~ -1, 
			is.na(sample_collection_end_datetime) | sample_collection_end_datetime == "" ~ -1, 
			.default = 1
		)
	)
	
	# Print the update df only if there is sample data.
	rfn <- paste0(UPDIR, "/update.sample.txt", sep="")
	write.table(update_df, file = rfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

ret_val <- nrow(update_df)
if (ret_val > 0) {
	if (nrow(update_df %>% filter(validation_key == -1)) > 0) {
		ret_val <- -1 * ret_val
	}
}

cat(ret_val)

