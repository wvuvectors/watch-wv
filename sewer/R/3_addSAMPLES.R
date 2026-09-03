#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 directory path on STDIN:
#   UPDIR is the directory that contains the update files for new batch data. UPDIR must 
# 	exist, and it must contain several items that are the outcome of script 1:
# 	1. A file "watchdb.update_files.txt", a tab-delim table with batch type, id, file name, 
# 		 and absolute path to the original data file.
# 	2. A folder "BATCH_FILES" that contains copies of the original batch files.
#

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) > 0) {
#   UPDIR <- args[1]
# }
UPDIR <- "test_example/PROCESSED_DATA/UPDATES/TEST"


# Holds the input rows for this data type.
input_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/watchdb.update_files.txt", sep=""), 
    quote="", 
    header=TRUE, 
    check.names=FALSE)
)


# Holds the available files for this data type.
update_df <- input_df %>% filter(batch_type == "Sample")


# Get the sample IDs for this update. We only need to handle the samples with run 
# data (from the quant_run_merged file).
run_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.quant_run_merged.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    check.names=FALSE)
)


#for (i in nrow(update_df)) {
i = 1

# Read the Sample csv file into a dataframe.
this_fn <- update_df$file_name[i]
this_fpath <- paste0(UPDIR, "/BATCH_FILES/", this_fn, sep="")
allsample_df <- fread(this_fpath)


allsample_df <- allsample_df %>% mutate(across(where(is.character), ~ na_if(.,"")))

# Cleanup the keys so they are consistent.
batch_keys <- colnames(allsample_df)
batch_keys <- str_replace_all(batch_keys, " ", "_")
batch_keys <- str_replace_all(batch_keys, "\\.", "")
batch_keys <- str_replace_all(batch_keys, "\\)", "")
batch_keys <- str_replace_all(batch_keys, "\\(", "")
batch_keys <- str_replace_all(batch_keys, "\\[", "")
batch_keys <- str_replace_all(batch_keys, "\\]", "")
batch_keys <- str_replace_all(batch_keys, "µ", "u")
batch_keys <- str_replace_all(batch_keys, "/", "_per_")
batch_keys <- tolower(batch_keys)
colnames(allsample_df) <- batch_keys

sample_df <- allsample_df %>% select(
	sample_id = asset_tag_id, 
	location_id = location, 
	sample_collection_method, 
	sample_collection_by, 
	sample_event = event_type,
	sample_qc = sample_qc_check,
	sample_collection_start_datetime = sample_collection_start,
	sample_collection_end_datetime = sample_collection_end,
	sample_recovered_datetime = sample_retrieved_date_per_time,
	sample_flow = sample_flow_mgd,
	sample_received_by,
	sample_received_date,
	sample_comment = comments)


cleaned_df <- sample_df %>% filter(
  !is.na(sample_collection_start_datetime) & 
    !is.na(sample_collection_end_datetime) & 
    !is.na(location_id)
)

# Get the sample data for the samples in the current update.
update_df <- cleaned_df %>% filter(sample_id %in% run_df$sample_id)


# Only validate and print if there is update data.
if (nrow(update_df) > 0) {

	# Validate the input data against required columns.
	# Add a column to each output df to hold a validation key.
	# 	1 = all data present and proeprly formed.
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
	
	# Print the update df only if there is run data.
	rfn <- paste0(UPDIR, "/update.sample.txt", sep="")
	write.table(update_df, file = rfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

# Return the number of rows in the update df so the controller script can do the right thing.
ret_val <- nrow(update_df)
if (ret_val > 0) {
	if (nrow(update_df %>% filter(validation_key == -1)) > 0) {
		ret_val <- -1 * ret_val
	}
}

#cat(ret_val)

