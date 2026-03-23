#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 directory path on STDIN:
#   UPDIR is the directory that contains the update files for new batch data. UPDIR must 
# 	exist, and it must contain several items that are the outcome of script 1:
# 	1. A file "update.batch_files.txt", a tab-delim table with batch type, id, file name, 
# 		 and absolute path to the original data file.
# 	2. A folder "batches" that contains copies of the original batch files.
#

# f <- file("stdin")
# open(f)
# while(length(line <- readLines(f, n = 1)) > 0) {
#  fpaths <- strsplit(line, " ")[[1]]
# }
# close(f)

fpaths <- c("../sewer/tmp/")
UPDIR <- fpaths[1]

# Holds the input rows for this data type.
input_df <- as.data.frame(
	read.csv(
		paste0(UPDIR, "/batches/Samples.csv", sep=""), 
		quote="", 
		header=TRUE, 
		check.names=FALSE)
)

input_df <- input_df %>% mutate(across(where(is.character), ~ na_if(.,"")))

cleaned_df <- input_df %>% filter(
  !is.na(`Sample Collection Start`) & 
  !is.na(`Sample Collection End`)
)

allsample_df <- cleaned_df %>% select(
	sample_id = `Asset Tag ID`, 
	location_id = `Location`, 
	sample_collection_method = `Sample Collection Method`, 
	sample_collection_by = `Sample Collection By`, 
	sample_event = `Event Type`,
	sample_qc = `Sample QC Check`,
	sample_collection_start_datetime = `Sample Collection Start`,
	sample_collection_end_datetime = `Sample Collection End`,
	sample_recovered_datetime = `Sample Retrieved Date/Time`,
	sample_flow = `Sample Flow (MGD)`,
	sample_received_by = `Sample Received By`,
	sample_received_date = `Sample Received Date`,
	sample_comment = `Comments`)

#allsample_df$sample_event[is.na(allsample_df$sample_event)] <- "Routine Surveillance"

# Get the sample IDs present in the current update.
newsample_ids <- c()
id_files <- c("update.concentration.txt", "update.extraction.txt", "update.assay.txt")
for (fn in id_files) {
	fpath <- paste0(UPDIR, "/", fn, sep="")
	if (file.exists(fpath)) {
		id_df <- as.data.frame(
			read.table(fpath, sep = "\t", quote="", header=TRUE, check.names=FALSE)
		)
		newsample_ids <- c(newsample_ids, id_df$sample_id)
	}
}
newsample_ids <- unique(newsample_ids)

# Get the sample data for the samples in the current update.
update_df <- allsample_df %>% filter(sample_id %in% newsample_ids)


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
			is.na(sample_event) | sample_event == "" ~ -1, 
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

cat(ret_val)

