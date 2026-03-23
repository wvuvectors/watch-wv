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

#f <- file("stdin")
#open(f)
#while(length(line <- readLines(f, n = 1)) > 0) {
#  fpaths <- strsplit(line, " ")[[1]]
#}
#close(f)

fpaths <- c("../sewer/tmp/")
UPDIR <- fpaths[1]

input_df <- as.data.frame(
	read.table(
		paste0(UPDIR, "/update.batch_files.txt", sep=""), 
		quote="", 
		sep="\t", 
		header=TRUE, 
		check.names=FALSE)
)

# Holds the input rows for this data type.
update_df <- input_df %>% filter(batch_type == "concentration")

run_df <- data.frame(
  concentration_id = character(),
  sample_id = character(),
  concentration_well = character(),
  concentration_batch_id = character(),
  concentration_comment = character(),
  concentration_spike = character()
)

batch_df <- data.frame(
  concentration_batch_id = character(),
  concentration_batch_type = character(),
  concentration_batch_record_version = character(),
  concentration_date = character(),
  concentration_run_by = character(),
  concentration_machine = character(),
  concentration_batch_comment = character(),
  concentration_input_ml = character(),
  concentration_output_ml = character(),
  concentration_method = character(),
  concentration_method_lot_id = character()
)

# Loop over the new concentration batches in this update.
for (i in nrow(update_df)) {
	# Read the batch Excel file into a named list, with each sheet as a df.
  this_fn <- update_df$file_name[i]
  this_fpath <- paste0(UPDIR, "/batches/", this_fn, sep="")
	platef_in <- excel2df(this_fpath)
	
	# Extract each data sheet to make life easier.
	metadata_df <- as.data.frame(t(platef_in$Metadata))

	plate_df <- platef_in$Plate_Map %>% column_to_rownames(var = "...1")
	plate_df <- plate_df %>% select(where(~!all(is.na(.))))
	
	spike_df <- platef_in$Spike_Map %>% column_to_rownames(var = "...1")
	comment_df <- platef_in$Comment_Map %>% column_to_rownames(var = "...1")
	
	# Process the batch metadata.
	# This sheet contains a two column table, with keys in column 1 and vals in column 2.
	batchup_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(batchup_df) <- as.character(batchup_df[1, ])
	batchup_df <- batchup_df[-c(1:1), ]
	batchup_df <- batchup_df %>% slice(1:1)
	batchup_df <- batchup_df %>% select(where(~!all(is.na(.))))
	
	# Cleanup the keys so they are consistent.
	batch_keys <- colnames(batchup_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- tolower(batch_keys)
	batch_keys <- str_replace_all(batch_keys, "^", "concentration_")
	batch_keys <- str_replace_all(batch_keys, "^concentration_concentration_", "concentration_")
	colnames(batchup_df) <- batch_keys
	
	# Date munging, of course.
	batchup_df$concentration_date <- convertToDateTime(as.numeric(batchup_df$concentration_date))
	
	# Append the batch metadata to the main batch df.
	batch_df <- rbind(batch_df, batchup_df)

	# Init a df to hold the run data for this batch. Run data is found across multiple sheets. 
	runup_df <- data.frame(
	  concentration_id = character(),
	  sample_id = character(),
	  concentration_well = character(),
	  concentration_batch_id = character(),
	  concentration_comment = character(),
	  concentration_spike = character()
	)
	
	# Extract the batch id for convenience.
	batch_id <- batchup_df$concentration_batch_id[1]
	
	for (i in 1:nrow(plate_df)) {
		well_r <- toupper(rownames(plate_df)[i])
		for (j in 1:ncol(plate_df)) {
			well_c <- colnames(plate_df)[j]
			if (as.numeric(well_c) < 10) {
				well_c <- paste0(0, well_c, sep="")
			}
			sample_id <- plate_df[i, j]
			
			# Need a unique string to create a concentration id.
			# This uses the sample id and the current epoch time (seconds since 1/1/1970).
			# It ensures that even if a sample was run multiple times on the same plate, each will 
			# receive a UID.
			#date_str <- str_replace_all(batchup_df$concentration_date[1], "-", "")
			epoch_time <- as.numeric(now())
		
			if (!is.na(sample_id) & sample_id != "") {
			  runup_df <- add_row(
  			  runup_df,
	  		  concentration_id = paste0(sample_id, ".", epoch_time, sep=""),
		  	  sample_id = sample_id,
			    concentration_well = paste0(well_r, well_c, sep=""),
  				concentration_batch_id = batch_id,
  				concentration_comment = comment_df[i, j],
  				concentration_spike = spike_df[i, j]
  			)
			}
		}
	}
	
	# Add the run update df to the main run df.
	run_df <- rbind(run_df, runup_df)
}

# Only validate and print if there is run data.
if (nrow(run_df) > 0) {

	# Validate the input data against required columns.
	# Add a column to each output df to hold a validation key.
	# 	1 = all data present and proeprly formed.
	# 	0 = all required data present and properly formed, but at least one optional data field is missing or malformed.
	# 	-1 = at least one required data value missing or malformed.
	# If any rows are missing required data, alert the control script.
	# That is done by making the output value negative.
		
	run_df <- run_df %>% mutate(
		validation_key = case_when(
			is.na(sample_id) | sample_id == "" ~ -1, 
			is.na(concentration_batch_id) | concentration_batch_id == "" ~ -1, 
			.default = 1
		)
	)
	
	batch_df <- batch_df %>% mutate(
		validation_key = case_when(
			is.na(concentration_batch_id) | concentration_batch_id == "" ~ -1, 
			is.na(concentration_input_ml) | concentration_input_ml == "" ~ -1, 
			is.na(concentration_output_ml) | concentration_output_ml == "" ~ -1, 
			is.na(concentration_method) | concentration_method == "" ~ -1, 
			.default = 1
		)
	)

	# Print the run and batch df only if there is run data.
	dfn <- paste0(UPDIR, "/update.concentration.txt", sep="")
	write.table(run_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	bfn <- paste0(UPDIR, "/update.cbatch.txt", sep="")
	write.table(batch_df, file = bfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

# Return the number of rows in the run df so the controller script can do the right thing.
ret_val <- nrow(run_df)
if (ret_val > 0) {
	if (nrow(run_df %>% filter(validation_key == -1)) > 0 | nrow(batch_df %>% filter(validation_key == -1)) > 0) {
		ret_val <- -1 * ret_val
	}
}

cat(ret_val)

