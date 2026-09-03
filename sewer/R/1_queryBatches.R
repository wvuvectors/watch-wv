#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")


# Accepts 2 directory paths on STDIN:
#   1. INDIR contains the lab data files to query for new batch data.
#		2. OUTDIR is a path to the directory where the data for this run will be written. 
#		3. ARCHFILE is a file that contains the ids of all batches processed already. This 
#			 info is used to avoid redundant processing of data.
# Arguments must be passed in the above order to this script.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  INDIR <- args[1]
	OUTDIR <- args[2]
	ARCHFILE <- args[3]
}

#INDIR <- "test_example/RAW_DATA"
#OUTDIR <- "test_example/PROCESSED_DATA/UPDATES/TEST"
#ARCHFILE <- "test_example/PROCESSED_DATA/UPDATES/completed_batches.txt"

# A dataframe to hold the batch files that need to be processed in this run.
bf_to_process_df <- data.frame(
  batch_type = character(),
  batch_id = character(),
  path_to_original = character(),
  file_name = character()
)


# List the batch files in the current input directory.
bf_to_check_list <- list.files(
	path = paste0(INDIR,  sep=""),
	pattern = "\\.xlsx|m$",
	full.names = TRUE
)
  
# Retrieve the batches that have already been processed. We don't want to run these again.
bf_done_df <- as.data.frame(
	read.table(
		ARCHFILE, 
		quote="", 
		sep="\t", 
		header=TRUE, 
		check.names=FALSE)
)

# Loop over the batch files in the input dir and flag batch files that need to
# be processed. This is determined by the Batch ID entry in the Metadata sheet 
# of the file.
keycol <- "batch_id"
for (fpath in bf_to_check_list) {

	platef_in <- excel2df(fpath)
	metadata_df <- as.data.frame(t(platef_in$Metadata))
	metadata_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(metadata_df) <- as.character(metadata_df[1, ])
	metadata_df <- metadata_df[-c(1:1), ]

	# Test this batch id against the existing batch IDs in batch_df. If it doesn't 
	# exist already, add it to the list of files to be processed.
	this_type <- metadata_df$`Batch Type`[1]
	this_id <- metadata_df$`Batch ID`[1]
	test_df <- bf_done_df %>% filter(.data[[keycol]] == this_id)
	if (nrow(test_df) == 0) {
		add_df <- data.frame(
			batch_type = this_type,
			batch_id = this_id,
			path_to_original = fpath,
			file_name = basename(fpath)
		)
		bf_to_process_df <- rbind(bf_to_process_df, add_df)
	}
}


# If any quantification batches are found, need to add the corresponding result file.
# If there is no result file, remove the quantification batch from this run.
qbatches_df <- bf_to_process_df %>% filter(batch_type == "Quantification")

for (i in 1:nrow(qbatches_df)) {

	bid <- qbatches_df$batch_id[i]
  result_files_list <- list.files(
    path = paste0(INDIR, sep=""),
    pattern = paste0("^", bid, ".+\\.csv$", sep=""),
    full.names = TRUE
  )

  if (length(result_files_list) == 0) {
    bf_to_process_df <- bf_to_process_df %>% filter(!(batch_id == bid))
  } else {
    for (fpath in result_files_list) {
  		addr_df <- data.frame(
  			batch_type = "Result",
  			batch_id = bid,
  			path_to_original = fpath,
			  file_name = basename(fpath)
  		)
  		bf_to_process_df <- rbind(bf_to_process_df, addr_df)
    }
  }
}

if (nrow(bf_to_process_df) > 0) {
  outfn <- paste0(OUTDIR, "/update_files.txt", sep="")
  write.table(bf_to_process_df, file = outfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	
	for (i in 1:nrow(bf_to_process_df)) {
		src_file <- bf_to_process_df$path_to_original[i]
		fname <- bf_to_process_df$file_name[i]
		file.copy(from = src_file, to = paste0(OUTDIR, "/BATCH_FILES/", fname, sep=""))
	}
}

# Need to send the result count back to bash for further processing.
ret_val <- as.numeric(nrow(bf_to_process_df))
cat(ret_val)



