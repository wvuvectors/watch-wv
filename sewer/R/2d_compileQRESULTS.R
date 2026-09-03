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


targ2mix_df <- as.data.frame(
	read.table(
		paste0("resources/target2reaction_mix.txt", sep=""), 
		quote="", 
		sep="\t", 
		header=TRUE, 
		check.names=FALSE)
)

input_df <- as.data.frame(
	read.table(
		paste0(UPDIR, "/watchdb.update_files.txt", sep=""), 
		quote="", 
		sep="\t", 
		header=TRUE, 
		check.names=FALSE)
)

# Holds the available files for this data type.
update_df <- input_df %>% filter(batch_type == "Result")


# Holds data for the individual assays in this batch. Most of this data was already 
# extracted from the qbatch file. The rest will be updated from the incoming Result file.
# This df will be printed at the end of the compile.
run_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.quant_run.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    check.names=FALSE)
)
run_df$assay_target_copies_per_ul_reaction <- as.double(run_df$assay_target_copies_per_ul_reaction)

# Holds data for the individual controls in this batch. Most of this data was already 
# extracted from the qbatch file. The rest will be updated from the incoming Result file.
# This df will be printed at the end of the compile.
control_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.quant_control.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    check.names=FALSE)
)
control_df$assay_target_copies_per_ul_reaction <- as.double(control_df$assay_target_copies_per_ul_reaction)

# Just a reminder of the columns in these two dataframes. The current compile 
# will add values to the last three cols.
#   assay_id = character(),
#   assay_type = character(),
#   sample_id = character(),
#   assay_well = character(),
#   assay_block_id = character(),
#   assay_batch_id = character(),
# 	assay_target_name = character(),
# 	assay_target_abbreviation = character(),
# 	assay_target_fluorophore = character(),
# 	assay_target_copies_per_ul_reaction = double(),
# 	assay_accepted_partitions = double(),
# 	assay_positive_partitions = double()


# Loop over the Result files in this update.
#for (i in nrow(update_df)) {
	i = 1
	
	# Read the Results csv file into a dataframe.
  this_fn <- update_df$file_name[i]
  this_fpath <- paste0(UPDIR, "/BATCH_FILES/", this_fn, sep="")
	result_df <- fread(this_fpath, quote="")
	
	# Standardize the well ids to double digits.
	result_df <- result_df %>% mutate(
	  PlateRow = substr(Well, 1, 1),
	  PlateCol = as.numeric(substr(Well, 2, NULL))
	)
	result_df <- result_df %>% mutate(
	  PlateCol = if_else(PlateCol < 10, paste0(0, as.character(PlateCol), sep=""), as.character(PlateCol))
	)
	result_df <- result_df %>% unite("Well", PlateRow, PlateCol, sep = "", remove = FALSE)
	
	# Cleanup the keys so they are consistent.
	batch_keys <- colnames(result_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- str_replace_all(batch_keys, "\\.", "")
	batch_keys <- str_replace_all(batch_keys, "\\)", "")
	batch_keys <- str_replace_all(batch_keys, "\\(", "")
	batch_keys <- str_replace_all(batch_keys, "\\[", "")
	batch_keys <- str_replace_all(batch_keys, "\\]", "")
	batch_keys <- str_replace_all(batch_keys, "µ", "u")
	batch_keys <- str_replace_all(batch_keys, "/", "_per_")
	batch_keys <- tolower(batch_keys)
	colnames(result_df) <- batch_keys
  
	# Simplify the result_df to get rid of cols we don't need. First need to know: 
	# Is this a Qiacuity or QX result?
	if (result_df %>% colnames() %>% head(1) == "plate_name") {
    # This is a Qiacuity plate. Match to entries in run_df by assay_well, assay_batch_id, 
	  # and assay_target_abbreviation. Qiacuity does not output the fluorophore.
  	result_df <- result_df %>% select(
    	assay_batch_id = plate_name,
    	assay_target_abbreviation = target_name,
    	assay_target_copies_per_ul_reaction = conc_cp_per_ul_dpcr_reaction,
    	assay_accepted_partitions = partitions_valid,
	    assay_positive_partitions = partitions_positive,
	    assay_well = well
  	)
  	# Fill assay_batch_id column with value from first row.
  	result_df$assay_batch_id <- result_df$assay_batch_id[1]
  	
  	# Merge result_df into run_df and control_df.
  	merge_run_df <- rows_update(run_df, result_df, unmatched = "ignore", by = c("assay_well", "assay_batch_id", "assay_target_abbreviation"))
  	merge_control_df <- rows_update(control_df, result_df, unmatched = "ignore", by = c("assay_well", "assay_batch_id", "assay_target_abbreviation"))
  	
	} else if (result_df %>% colnames() %>% head(1) == "well") {
	  # This is a QX plate. Match the entries in run_df by assay_well, assay_batch_id, 
	  # and assay_target_fluorophore.
	  result_df <- result_df %>% select(
	    assay_target_fluorophore = dyenames,
	    assay_target_abbreviation = target,
	    assay_target_copies_per_ul_reaction = conccopies_per_uL,
	    assay_accepted_partitions = accepted_droplets,
	    assay_positive_partitions = positives,
	    assay_well = well
	  )
	  # Add assay_batch_id column using value from the file name. The batch id is 
	  # everything before the substring _analysis_.
	  fn_split_vec <- str_split_1(this_fn, pattern = fixed("_analysis_"))
	  result_df$assay_batch_id <- fn_split_vec[1]
	  
	  # Merge result_df into run_df and control_df.
	  merge_run_df <- rows_update(run_df, result_df, unmatched = "ignore", by = c("assay_well", "assay_batch_id", "assay_target_fluorophore"))
	  merge_control_df <- rows_update(control_df, result_df, unmatched = "ignore", by = c("assay_well", "assay_batch_id", "assay_target_fluorophore"))
	  
	} else {
	  # Unrecognizable result file. Set the nrows of result_df to 0.
	}
	
#}


if (nrow(merge_run_df) > 0) {

	# Simple validation. This checks if any result data is missing.
	merge_run_df <- merge_run_df %>% mutate(
		validation_key = case_when(
			is.na(assay_target_copies_per_ul_reaction) | assay_target_copies_per_ul_reaction == as.double(-1) ~ -1, 
			.default = 1
		)
	)

	# Print the merged run df only if there is result data.
	dfn <- paste0(UPDIR, "/update.quant_run_merged.txt", sep="")
	write.table(merge_run_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}


if (nrow(merge_control_df) > 0) {
	  
  # Simple validation. This checks if any result data is missing.
  merge_control_df <- merge_control_df %>% mutate(
    validation_key = case_when(
      is.na(assay_target_copies_per_ul_reaction) | assay_target_copies_per_ul_reaction == as.double(-1) ~ -1, 
      .default = 1
    )
  )
  
  # Print the merged control df only if there is result data.
  cfn <- paste0(UPDIR, "/update.quant_control_merged.txt", sep="")
	write.table(merge_control_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

# Return the number of rows in the run df so the controller script can do the right thing.
ret_val <- nrow(merge_run_df) + nrow(merge_control_df)
if (ret_val > 0) {
  if (nrow(merge_run_df %>% filter(validation_key == -1)) > 0 | nrow(merge_control_df %>% filter(validation_key == -1)) > 0) {
    ret_val <- -1 * ret_val
  }
}
 
# cat(ret_val)

