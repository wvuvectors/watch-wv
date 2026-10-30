#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 argument on STDIN:
#   UPDIR is the directory that contains the update files for run being processed. UPDIR must 
# 	exist, and it must contain a folder "run_files" that contains "fixed" copies of the 
#		original run files.
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	UPDIR <- args[1]
}
#UPDIR <- "../data/updates/2026-10-22"
RUNID <- basename(UPDIR)


# Holds data for the individual assays in this batch. Most of this data was already 
# extracted from the qbatch file. The rest will be updated from the incoming Result file.
# This df will be printed at the end of the compile.
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
assay_df$assay_target_copies_per_ul_reaction <- as.character(assay_df$assay_target_copies_per_ul_reaction)

# Holds data for the individual controls in this batch. Most of this data was already 
# extracted from the qbatch file. The rest will be updated from the incoming Result file.
# This df will be printed at the end of the compile.
control_df <- as.data.frame(
  read.table(
    paste0(UPDIR, "/update.qcontrol.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
		fill=TRUE,
		comment.char="",
    check.names=FALSE)
)
control_df$assay_target_copies_per_ul_reaction <- as.character(control_df$assay_target_copies_per_ul_reaction)

# Just a reminder of the columns in these two dataframes. The current compile 
# will add values to the last three cols.
#   assay_id = character(),
#   assay_type = character(),
#   sample_id = character(),
#   assay_well = character(),
#   assay_block_id = character(),
#   quantification_batch_id = character(),
# 	assay_target_name = character(),
# 	assay_target_abbreviation = character(),
# 	assay_target_fluorophore = character(),
# 	assay_target_copies_per_ul_reaction = double(),
# 	assay_accepted_partitions = double(),
# 	assay_positive_partitions = double()



runfile_list <- list.files(
	path = paste0(UPDIR, "/run_files", sep=""),
	pattern = "\\.csv$",
	full.names = TRUE
)

# Loop over the csv files in this update.
for (fpath in runfile_list) {
	#fpath <- runfile_list[1]
	# Read the csv file into a dataframe.
	result_df <- fread(fpath, quote="")
	
	if ("Asset Tag ID" %in% names(result_df)) {
		# This is a Samples csv and is ignored.
		next
	}
	
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
	if ("plate_name" %in% names(result_df)) {
    # This is a Qiacuity plate. Match to entries in assay_df by assay_well, assay_batch_id, 
	  # and assay_target_abbreviation. Qiacuity does not output the fluorophore.
  	result_df <- result_df %>% select(
    	quantification_batch_id = plate_name,
    	assay_target_abbreviation = target_name,
    	assay_target_copies_per_ul_reaction = conc_cp_per_ul_dpcr_reaction,
    	assay_accepted_partitions = partitions_valid,
	    assay_positive_partitions = partitions_positive,
	    assay_well = well
  	)
  	# Fill assay_batch_id column with value from first row.
  	result_df$quantification_batch_id <- result_df$quantification_batch_id[1]
  	
  	# Merge result_df into assay_df and control_df.
  	merge_assay_df <- rows_update(assay_df, result_df, unmatched = "ignore", by = c("assay_well", "quantification_batch_id", "assay_target_abbreviation"))
  	merge_control_df <- rows_update(control_df, result_df, unmatched = "ignore", by = c("assay_well", "quantification_batch_id", "assay_target_abbreviation"))
  	
	} else {
	  # This is a QX plate. Match the entries in assay_df by assay_well, assay_batch_id, 
	  # and assay_target_fluorophore.
	  result_df <- result_df %>% select(
 			sample_description_1,
			sample_description_2,
			sample_description_3,
			sample_description_4,
	    assay_target_fluorophore = dyenames,
	    assay_target_abbreviation = target,
	    assay_target_copies_per_ul_reaction = conccopies_per_ul,
	    assay_accepted_partitions = accepted_droplets,
	    assay_positive_partitions = positives,
	    assay_well = well
	  )
	  
		result_df <- result_df %>% mutate(
			quantification_batch_id = case_when(
				all(!is.na(precheck_df$sample_description_1)) ~ precheck_df$sample_description_1[1],
				all(is.na(precheck_df$sample_description_1)) ~ precheck_df$sample_description_2[1],
				all(is.na(precheck_df$sample_description_2)) ~ precheck_df$sample_description_3[1],
				all(is.na(precheck_df$sample_description_3)) ~ precheck_df$sample_description_4[1],
				.default = "UNK"
			)
		)
		result_df <- result_df %>% select(!starts_with("sample_description_"))
		
	  # Merge result_df into assay_df and control_df.
	  merge_assay_df <- rows_update(assay_df, result_df, unmatched = "ignore", by = c("assay_well", "quantification_batch_id", "assay_target_fluorophore"))
	  merge_control_df <- rows_update(control_df, result_df, unmatched = "ignore", by = c("assay_well", "quantification_batch_id", "assay_target_fluorophore"))
	  
	}
}


if (nrow(merge_assay_df) > 0) {

	# Simple validation. This checks if any result data is missing.
	merge_assay_df <- merge_assay_df %>% mutate(
		validation_key = case_when(
			is.na(assay_target_copies_per_ul_reaction) | as.numeric(assay_target_copies_per_ul_reaction) == as.numeric(-1) ~ -1, 
			.default = 1
		)
	)

	# Print the merged assay df only if there is result data.
	dfn <- paste0(UPDIR, "/update.qassay_merged.txt", sep="")
	write.table(merge_assay_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}


if (nrow(merge_control_df) > 0) {
	  
  # Simple validation. This checks if any result data is missing.
  merge_control_df <- merge_control_df %>% mutate(
    validation_key = case_when(
      is.na(assay_target_copies_per_ul_reaction) | as.numeric(assay_target_copies_per_ul_reaction) == as.numeric(-1) ~ -1, 
      .default = 1
    )
  )
  
  # Print the merged control df only if there is result data.
  cfn <- paste0(UPDIR, "/update.qcontrol_merged.txt", sep="")
	write.table(merge_control_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

# Return the number of rows in the run df so the controller script can do the right thing.
ret_val <- nrow(merge_assay_df) + nrow(merge_control_df)
if (ret_val > 0) {
  if (nrow(merge_assay_df %>% filter(validation_key == -1)) > 0 | nrow(merge_control_df %>% filter(validation_key == -1)) > 0) {
    ret_val <- -1 * ret_val
  }
}
cat(ret_val)

