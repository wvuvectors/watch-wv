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

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) > 0) {
#  UPDIR <- args[1]
#}
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
update_df <- input_df %>% filter(batch_type == "Quantification")

# Holds the batch metadata.
batch_df <- data.frame(
  batch_type = character(),
  quantification_batch_record_version = character(),
  quantification_batch_machine = character(),
  quantification_batch_analysis_software_version = character(),
  quantification_batch_amplification_method = character(),
  quantification_batch_amplification_method_lot_id = character(),
  quantification_batch_quantification_type = character(),
  quantification_batch_date = character(),
  quantification_batch_run_by = character(),
  quantification_batch_id = character(),
  quantification_batch_comment = character()
)


# Holds data for the assay blocks in this batch.
block_df <- data.frame(
  assay_block_id = character(),
  assay_block_batch_id = character(),
  assay_block_method = character(),
  assay_block_method_lot_id = character(),
  assay_block_sample_input_vol_ul = character(),
  assay_block_reaction_vol_ul = character(),
  assay_block_comment = character()
)

# Holds data for the individual assays in this batch. Some of this data is extracted from 
# the results file, the rest from the Excel "plate" file.
run_df <- data.frame(
  assay_id = character(),
  assay_type = character(),
  sample_id = character(),
  assay_well = character(),
  assay_block_id = character(),
  assay_batch_id = character(),
	assay_target_name = character(),
	assay_target_abbreviation = character(),
	assay_target_fluorophore = character(),
	assay_target_copies_per_ul_reaction = double(),
	assay_accepted_partitions = double(),
	assay_positive_partitions = double()
)


control_df <- data.frame(
  assay_id = character(),
  assay_type = character(),
  sample_id = character(),
  assay_well = character(),
  assay_block_id = character(),
  assay_batch_id = character(),
  assay_target_name = character(),
  assay_target_abbreviation = character(),
  assay_target_fluorophore = character(),
  assay_target_copies_per_ul_reaction = double(),
  assay_accepted_partitions = double(),
  assay_positive_partitions = double()
)

# Loop over the new quantification batches in this update.
#for (i in nrow(update_df)) {
	i = 1
	# Read the batch Excel file into a named list, with each sheet as a df.
  this_fn <- update_df$file_name[i]
  this_fpath <- paste0(UPDIR, "/BATCH_FILES/", this_fn, sep="")
	platef_in <- excel2df(this_fpath)
	
	# Process the Metadata sheet.
	# This sheet contains a two column table, with keys in column 1 and vals in column 2.
	metadata_pre_df <- as.data.frame(platef_in[["Metadata"]])
	metadata_df <- as.data.frame(t(metadata_pre_df[, c(1,2)]))

	batchup_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(batchup_df) <- as.character(batchup_df[1, ])
	batchup_df <- batchup_df[, -c(1)]
	batchup_df <- batchup_df[-c(1), ]
	#batchup_df <- batchup_df %>% slice(1:1)
	#batchup_df <- batchup_df %>% select(where(~!all(is.na(.))))
	
	# Cleanup the keys so they are consistent.
	batch_keys <- colnames(batchup_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- str_replace_all(batch_keys, "\\)", "")
	batch_keys <- str_replace_all(batch_keys, "\\(", "")
	batch_keys <- str_replace_all(batch_keys, "µ", "u")
	batch_keys <- tolower(batch_keys)
	batch_keys <- str_replace_all(batch_keys, "^", "quantification_batch_")
	batch_keys <- str_replace_all(batch_keys, "^quantification_batch_batch_", "quantification_batch_")
	batch_keys <- str_replace_all(batch_keys, "^quantification_batch_type", "batch_type")
	colnames(batchup_df) <- batch_keys
	
	# Date munging, of course.
	# This function is part of the openxlsx package. It should take an Excel DateTime (which is numeric)
	# and convert it to a POSIXct object in the format YYY-MM-DD.
	batchup_df$quantification_batch_date <- convertToDateTime(as.numeric(batchup_df$quantification_batch_date))
	
	# Append the batch metadata to the main batch df.
	batch_df <- rbind(batch_df, batchup_df)
	
	# Extract the batch id for easier use later.
	batch_id <- batchup_df$quantification_batch_id[1]

	# Each batch consists of one or more assay blocks, each in its own Excel sheet whose 
	# name begins with "ASSAY". Loop through all of these blocks and update the block_df 
	# and run_df dataframes.
	for (pname in names(platef_in)) {
		print(pname)
#		pname <- "ASSAY Respiratory Multiplex"
		if (!startsWith(pname, "ASSAY")) {
			next
		}
		block_plate_df <- platef_in[[pname]]

		# Extract the block metadata.
		block_metadata_pre_df <- block_plate_df[c(1:5), c(2:3)] 
		block_metadata_df <- as.data.frame(t(block_metadata_pre_df[, c(1,2)]))
	
		blockup_df <- rownames_to_column(block_metadata_df, var = "V0")
		colnames(blockup_df) <- as.character(blockup_df[1, ])
		blockup_df <- blockup_df[, -c(1)]
		blockup_df <- blockup_df[-c(1), ]
		#blockup_df <- blockup_df %>% slice(1:1)
		#blockup_df <- blockup_df %>% select(where(~!all(is.na(.))))
		
		# Cleanup the keys so they are consistent.
		block_keys <- colnames(blockup_df)
		block_keys <- str_replace_all(block_keys, " ", "_")
		block_keys <- str_replace_all(block_keys, "\\)", "")
		block_keys <- str_replace_all(block_keys, "\\(", "")
		block_keys <- str_replace_all(block_keys, "µ", "u")
		block_keys <- tolower(block_keys)
		block_keys <- str_replace_all(block_keys, "^", "assay_block_")
		block_keys <- str_replace_all(block_keys, "assay_block_assay_", "assay_block_")
		colnames(blockup_df) <- block_keys
		
		# Add a block id, which is just the assay method, and the batch id.
		blockup_df$assay_block_id <- blockup_df$assay_block_method
		blockup_df$assay_block_batch_id <- batch_id
		
		# Append the block metadata to the main block df.
		block_df <- rbind(block_df, blockup_df)
		
		# Extract the block id for easier use later.
		block_id <- blockup_df$assay_block_id[1]

		# Extract the block target data.
		block_target_df <- block_plate_df[c(8:13), c(2:6)]
		colnames(block_target_df) <- tolower(block_target_df[1, ])
		block_target_df <- block_target_df[-1, ]
		block_target_keys <- colnames(block_target_df)
		block_target_keys <- str_replace_all(block_target_keys, " ", "_")
		colnames(block_target_df) <- block_target_keys
		block_target_df <- block_target_df %>% filter(!if_all(everything(), is.na))

		# Extract the block plate map and assign sample to reaction well.
		plate_map_df <- block_plate_df[c(17:24), c(2:13)] 
		colnames(plate_map_df) <- c(1:12)
		rownames(plate_map_df) <- LETTERS[1:8]
		
		for (i in 1:nrow(plate_map_df)) {
			well_r <- toupper(rownames(plate_map_df)[i])
			for (j in 1:ncol(plate_map_df)) {
				well_c <- colnames(plate_map_df)[j]
				if (as.numeric(well_c) < 10) {
					well_c <- paste0(0, well_c, sep="")
				}
				sample_id <- plate_map_df[i, j]
				if (is.na(sample_id) | sample_id == "") {
					next
				}
				well_id <- paste0(well_r, well_c, sep="")
				
				# Each reaction well contains a number of assays determined by the block target table.
				# Those are keyed on target in the block df, target_abbreviation in the 
				# dPCR result file, and fluorophore in the ddPCR result file. Use targ2rxn_df 
				# to map these together.
				# Every assay is added as a unique row into the run dataframe.
				
				for (target_name in block_target_df$target) {
				
					# Need a unique string to create an assay or control id.
					# This uses the sample id and the current epoch time (seconds since 1/1/1970).
					# It ensures that even if a sample was run multiple times on the same plate, 
					# each assay of that sample will receive a UID.
					epoch_time <- as.numeric(now())
						
					my_target <- block_target_df %>% filter(target == target_name)
					my_targ_map <- targ2mix_df %>% filter(target == target_name)
					targ_abbrev <- my_targ_map$target_abbreviation[1]
					targ_fluor <- my_target$fluorophore[1]
					
					if (str_detect(well_id, my_target$negative_control_wells)) {
						# Add to controls df as a negative ctl.
						control_df <- control_df %>% add_row(
							assay_id = paste0("pnc.", epoch_time, sep=""),
							sample_id = "CTL",
							assay_type = "pcr_negative_control",
							assay_well = well_id,
							assay_block_id = block_id,
							assay_batch_id = batch_id,
							assay_target_name = target_name,
							assay_target_abbreviation = targ_abbrev,
							assay_target_fluorophore = targ_fluor,
							assay_target_copies_per_ul_reaction = -1,
							assay_accepted_partitions = -1,
							assay_positive_partitions = -1
						)
					} else if (str_detect(well_id, my_target$positive_control_wells)) {
						# Add to controls df as a positive ctl.
						control_df <- control_df %>% add_row(
						  assay_id = paste0("ppc.", epoch_time, sep=""),
						  sample_id = "CTL",
						  assay_type = "pcr_positive_control",
						  assay_well = well_id,
						  assay_block_id = block_id,
						  assay_batch_id = batch_id,
						  assay_target_name = target_name,
						  assay_target_abbreviation = targ_abbrev,
						  assay_target_fluorophore = targ_fluor,
						  assay_target_copies_per_ul_reaction = -1,
						  assay_accepted_partitions = -1,
						  assay_positive_partitions = -1
						)
					} else {
						# Add to the run df as an assay.
						
						# The columns here must match the columns in run_df. (We don't have the actual 
						# results for this assay yet, so we'll just add blanks for the result data and 
						# add it later.)
						run_df <- run_df %>% add_row(
							assay_id = paste0(sample_id, ".", epoch_time, sep=""),
							sample_id = sample_id,
							assay_type = "sample",
							assay_well = well_id,
							assay_block_id = block_id,
							assay_batch_id = batch_id,
							assay_target_name = target_name,
							assay_target_abbreviation = targ_abbrev,
							assay_target_fluorophore = targ_fluor,
							assay_target_copies_per_ul_reaction = -1,
							assay_accepted_partitions = -1,
							assay_positive_partitions = -1
						)
					}
				}
			}
		}
	}
#}

if (nrow(run_df) > 0) {

	run_df <- run_df %>% mutate(
		validation_key = case_when(
			is.na(assay_well) | assay_well == "" ~ -1, 
			is.na(sample_id) | sample_id == "" ~ -1, 
			.default = 1
		)
	)

	# Print the run, control, block, and batch df only if there is run data.
	dfn <- paste0(UPDIR, "/update.quant_run.txt", sep="")
	write.table(run_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	cfn <- paste0(UPDIR, "/update.quant_control.txt", sep="")
	write.table(control_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	lfn <- paste0(UPDIR, "/update.quant_block.txt", sep="")
	write.table(block_df, file = lfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	bfn <- paste0(UPDIR, "/update.qbatch.txt", sep="")
	write.table(batch_df, file = bfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)


}

# Return the number of rows in the run df so the controller script can do the right thing.
ret_val <- nrow(run_df)
if (ret_val > 0) {
	if (nrow(run_df %>% filter(validation_key == -1)) > 0) {
		ret_val <- -1 * ret_val
	}
}

#cat(ret_val)

