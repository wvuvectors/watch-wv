#! /usr/bin/env Rscript --vanilla

library(tidyverse)
library(dplyr)
library(data.table)
library(DT)
library(zoo)
library(rlang)
library(glue)

library(readxl)
library(openxlsx)

library(scales)
library(lubridate)

excel2df <- function(fname) { 
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x, col_types = "text"))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # return the data frame
  data_frame
} 

f <- file("stdin")
open(f)
while(length(line <- readLines(f, n = 1)) > 0) {
 fpaths <- strsplit(line, " ")[[1]]
}
close(f)

#fpaths <- c("../patchr/tmp/")
UPDIR <- fpaths[1]

input_df <- as.data.frame(
	read.table(
		paste0(UPDIR, "/update.batch_files.txt", sep=""), 
		quote="", 
		sep="\t", 
		header=TRUE, 
		check.names=FALSE)
)

update_df <- input_df %>% filter(batch_type == "assay")

# These are extracted from the results file associated with the current APLATE.
# assay_target_copies_per_ul_reaction
# assay_accepted_droplets
# assay_positive_droplets
assay_df <- data.frame(
  assay_id = character(),
  sample_id = character(),
  assay_well = character(),
  assay_batch_id = character(),
  assay_input_ul = character(),
  assay_comment = character(),
	assay_target = character(),
	assay_target_genetic_locus = character(),
	assay_target_template = character(),
	assay_target_macromolecule = character(),
	assay_target_fluorophore = character(),
	assay_target_copies_per_ul_reaction = character(),
	assay_accepted_droplets = character(),
	assay_positive_droplets = character()
)
control_df <- data.frame(
	control_id = character(),
	control_type = character(),
	assay_well = character(),
	assay_batch_id = character(),
	assay_input_ul = character(),
	assay_comment = character(),
	assay_target = character(),
	assay_target_genetic_locus = character(),
	assay_template = character(),
	assay_target_macromolecule = character(),
	assay_target_fluorophore = character(),
	assay_target_copies_per_ul_reaction = character(),
	assay_accepted_droplets = character(),
	assay_positive_droplets = character()
)
batch_df <- data.frame(
  assay_batch_id = character(),
  assay_batch_type = character(),
  assay_batch_record_version = character(),
  assay_date = character(),
  assay_run_by = character(),
  assay_machine = character(),
  assay_analysis_software_version = character(),
  assay_reaction_ul = character(),
  assay_amplification_method = character(),
  assay_amplification_method_lot_id = character(),
  assay_quantification_type = character(),
  assay_method = character(),
  assay_method_lot_id = character(),
  assay_batch_comment = character()
)

for (i in nrow(update_df)) {
  this_fn <- update_df$file_name[i]
  this_fpath <- paste0(UPDIR, "/batches/", this_fn, sep="")

	platef_in <- excel2df(this_fpath)
	plate_df <- platef_in$Plate_Map %>% column_to_rownames(var = "...1")
	plate_df <- plate_df %>% select(where(~!all(is.na(.))))
	
	comment_df <- platef_in$Comment_Map %>% column_to_rownames(var = "...1")
	storage_df <- platef_in$Storage_Map %>% column_to_rownames(var = "...1")
	vol_override_df <- platef_in$Volume_Override_Map %>% column_to_rownames(var = "...1")
	
	# Read in the metadata for this plate and transpose it.
	metadata_df <- as.data.frame(t(platef_in$Metadata))
	metadata_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(metadata_df) <- as.character(metadata_df[1, ])
	metadata_df <- metadata_df[-c(1:1), ]
	# We just want the first two rows of the transposed metadata.
	metadata_df <- metadata_df %>% slice(1:1)
	# Remove any columns with only NAs.
	metadata_df <- metadata_df %>% select(where(~!all(is.na(.))))
	# Simplify the column headers.
	batch_keys <- colnames(metadata_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- tolower(batch_keys)
	batch_keys <- str_replace_all(batch_keys, "^", "assay_")
	batch_keys <- str_replace_all(batch_keys, "^assay_assay_", "assay_")
	colnames(metadata_df) <- batch_keys
	# Format the data correctly.
	metadata_df$assay_date <- convertToDateTime(as.numeric(metadata_df$assay_date))
	
	# Add this plate's metadata to the main batch dataframe.
	batchup_df <- metadata_df %>% select(
		assay_batch_id,
		assay_batch_type,
		assay_batch_record_version,
		assay_date,
		assay_run_by,
		assay_machine,
		assay_analysis_software_version,
		assay_reaction_ul,
		assay_amplification_method,
		assay_amplification_method_lot_id,
		assay_quantification_type,
		assay_method,
		assay_method_lot_id,
		assay_batch_comment
	)
	batch_df <- rbind(batch_df, batchup_df)
	
	# Extract the batch id for easier use later.
	batch_id <- metadata_df$assay_batch_id[1]
	# Create a suffix for assay ids based on the plate run date.
	date_str <- str_replace_all(metadata_df$assay_date[1], "-", "")

	# Read in the targets for this plate.
	target_df <- as.data.frame(platef_in$Targets)
	# Simplify the column headers.
	target_keys <- colnames(target_df)
	target_keys <- str_replace_all(target_keys, " ", "_")
	target_keys <- tolower(target_keys)
	target_keys <- str_replace_all(target_keys, "^", "assay_")
	target_keys <- str_replace_all(target_keys, "^assay_assay_", "assay_")
	colnames(target_df) <- target_keys
	
	# Init the assay update table.
	assayup_df <- data.frame(
		assay_id = character(),
		sample_id = character(),
		assay_well = character(),
		assay_batch_id = character(),
		assay_input_ul = character(),
		assay_comment = character(),
		assay_target = character(),
		assay_target_genetic_locus = character(),
		assay_template = character(),
		assay_target_macromolecule = character(),
		assay_target_fluorophore = character()
	)
#		assay_target_copies_per_ul_reaction = character(), # Add this col in the results loop
#		assay_accepted_droplets = character(), # Add this col in the results loop
#		assay_positive_droplets = character() # Add this col in the results loop
	
	# Save any control wells into a separate df to update the main controls_df.
	controlup_df <- data.frame(
		control_id = character(),
		control_type = character(),
		assay_well = character(),
		assay_batch_id = character(),
		assay_input_ul = character(),
		assay_comment = character(),
		assay_target = character(),
		assay_target_genetic_locus = character(),
		assay_template = character(),
		assay_target_macromolecule = character(),
		assay_target_fluorophore = character()
	)
#		assay_target_copies_per_ul_reaction = character(), # Add this col in the results loop
#		assay_accepted_droplets = character(), # Add this col in the results loop
#		assay_positive_droplets = character() # Add this col in the results loop

	
	# Loop through the targets in target_df and combine with the sample ids from plate_df
	for (k in 1:nrow(target_df)) {
		# Init a few incremental counters.
		ctl_incr <- 1

		# Get the assay target and target locus.
	  targ_combo <- target_df[k, "assay_target"]
	  patt <- "^(.+) \\((.+)\\)"
	  target_vec <- str_match(targ_combo, patt)
	  this_target <- targ_combo
	  this_locus <- "Unk" 
	  if (length(target_vec) > 1 & !is.na(target_vec[2])) {
	    this_target <- target_vec[2]
	  }
	  if (length(target_vec) > 2 & !is.na(target_vec[3])) {
	    this_locus <- target_vec[3]
	  }
	  
	  negative_controls <- c()
	  if (!is.na(target_df[k, "assay_negative_control_wells"])) {
      nstr <- target_df[k, "assay_negative_control_wells"]
	    negative_controls <- unlist(str_split(nstr, ","))
	  }
	  positive_controls <- c()
	  if (!is.na(target_df[k, "assay_positive_control_wells"])) {
      pstr <- target_df[k, "assay_positive_control_wells"]
	    positive_controls <- unlist(str_split(pstr, ","))
	  }
	  
	  # Each cell of the plate_df contains a sample id. The plate well is a combination of 
	  # the row and column names. That well corresponds to the results file, so we need it.
	  for (i in 1:nrow(plate_df)) {
			well_r <- toupper(rownames(plate_df)[i])
			for (j in 1:ncol(plate_df)) {
				well_c <- colnames(plate_df)[j]
				# This block ensures the well label is formatted so it aligns with the results file, 
				# by standardizing row capitalization and making all column numbers exactly two digits.
				if (as.numeric(well_c) < 10) {
					well_c <- paste0(0, well_c, sep="")
				}
				well_id <- paste0(well_r, well_c, sep="")
				
				# Get the sample id for this well.
				sample_id <- plate_df[i, j]
				
				# If there is a custom input vol for this well, get it from the vol override plate. 
				# Otherwise we use the default from the metadata plate.
				input_v <- metadata_df[1, "assay_input_ul"]
				if (!is.na(vol_override_df[i, j])) {
					input_v <- vol_override_df[i, j]
				}
				
				this_assay_id <- paste0(sample_id, ".", date_str, ".", as.numeric(k), sep="")
				this_control_id <- paste0(batch_id, ".", date_str, ".", target_df[k, "assay_target_fluorophore"], ".", ctl_incr, sep="")
				# If this well is a control well, store it in the controls df. Otherwise, store it 
				# in the assay df.
				if (well_id %in% negative_controls) {
					controlup_df <- add_row(
						controlup_df,
						control_id = this_control_id,
						control_type = "negative",
						assay_well = well_id,
						assay_batch_id = batch_id,
						assay_input_ul = input_v,
						assay_comment = comment_df[i, j],
						assay_target = this_target,
						assay_target_genetic_locus = this_locus,
						assay_template = target_df[k, "assay_negative_control_template"],
						assay_target_macromolecule = target_df[k, "assay_negative_control_macromolecule"],
						assay_target_fluorophore = target_df[k, "assay_target_fluorophore"]
					)
					ctl_incr <- ctl_incr+1
				} else if (well_id %in% positive_controls) {
					controlup_df <- add_row(
						controlup_df,
						control_id = this_control_id,
						control_type = "positive",
						assay_well = well_id,
						assay_batch_id = batch_id,
						assay_input_ul = input_v,
						assay_comment = comment_df[i, j],
						assay_target = this_target,
						assay_target_genetic_locus = this_locus,
						assay_template = target_df[k, "assay_positive_control_template"],
						assay_target_macromolecule = target_df[k, "assay_positive_control_macromolecule"],
						assay_target_fluorophore = target_df[k, "assay_target_fluorophore"]
					)
					ctl_incr <- ctl_incr+1
				} else {
					assayup_df <- add_row(
						assayup_df,
						assay_id = this_assay_id,
						sample_id = sample_id,
						assay_well = well_id,
						assay_batch_id = batch_id,
						assay_input_ul = input_v,
						assay_comment = comment_df[i, j],
						assay_target = this_target,
						assay_target_genetic_locus = this_locus,
						assay_template = "sample", 
						assay_target_macromolecule = target_df[k, "assay_target_macromolecule"],
						assay_target_fluorophore = target_df[k, "assay_target_fluorophore"]
					)
				}
			}
		}
	}
	
	# Get the results from the result file for this batch. Use these data to complete the 
	# assayup_df and controlup_df tables.
	csv_pattern <- paste0(batch_id, "_.+\\.csv$", sep="")
	csv_files <- list.files(
		path = paste0(UPDIR, "/batches/", sep=""),
		pattern = csv_pattern,
		full.names = TRUE
	)
	if (length(csv_files) > 0) {
		for (f in csv_files) {
			result_df <- as.data.frame(
				read.csv(
					csv_files[1], 
					quote="", 
					header=TRUE, 
					check.names=FALSE)
			)
			result_df <- result_df %>% select(
				assay_well = `Well`, 
				assay_target_fluorophore = `DyeName(s)`, 
				assay_target_copies_per_ul_reaction = `Conc(copies/uL)`, 
				assay_accepted_droplets = `Accepted Droplets`, 
				assay_positive_droplets = `Positives`
			)
			assay_result_df <- result_df %>% filter(assay_well %in% assayup_df$assay_well)
			if (nrow(assay_result_df)) {
				assayup_df <- left_join(assayup_df, assay_result_df, by = c("assay_well", "assay_target_fluorophore"))
			}
			control_result_df <- result_df %>% filter(assay_well %in% controlup_df$assay_well)
			if (nrow(control_result_df)) {
				controlup_df <- left_join(controlup_df, control_result_df, by = c("assay_well", "assay_target_fluorophore"))
			}
		}
	}

	
	assay_df <- rbind(assay_df, assayup_df)
	if (nrow(controlup_df)) {
		control_df <- rbind(control_df, controlup_df)
	}
}

if (nrow(batch_df) > 0) {
	bfn <- paste0(UPDIR, "update.abatch.txt", sep="")
	write.table(batch_df, file = bfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	dfn <- paste0(UPDIR, "update.assay.txt", sep="")
	write.table(assay_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	cfn <- paste0(UPDIR, "update.control.txt", sep="")
	write.table(control_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

}

ret_val <- nrow(assay_df)
cat(ret_val)

