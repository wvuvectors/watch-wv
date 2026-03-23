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

f <- file("stdin")
open(f)
while(length(line <- readLines(f, n = 1)) > 0) {
 fpaths <- strsplit(line, " ")[[1]]
}
close(f)

#fpaths <- c("../sewer/tmp/")
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
update_df <- input_df %>% filter(batch_type == "assay")

# These are extracted from the results file associated with the current APLATE.
# assay_target_copies_per_ul_reaction
# assay_accepted_droplets
# assay_positive_droplets
run_df <- data.frame(
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


# Loop over the new concentration batches in this update.
for (i in nrow(update_df)) {
	# Read the batch Excel file into a named list, with each sheet as a df.
  this_fn <- update_df$file_name[i]
  this_fpath <- paste0(UPDIR, "/batches/", this_fn, sep="")
	platef_in <- excel2df(this_fpath)
	
	# Extract each data sheet to make life easier.
	metadata_df <- as.data.frame(t(platef_in$Metadata))

	targets_df <- platef_in$Targets
	targets_df <- targets_df %>% select(where(~!all(is.na(.))))

	plate_df <- platef_in$Plate_Map %>% column_to_rownames(var = "...1")
	plate_df <- plate_df %>% select(where(~!all(is.na(.))))
	
	comment_df <- platef_in$Comment_Map %>% column_to_rownames(var = "...1")
	storage_df <- platef_in$Storage_Map %>% column_to_rownames(var = "...1")
	vol_override_df <- platef_in$Volume_Override_Map %>% column_to_rownames(var = "...1")
	
	# Process the batch metadata.
	# This sheet contains a two column table, with keys in column 1 and vals in column 2.
	batchup_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(batchup_df) <- as.character(batchup_df[1, ])
	batchup_df <- batchup_df[-c(1:1), ]
	# We just want the first two rows of the transposed metadata.
	batchup_df <- batchup_df %>% slice(1:1)
	# Remove any columns with only NAs.
	batchup_df <- batchup_df %>% select(where(~!all(is.na(.))))
	# Simplify the column headers.
	batch_keys <- colnames(batchup_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- tolower(batch_keys)
	batch_keys <- str_replace_all(batch_keys, "^", "assay_")
	batch_keys <- str_replace_all(batch_keys, "^assay_assay_", "assay_")
	colnames(batchup_df) <- batch_keys
	# Format the data correctly.
	batchup_df$assay_date <- convertToDateTime(as.numeric(batchup_df$assay_date))
	
	# Add this plate's metadata to the main batch dataframe.
	batchup_df <- batchup_df %>% select(
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
	batch_id <- batchup_df$assay_batch_id[1]

	# Process the targets for this plate.
	target_keys <- colnames(target_df)
	target_keys <- str_replace_all(target_keys, " ", "_")
	target_keys <- tolower(target_keys)
	target_keys <- str_replace_all(target_keys, "^", "assay_")
	target_keys <- str_replace_all(target_keys, "^assay_assay_", "assay_")
	colnames(target_df) <- target_keys
	
	# Init a run update table.
	runup_df <- data.frame(
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
				input_v <- batchup_df[1, "assay_input_ul"]
				if (!is.na(vol_override_df[i, j])) {
					input_v <- vol_override_df[i, j]
				}
				
				# Need a unique string to create an assay id.
				# This uses the sample id and the current epoch time (seconds since 1/1/1970).
				# It ensures that even if a sample was run multiple times on the same plate, each will 
				# receive a UID.
				#date_str <- str_replace_all(batchup_df$concentration_date[1], "-", "")
				epoch_time <- as.numeric(now())

				this_assay_id <- paste0(sample_id, ".", epoch_time, ".", as.numeric(k), sep="")
				this_control_id <- paste0(batch_id, ".", epoch_time, ".", target_df[k, "assay_target_fluorophore"], ".", ctl_incr, sep="")

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
					runup_df <- add_row(
						runup_df,
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
			assay_result_df <- result_df %>% filter(assay_well %in% runup_df$assay_well)
			if (nrow(assay_result_df)) {
				runup_df <- left_join(runup_df, assay_result_df, by = c("assay_well", "assay_target_fluorophore"))
			}
			control_result_df <- result_df %>% filter(assay_well %in% controlup_df$assay_well)
			if (nrow(control_result_df)) {
				controlup_df <- left_join(controlup_df, control_result_df, by = c("assay_well", "assay_target_fluorophore"))
			}
		}
	}

	
	run_df <- rbind(run_df, runup_df)
	if (nrow(controlup_df)) {
		control_df <- rbind(control_df, controlup_df)
	}
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
			is.na(assay_batch_id) | assay_batch_id == "" ~ -1, 
			is.na(assay_input_ul) | assay_input_ul == "" ~ -1, 
			is.na(assay_target) | assay_target == "" ~ -1, 
			is.na(assay_target_genetic_locus) | assay_target_genetic_locus == "" ~ -1, 
			is.na(assay_target_fluorophore) | assay_target_fluorophore == "" ~ -1, 
			is.na(assay_target_copies_per_ul_reaction) | assay_target_copies_per_ul_reaction == "" ~ -1, 
			.default = 1
		)
	)
	
	batch_df <- batch_df %>% mutate(
		validation_key = case_when(
			is.na(assay_batch_id) | assay_batch_id == "" ~ -1, 
			is.na(assay_reaction_ul) | assay_reaction_ul == "" ~ -1, 
			is.na(assay_amplification_method) | assay_amplification_method == "" ~ -1, 
			is.na(assay_method) | assay_method == "" ~ -1, 
			.default = 1
		)
	)

	# Print the run and batch df only if there is run data.
	rfn <- paste0(UPDIR, "/update.assay.txt", sep="")
	write.table(run_df, file = rfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	bfn <- paste0(UPDIR, "/update.abatch.txt", sep="")
	write.table(batch_df, file = bfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	cfn <- paste0(UPDIR, "/update.control.txt", sep="")
	write.table(control_df, file = cfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
}

# Return the number of rows in the run df so the controller script can do the right thing.
ret_val <- nrow(run_df)
if (ret_val > 0) {
	if (nrow(run_df %>% filter(validation_key == -1)) > 0 | nrow(batch_df %>% filter(validation_key == -1)) > 0) {
		ret_val <- -1 * ret_val
	}
}

cat(ret_val)

