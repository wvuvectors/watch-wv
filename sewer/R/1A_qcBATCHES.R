#! /usr/bin/env Rscript

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 argument on STDIN:
#   INDIR is the directory that contains the run data.
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	INDIR <- args[1]
}
#INDIR <- "test/INPUT_FILES/2026-42-43"
RUNID <- basename(INDIR)

now <- Sys.time()

qcsum_df <- data.frame(
  batch_id = character(),
  batch_type = character(),
  batch_subtype = character(),
  batch_date = character(),
  num_samples_in_batch = integer(),
  num_validated_samples_in_batch = integer(),
  unknown_sample_list = character(),
  valid_batch_id = character(),
  valid_batch_date = character(),
  valid_sample_set = character(),
  valid_controls_pos = character(),
  valid_controls_neg = character()
)

batch_names <- c("Sample", "Concentration", "Extraction", "Quantification")
batch_list <- setNames(vector("list", length(batch_names)), batch_names)


# Loop over the Excel files in this update.
#
runfile_list <- list.files(
	path = INDIR,
	pattern = "\\.xlsx|m$",
	full.names = TRUE
)

for (fpath in runfile_list) {
	# Read the batch file into a named list, with each sheet as a df.
	platef_in <- excel2df(fpath)
	
	# Test if this is an AssetTiger (Sample) file, or a batch file.
	#
	# Read the first sheet.
	sheet1_df <- as.data.frame(platef_in[1])
	colnames(sheet1_df) <- as.character(sheet1_df[1, ])
	sheet1_df <- sheet1_df[-1, ]

	if ("Asset Tag ID" %in% colnames(sheet1_df)) {
		# This is a Sample file.
		
		# Cleanup the keys so they are consistent.
		batch_keys <- colnames(sheet1_df)
		batch_keys <- str_replace_all(batch_keys, " ", "_")
		batch_keys <- str_replace_all(batch_keys, "\\.", "")
		batch_keys <- str_replace_all(batch_keys, "\\)", "")
		batch_keys <- str_replace_all(batch_keys, "\\(", "")
		batch_keys <- str_replace_all(batch_keys, "\\[", "")
		batch_keys <- str_replace_all(batch_keys, "\\]", "")
		batch_keys <- str_replace_all(batch_keys, "µ", "u")
		batch_keys <- str_replace_all(batch_keys, "/", "_per_")
		batch_keys <- tolower(batch_keys)
		colnames(sheet1_df) <- batch_keys
		
		# Need a list of Asset IDs to check against the run batches for missing or invalid IDs.
		#
		s_df <- sheet1_df %>% select(
			sample_id = asset_tag_id,
			sample_qc_check,
			event_type
		)
		batch_list[["Sample"]] <- s_df
		
		#########################################
		qcsum_df <- qcsum_df %>% add_row(
			batch_id = "NA",
			batch_type = "Sample",
			batch_subtype = "NA",
			batch_date = "NA",
			num_samples_in_batch = nrow(sheet1_df),
			num_validated_samples_in_batch = nrow(sheet1_df),
			unknown_sample_list = "",
			valid_batch_id = "NA",
			valid_batch_date = "NA",
			valid_sample_set = "PASS",
			valid_controls_pos = "NA",
			valid_controls_neg = "NA"
		)
		#########################################

	} else if ("Metadata" %in% names(platef_in)) {
		# Any valid batch Excel file will have a valid Metadata sheet.

		# Read the Metadata sheet.
		# This sheet contains a two column table, with keys in column 1 and vals in column 2.
		metadata_pre_df <- as.data.frame(platef_in[["Metadata"]])
		metadata_df <- as.data.frame(t(metadata_pre_df[, c(1,2)]))
		metadata_df <- metadata_df %>% select(where(~!all(is.na(.))))
		
		batchup_df <- rownames_to_column(metadata_df, var = "V0")
		colnames(batchup_df) <- as.character(batchup_df[1, ])
		batchup_df <- batchup_df[, -c(1)]
		batchup_df <- batchup_df[-c(1), ]
	
		# Cleanup the keys so they are consistent.
		batch_keys <- colnames(batchup_df)
		batch_keys <- str_replace_all(batch_keys, " ", "_")
		batch_keys <- str_replace_all(batch_keys, "\\.", "")
		batch_keys <- str_replace_all(batch_keys, "\\)", "")
		batch_keys <- str_replace_all(batch_keys, "\\(", "")
		batch_keys <- str_replace_all(batch_keys, "\\[", "")
		batch_keys <- str_replace_all(batch_keys, "\\]", "")
		batch_keys <- str_replace_all(batch_keys, "µ", "u")
		batch_keys <- str_replace_all(batch_keys, "/", "_per_")
		batch_keys <- tolower(batch_keys)
		colnames(batchup_df) <- batch_keys
		
		this_type <- batchup_df$batch_type[1]
		this_id <- batchup_df$batch_id[1]
		this_date <- batchup_df$date[1]
		
		# This function is part of the openxlsx package. It should take an Excel DateTime (which is numeric)
		# and convert it to a POSIXct object in the format YYY-MM-DD.
		this_date_obj <- convertToDateTime(as.numeric(this_date))
		this_date_str <- format(this_date_obj, "%m-%d-%Y")
	
		
		if (difftime(now, this_date_obj, units = "weeks") < 2) {
			this_date_check <- "PASS"
		} else {
			this_date_check <- "FAIL"
		}
		
		qcsum_df <- qcsum_df %>% add_row(
			batch_id = this_id,
			batch_type = this_type,
			batch_subtype = "NA",
			batch_date = this_date_str,
			num_samples_in_batch = NA,
			num_validated_samples_in_batch = NA,
			unknown_sample_list = "",
			valid_batch_id = if_else(is.na(batchup_df$batch_id[1]), "FAIL", "PASS", "FAIL"),
			valid_batch_date = this_date_check,
			valid_sample_set = "NA",
			valid_controls_pos = "NA",
			valid_controls_neg = "NA"
		)
	}
}

fn <- paste0(INDIR, "/_qcbatches.", RUNID, ".summary.txt", sep="")
write.table(qcsum_df, file = fn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

identical(batch_names, qcsum_df$batch_type)

ret_val <- 0
if (!identical(batch_names, qcsum_df$batch_type)) {
	ret_val <- ret_val + 1
}
if (as.numeric(qcsum_df %>% filter(valid_batch_id == "FAIL") %>% count()) > 0) {
	ret_val <- ret_val + 2
}

# 0 = no errors
# 1 = mssing at least one required batch file
# 2 = batch ID is missing
# 3 = both missing batch id and batch file
#
cat(ret_val)
