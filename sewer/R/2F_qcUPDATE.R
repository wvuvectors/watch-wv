#! /usr/bin/env Rscript

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 argument on STDIN:
#   UPDIR is the directory containing the update files.
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	UPDIR <- args[1]
}
#UPDIR <- "../data/updates/2026-10-22"

# Init a dataframe to hold errors.
#
errors_df <- data.frame(
	table_name = character(),
	error_priority = character(),	# warning, error
	error_detail = character(),	# sample ID missing; sample failed; sampel date misalignment
	problem_value = character()
)

# Read in the AssetTiger database to check for invalid sample IDs.
#
idcheck_df <- data.frame(
	sample_id = character(),
	sample_qc_check = character(),
	event_type = character()
)

runfile_list <- list.files(
	path = paste0(UPDIR, "/run_files", sep=""),
	pattern = "\\.xlsx|m$",
	full.names = TRUE
)

for (fpath in runfile_list) {
	# Read the batch file into a named list, with each sheet as a df.
	platef_in <- excel2df(fpath)
	
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
		idcheck_df <- sheet1_df %>% select(
			sample_id = asset_tag_id,
			sample_qc_check,
			sample_flow_mgd,
			event_type
		)
	}
}

# We also need info from the resources Excel file, specifically the targets sheet.
#
#resource_tables <- c("location", "wwtp", "county", "target", "lab", "loci", "methods")
resource_tables <- c("target")
platef_in <- excel2df("resources/watchdb.all_tables.xlsx")
resource_list <- list()
for (tname in resource_tables) {
	r_df <- as.data.frame(platef_in[[tname]])
	colnames(r_df) <- as.character(r_df[1, ])
	r_df <- r_df[-1, ]
	resource_list[[tname]] <- r_df
}
target_df <- resource_list[["target"]]



# Now read in the update tables.
#
data_tables <- c(
	"cbatch", "concentration", 
	"ebatch", "extraction",
	"qbatch", "qblock", "qcontrol", "qassay", 
	"sample"
)

data_list <- list()
for (dname in data_tables) {
	d_df <- as.data.frame(
	read.table(
    paste0(UPDIR, "/update." , dname, ".txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
		fill=TRUE,
		comment.char="", 
    check.names=FALSE)
)
	data_list[[dname]] <- d_df
}


#
# STEP 1.
# Check that the assay & control target names match the target list.
#
for (dname in c("qassay", "qcontrol")) {
	mismatched_df <- anti_join(data_list[[dname]], target_df, by = c("assay_target_name" = "target_id"))
	if (nrow(mismatched_df) > 0) {
		mismatched_df <- mismatched_df %>% 
			select(problem_value = assay_target_name) %>% 
			mutate(
				table_name = dname,
				error_priority = "urgent error",
				error_detail = "unknown target"
			)
		errors_df <- rbind(errors_df, mismatched_df)
	}
}

#
# STEP 2.
# Check that sample IDs in the data tables exist in AT.
#
for (dname in data_tables) {
	table_df <- data_list[[dname]]
	
	if (!"sample_id" %in% colnames(table_df) | dname == "qcontrol") {
	  next
	}
	
	missing_df <- anti_join(table_df, idcheck_df, by = "sample_id")
	if (nrow(missing_df) > 0) {
		missing_df <- missing_df %>% 
			select(problem_value = sample_id) %>% 
			mutate(
				table_name = dname,
				error_priority = "possible error",
				error_detail = "sample ID not in database"
			)
		errors_df <- rbind(errors_df, missing_df)
	}
}

sample_df <- data_list[["sample"]]

#
# STEP 3.
# Check for samples that failed lab QC.
#
labqc_fail_df <- sample_df %>% filter(sample_qc == "FAIL")
if (nrow(labqc_fail_df) > 0) {
	labqc_fail_df <- labqc_fail_df %>% 
		select(problem_value = sample_id) %>% 
		mutate(
			table_name = "sample",
			error_priority = "warning",
			error_detail = "sample failed lab QC"
		)
	errors_df <- rbind(errors_df, labqc_fail_df)
}

#
# STEP 4.
# Check for collection date problems (sample table).
#
sample_df$collection_start_datetime <- mdy_hm(sample_df$sample_collection_start_datetime)
sample_df$collection_end_datetime <- mdy_hm(sample_df$sample_collection_end_datetime)

dtexist_df <- sample_df %>% filter(is.na(collection_start_datetime) | is.na(collection_end_datetime))

if (nrow(dtexist_df) > 0) {

	dtexist_df <- dtexist_df %>% 
		select(problem_value = sample_id) %>% 
		mutate(
			table_name = "sample",
			error_priority = "urgent error",
			error_detail = "missing a collection start or stop timestamp"
		)
	errors_df <- rbind(errors_df, dtexist_df)

} else {

	dtcheck_df <- sample_df %>% 
		select(sample_id, collection_start_datetime, collection_end_datetime) %>% 
		mutate(collection_window = difftime(collection_end_datetime, collection_start_datetime, units = "hours"))
	
	dtcheck1_df <- dtcheck_df %>% filter(collection_window < 0)
	if (nrow(dtcheck1_df) > 0) {
		dtcheck1_df <- dtcheck1_df %>% 
			select(problem_value = sample_id) %>% 
			mutate(
				table_name = "sample",
				error_priority = "urgent error",
				error_detail = "collection end comes before collection start"
			)
		errors_df <- rbind(errors_df, dtcheck1_df)
	}
	
	dtcheck2_df <- dtcheck_df %>% filter(collection_window >= 0 & collection_window <= 2)
	if (nrow(dtcheck2_df) > 0) {
		dtcheck2_df <- dtcheck2_df %>% 
			select(problem_value = sample_id) %>% 
			mutate(
				table_name = "sample",
				error_priority = "urgent error",
				error_detail = "collection window is very small"
			)
		errors_df <- rbind(errors_df, dtcheck2_df)
	}

}


# 
# STEP 5. Check for any validation keys that failed.
#
for (dname in data_tables) {
	table_df <- data_list[[dname]]
	col_name <- "batch_validation_key"
	if (dname == "qassay" | dname == "qcontrol") {
		col_name <- "assay_validation_key"
	}
	
	if (! col_name %in% colnames(table_df)) {
	  next
	}
	
	table_df <- table_df %>% filter(.data[[col_name]] == "FAIL")
	if (nrow(table_df) > 0) {
		fails_df <- data.frame(
		  problem_value = dname, 
			table_name = dname,
			error_priority = "attention",
			error_detail = paste0("validation fail count: ", nrow(table_df), sep="")
		)
		errors_df <- rbind(errors_df, fails_df)
	}
}


# 
# STEP 6. Check for NTC amplification.
#
table_df <- data_list[["qcontrol"]]
ntc_df <- table_df %>% filter(str_detect(assay_type, "negative_control") & assay_positive_partitions >= 5)
if (nrow(ntc_df) > 0) {
	ntc_df <- ntc_df %>% 
		select(problem_value = assay_id) %>% 
		mutate(
			table_name = dname,
			error_priority = "warning",
			error_detail = "NTC value exceeds maximum"
		)
	errors_df <- rbind(errors_df, ntc_df)
}


# 
# STEP FINAL. Generate a few data plots to check controls and assays visually.
#

# for (i in 1:length(target_df$target_id)) {
# 	targ <- target_df$target_id[i]
# 	plot_df <- cleaned_df %>% filter(target == targ)
# 	
# 	if (nrow(plot_df) > 0) {
# 		p1 <- ggplot(plot_df, aes()) + 
# 			geom_point(aes(x=as.factor(epi_date), y=target_copies_flow_normalized_per_capita, color=as.factor(location_id), fill=as.factor(location_id))) + 
# 			#scale_x_continuous(limits=c(-5, 5)) + 
# 			labs(y = "Target copies/mL", title=paste0(targ, sep=""))
# 
# 		ggsave(paste0(str_replace_all(targ, "\\/", "_"), ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = paste0("analyses/", this_epi_date, sep=""))
# 	}
# }
# 
# 


ret_val <- 0
outfn <- file.path(UPDIR, "qc", "summary.txt")
noerr <- c("No errors found during update quality control.")

if (nrow(errors_df) == 0) {
	writeLines(text = noerr, con = outfn)
} else {
	# Print the error table.
	write.table(errors_df, file = outfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	ret_val <- nrow(errors_df)
}

cat(ret_val)
