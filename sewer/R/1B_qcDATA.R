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

qcsum_df <- data.frame(
	batch_id = character(),
	total_wells = integer(),
	num_pass_partitions_valid = integer(),
	pct_pass_partitions_valid = numeric(),
	valid_partition_ratio = character(),
	valid_batch_id = character()
)

# Loop over the data (csv) files in this update.
#
runfile_list <- list.files(
	path = INDIR,
	pattern = "\\.csv$",
	full.names = TRUE
)

for (fpath in runfile_list) {
	# Read the csv file into a dataframe.
	data_df <- fread(fpath, quote="")
	
	if ("Asset Tag ID" %in% names(data_df)) {
		# This is a legacy Samples csv and is ignored.
		next
	}
	
	# Cleanup the keys so they are consistent.
	batch_keys <- colnames(data_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- str_replace_all(batch_keys, "\\.", "")
	batch_keys <- str_replace_all(batch_keys, "\\)", "")
	batch_keys <- str_replace_all(batch_keys, "\\(", "")
	batch_keys <- str_replace_all(batch_keys, "\\[", "")
	batch_keys <- str_replace_all(batch_keys, "\\]", "")
	batch_keys <- str_replace_all(batch_keys, "µ", "u")
	batch_keys <- str_replace_all(batch_keys, "/", "_per_")
	batch_keys <- tolower(batch_keys)
	colnames(data_df) <- batch_keys
  
	data_df <- data_df %>% mutate(
	  prow = substr(well, 1, 1),
	  pcol = as.numeric(substr(well, 2, NULL))
	)

		# Is this a Qiacuity or QX result?
	if ("plate_name" %in% names(data_df)) {
    # This is a Qiacuity plate.
    precheck_df <- data_df %>% select(
    	partitions_valid,
    	batch_id = plate_name,
    	partitions_positive,
    	channel = target_name,
    	well,
    	prow,
    	pcol
    )
  } else {
    precheck_df <- data_df %>% select(
			sample_description_1,
			sample_description_2,
			sample_description_3,
			sample_description_4,
			partitions_valid = accepted_droplets,
			partitions_positive = positives,
			channel = dyenames,
			well, 
			prow, 
			pcol
		)
		precheck_df <- precheck_df %>% mutate(
			batch_id = case_when(
				all(!is.na(precheck_df$sample_description_1)) ~ precheck_df$sample_description_1[1],
				all(is.na(precheck_df$sample_description_1)) ~ precheck_df$sample_description_2[1],
				all(is.na(precheck_df$sample_description_2)) ~ precheck_df$sample_description_3[1],
				all(is.na(precheck_df$sample_description_3)) ~ precheck_df$sample_description_4[1],
				.default = "UNK"
			)
		)
		precheck_df <- precheck_df %>% select(!starts_with("sample_description_"))
	}
	
	# Assemble the QC row for this data file.
	#
	total_wells_on_plate <- as.integer(nrow(precheck_df))
	accepted_pass <- as.numeric(precheck_df %>% filter(as.numeric(partitions_valid) > 12000) %>% count())
	
	qcsum_df <- qcsum_df %>% add_row(
		batch_id = precheck_df$batch_id[1],
		total_wells = total_wells_on_plate,
		num_pass_partitions_valid = accepted_pass,
		pct_pass_partitions_valid = 100*(accepted_pass / total_wells_on_plate), 
		valid_partition_ratio = if_else(pct_pass_partitions_valid > 95, "PASS", "FAIL", "FAIL"), 
		valid_batch_id = if_else(precheck_df$batch_id[1] == "UNK", "FAIL", "PASS", "FAIL")
	)
	
	
	# Print some useful visualizations.
	#
	
	data_df <- complete(data_df, prow, pcol, fill = list(accepted_droplets = 0, positives = 0))
	data_df <- data_df %>% mutate(
	  pcol = if_else(pcol < 10, paste0(0, as.character(pcol), sep=""), as.character(pcol))
	)
	precheck_df <- precheck_df %>% mutate(
	  channel = factor(channel, levels = c("ATTO590", "ROX", "Cy5.5", "Cy5", "HEX", "FAM"))
  )

	p1 <- ggplot(precheck_df %>% filter(channel == "FAM"), aes()) + 
		geom_col(aes(x=as.factor(pcol), y=partitions_valid, fill=partitions_valid)) + 
	  scale_fill_continuous(
	    type = "gradient",
	    high = "gold",
	    low = "darkgreen",
	    limits = c(0, 25000),
	    guide = guide_colourbar(reverse = FALSE)
	  ) + 
	  geom_hline(yintercept = 10000, color = "black", linetype = "dashed", linewidth = 1) + 
		geom_hline(yintercept = 12000, color = "black", linetype = "dotted", linewidth = 1) + 
    facet_wrap(~ prow, nrow = 1, scales = "free_x", strip.position = "bottom") +
    theme_minimal() +
    theme(
	    panel.spacing.x = unit(0, "lines"),          # Removes spaces between groups
	    strip.placement = "outside",                 # Places the group label below the axis text
	    strip.background = element_blank(),          # Cleans up the background
	    panel.border = element_rect(fill = NA, color = "gray80"), # Optional: wraps groups in box
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
	  ) +
		labs(y = "Accepted partition count", x = NULL)

	ggsave(paste0("_qcdata.", precheck_df$batch_id[1], ".accepted_partitions.png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = INDIR)

	p2 <- ggplot(precheck_df, aes()) + 
		geom_col(aes(x=as.factor(well), y=partitions_positive, fill=partitions_positive)) + 
	  scale_fill_continuous(
	    type = "gradient",
	    high = "gold",
	    low = "darkgreen",
	    #limits = c(0, 25000),
	    guide = guide_colourbar(reverse = FALSE)
	  ) + 
	  facet_wrap(~channel, ncol = 1) +
	  theme_minimal() +
	  theme(
	    #panel.spacing.x = unit(0, "lines"),          # Removes spaces between groups
	    #strip.placement = "outside",                 # Places the group label below the axis text
	    strip.background = element_blank(),          # Cleans up the background
	    panel.border = element_rect(fill = NA, color = "gray80"), # Optional: wraps groups in box
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
	  ) +
	  labs(y = "Positive partition count", x = "Plate well")

	ggsave(paste0("_qcdata.", precheck_df$batch_id[1], ".positive_partitions.png", sep=""), plot = p2, width = 3306, height = 1892, units = c("px"), path = INDIR)

 }

fn <- paste0(INDIR, "/_qcdata.summary.txt", sep="")
write.table(qcsum_df, file = fn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

ret_val <- 0
if (as.numeric(qcsum_df %>% filter(valid_partition_ratio == "FAIL") %>% count()) > 0) {
	ret_val <- ret_val + 1
}
if (as.numeric(qcsum_df %>% filter(valid_batch_id == "FAIL") %>% count()) > 0) {
	ret_val <- ret_val + 2
}

# 0 = no errors
# 1 = too many wells failed accepted partitions test
# 2 = batch ID is missing
# 3 = both missing batch id and too many wells failed accepted partitions test
#
cat(ret_val)

