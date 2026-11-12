#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")


# Accepts 2 directory paths on STDIN:
#		DBDIR is the directory containing the current tables, which are going to be replaced.
#   MERGEDIR is the directory containing the merged output.
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	DBDIR <- args[1]
	MERGEDIR <- args[2]
}
#DBDIR <- "../data/latest"
#MERGEDIR <- "../data/merges/2026-10-15"


# A vector of table file suffixes, to make looping easier.
suffix_vec <- c(
	"cbatch", "concentration", 
	"ebatch", "extraction",
	"qbatch", "qblock", "qcontrol", "qassay", 
	"sample", "result"
)

# Read in each merge table from MERGEDIR.
for (i in 1:length(suffix_vec)) {
	table_suffix <- suffix_vec[i]
	merged_df <- as.data.frame(
		read.table(
			paste0(MERGEDIR, "/merge.", table_suffix , ".txt", sep=""), 
			quote="", 
			sep="\t", 
			header=TRUE, 
			fill=TRUE,
			comment.char="",
			check.names=FALSE)
	)
	if ("source_table" %in% names(merged_df)) {
		merged_df <- merged_df %>% select(-source_table)
	}
	merged_df <- merged_df %>% mutate(across(everything(), as.character))
	
	# Write the merge tables to DBDIR.
	#
	outfn <- file.path(DBDIR, paste0("watchdb.", table_suffix, ".txt", sep=""))
	write.table(merged_df, file = outfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	
}

ret_val <- 0
cat(ret_val)

