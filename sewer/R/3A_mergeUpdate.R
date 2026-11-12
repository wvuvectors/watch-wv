#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")


# Accepts 3 directory paths on STDIN:
#   UPDIR is the readable directory that contains the update files for the run being processed.
#		DBDIR is the readable directory containing the current tables, which are going to be updated.
#   MERGEDIR is the writeable destination directory for the merged output.
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	UPDIR <- args[1]
	DBDIR <- args[2]
	MERGEDIR <- args[3]
}
#UPDIR <- "../data/updates/2026-10-15"
#DBDIR <- "../data/latest"
#MERGEDIR <- "../data/merges/2026-10-15"


# A vector of DB file suffixes, to make looping easier.
suffix_vec <- c(
	"cbatch", "concentration", 
	"ebatch", "extraction",
	"qbatch", "qblock", "qcontrol", "qassay", 
	"sample"
)


flagdups_vec <- c("concentration", "extraction", "sample")

duplicates_df <- data.frame(
	source_table = character(),
	table_name = character(),
	sample_id = character()
)

# Read in each matching pair of tables from UPDIR and DBDIR.
for (i in 1:length(suffix_vec)) {
#  table_suffix <- "cbatch"
  table_suffix <- suffix_vec[i]
  latest_df <- as.data.frame(
    read.table(
      paste0(DBDIR, "/watchdb.", table_suffix , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      fill=TRUE,
      comment.char="",
      check.names=FALSE)
  )
  latest_df <- latest_df %>% mutate(across(everything(), as.character))

  update_df <- as.data.frame(
    read.table(
      paste0(UPDIR, "/update.", table_suffix , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      fill=TRUE,
      comment.char="",
      check.names=FALSE)
  )
  update_df <- update_df %>% mutate(across(everything(), as.character))
  
  # Merge the update and the watchdb tables.
  if (nrow(latest_df) == 0) {
    merged_df <- update_df
    merged_df$source_table <- "inbound"
  } else if (nrow(update_df) == 0) {
    merged_df <- latest_df
    merged_df$source_table <- "existing"
  } else {
    merged_df <- bind_rows(existing=latest_df, inbound=update_df, .id = "source_table")
  }
  
  # If there are duplicate sample_ids in unallowed tables, flag them for verification.
  # Currently we do NOT remove duplicates.
  # Just check between existing and incoming tables; if they are already in the db, ignore.
  #
  if (table_suffix %in% flagdups_vec) {
    dups_in_update_df <- update_df %>% group_by(sample_id) %>% filter(n() > 1) %>% ungroup()
    dups_in_latest_df <- latest_df %>% group_by(sample_id) %>% filter(n() > 1) %>% ungroup()
    dups_in_both_df <- update_df %>% filter(sample_id %in% latest_df$sample_id)
    if (nrow(dups_in_both_df) > 0) {
      dups_df <- dups_in_both_df %>% select(sample_id) %>% 
      	mutate(
      		source_table = "BOTH", 
      		table_name = table_suffix
      	)
      duplicates_df <- bind_rows(duplicates_df, distinct(dups_df))
    }
    if (nrow(dups_in_update_df) > 0) {
      dups_df <- dups_in_update_df %>% select(sample_id) %>% 
      	mutate(
      		source_table = "update", 
      		table_name = table_suffix
      	)
      duplicates_df <- bind_rows(duplicates_df, distinct(dups_df))
    }
#     if (nrow(dups_in_latest_df) > 0) {
#       dups_df <- dups_in_latest_df %>% select(sample_id) %>% 
#       	mutate(
#       		source_table = "latest", 
#       		table_name = table_suffix
#       	)
#       duplicates_df <- bind_rows(duplicates_df, distinct(dups_df))
#     }
  }
  
  # Write the merge tables to MERGEDIR.
  #
  if (nrow(merged_df) > 0) {
    incrfn <- paste0(MERGEDIR, "/merge.", table_suffix, ".txt", sep="")
    write.table(merged_df, file = incrfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
  }
  
}

if (nrow(duplicates_df) > 0) {
  dupfn <- paste0(MERGEDIR, "/merge.DUPLICATE_SAMPLES.txt", sep="")
  write.table(duplicates_df, file = dupfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
}

# Need to send the duplicate count back to bash so the control script can do the right thing.
ret_val <- as.numeric(nrow(duplicates_df))
cat(ret_val)

