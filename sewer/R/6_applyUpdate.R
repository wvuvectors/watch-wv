#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")


# Accepts 3 directory paths on STDIN:
#   UPDIR is the readable directory that contains the update files for the run being processed.
#		DBDIR is the readable directory containing the current tables, which are going to be updated.
#   INCDIR is the writeable destination directory for the merged output.
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
	UPDIR <- args[1]
	DBDIR <- args[2]
	INCDIR <- args[3]
}
#UPDIR <- "test_example/updates/2026-42-42"
#DBDIR <- "test_example/latest"
#INCDIR <- "test_example/incremental/2026-42-42"


# A vector of file suffixes, to make looping easier.
suffix_vec <- c(
	"cbatch", "concentration", 
	"ebatch", "extraction",
	"qbatch", "qblock", "qcontrol", "qassay", 
	"sample"
)


flagdups_vec <- c("concentration", "extraction", "sample")

duplicates_df <- data.frame(
  sample_id = character(),
  duplicated_in = character()
)

# Read in each matching pair of tables from UPDIR and DBDIR.
for (i in 1:length(suffix_vec)) {
  
  latest_df <- as.data.frame(
    read.table(
      paste0(DBDIR, "/watchdb.", suffix_vec[i] , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      fill=TRUE,
      comment.char="",
      check.names=FALSE)
  )
  update_df <- as.data.frame(
    read.table(
      paste0(UPDIR, "/update.", suffix_vec[i] , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      fill=TRUE,
      comment.char="",
      check.names=FALSE)
  )
  
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
  if (suffix_vec[i] %in% flagdups_vec) {
    dups_df <- merged_df %>% group_by(sample_id) %>% filter(n() > 1) %>% ungroup()
    if (nrow(dups_df) > 0) {
      dups_df <- dups_df %>% select(sample_id)
      dups_df$duplicated_in <- suffix_vec[i]
      duplicates_df <- bind_rows(duplicates_df, dups_df)
    }
  }
  
  # Write the final database table and the incremental table to INCDIR. The latter retains 
  # "source_table" for easier identification of potential issues with the update.
  if (nrow(merged_df) > 0) {
    dbfn <- paste0(INCDIR, "/watchdb.", suffix_vec[i], ".txt", sep="")
    write.table(merged_df %>% select(!source_table), file = dbfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

    incrfn <- paste0(INCDIR, "/merge.", suffix_vec[i], ".txt", sep="")
    write.table(merged_df, file = incrfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
  }
  
}

if (nrow(duplicates_df) > 0) {
  dupfn <- paste0(INCDIR, "/merge.DUPLICATE_SAMPLES.txt", sep="")
  write.table(duplicates_df, file = dupfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
}

# Need to send the duplicate count back to bash so the control script can do the right thing.
ret_val <- as.numeric(nrow(duplicates_df))
cat(ret_val)

