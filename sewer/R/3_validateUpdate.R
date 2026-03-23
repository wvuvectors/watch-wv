#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")


# Validate the update in UPDIR against the existing tables in DBDIR.
f <- file("stdin")
open(f)
while(length(line <- readLines(f, n = 1)) > 0) {
 fpaths <- strsplit(line, " ")[[1]]
}
close(f)

#fpaths <- c("../patchr/tmp", "../patchr/data/latest")
UPDIR <- fpaths[1]
DBDIR <- fpaths[2]

# A vector of file suffixes, to make looping easier.
suffix_vec <- c(
	"cbatch", "concentration", 
	"ebatch", "extraction",
	"abatch", "assay", 
	"control", "sample"
)

# Read in each matching pair of tables from UPDIR and DBDIR.
for (i in 1:length(suffix_vec)) {
  
  latest_df <- as.data.frame(
    read.table(
      paste0(DBDIR, "/watchdb.", suffix_vec[i] , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      check.names=FALSE)
  )
  update_df <- as.data.frame(
    read.table(
      paste0(UPDIR, "/update.", suffix_vec[i] , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      check.names=FALSE)
  )
  
	# If there are IDs in the update that already exist in the database, pull them out to 
	# the duplicates df.
	#
	# CODE HERE!
	#
  
}

if (nrow(duplicates_df) > 0) {
  outfn <- paste0(UPDIR, "/update.DUPLICATES.txt", sep="")
  write.table(duplicates_df, file = outfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
}

# Need to send the duplicate count back to bash so the control script can do the right thing.
ret_val <- as.numeric(nrow(duplicates_df))
cat(ret_val)

