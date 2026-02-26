#! /usr/bin/env Rscript --vanilla

library(tidyverse)
library(dplyr)
library(data.table)
library(DT)
library(zoo)
library(rlang)
library(glue)

library(readxl)

library(scales)
library(lubridate)


excel2df <- function(fname) { 
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # return the data frame
  data_frame
} 

#The update (output), database, and input dirs are passed on STDIN. The sub-folders of the 
# input dir are hard-coded in the map_df dataframe.
f <- file("stdin")
open(f)
while(length(line <- readLines(f, n = 1)) > 0) {
 fpaths <- strsplit(line, " ")[[1]]
}
close(f)

#fpaths <- c("../patchr/tmp", "../patchr/data/latest", "/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/DRISCOLL_LAB/2 PROJECTS/WaTCH/TESTING_LAB/DATA_PCR")
UPDIR <- fpaths[1]
DBDIR <- fpaths[2]
INDIR <- fpaths[3]

# A map of file name labels to folder names.
map_df <- data.frame(
  batch_label = c("cbatch", "ebatch", "abatch"),
  in_folder = c("1 CB_PLATES", "2 EB_PLATES", "3 AB_PLATES"),
  batch_type = c("concentration", "extraction", "assay")
)

# A dataframe to hold the batch files that need to be processed in this run.
process_df <- data.frame(
  batch_type = character(),
  batch_id = character(),
  src_path = character(),
  file_name = character()
)


# Loop over the batch types stored in map_df.
for (i in 1:length(map_df$batch_label)) {
  
  # Retrieve the batches of this type that have already been processed. We don't 
  # want to run these again.
  batch_df <- as.data.frame(
    read.table(
      paste0(DBDIR, "/watchdb.", map_df$batch_label[i] , ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      check.names=FALSE)
  )
  # Assemble the key to extract the batch IDs from this file.
  keycol <- paste0(map_df$batch_type[i], "_batch_id", sep="")
  
  # List up the files for this batch type in the current input directory.
  fpath_list <- list.files(
    path = paste0(INDIR, "/", map_df$in_folder[i] ,"/", sep=""),
    pattern = "\\.xlsx$",
    full.names = TRUE
  )

  # Loop over the batch files in the input dir and flag batch files that need to
  # be processed. This is determined by the Batch ID entry in the Metadata sheet 
  # of the file.
  for (fpath in fpath_list) {

    platef_in <- excel2df(fpath)
	  metadata_df <- as.data.frame(t(platef_in$Metadata))
	  metadata_df <- rownames_to_column(metadata_df, var = "V0")
	  colnames(metadata_df) <- as.character(metadata_df[1, ])
	  metadata_df <- metadata_df[-c(1:1), ]

	  # Test this batch id against the existing batch IDs in batch_df. If it doesn't 
	  # exist already, add it to the list of files to be processed.
	  this_id <- metadata_df$`Batch ID`[1]
	  test_df <- batch_df %>% filter(.data[[keycol]] == this_id)
	  if (nrow(test_df) == 0) {
	    up_df <- data.frame(
	      batch_type = map_df$batch_type[i],
	      batch_id = this_id,
	      src_path = fpath,
			  file_name = basename(fpath)
	    )
      process_df <- rbind(process_df, up_df)
	  }
  }
}

# If any assay batches are found, need to add the corresponding result file.
# If there is no result file, remove the assay batch from this run.
assays_df <- process_df %>% filter(batch_type == "assay")
for (i in 1:nrow(assays_df)) {
	bid <- assays_df$batch_id[i]
  fpath_list <- list.files(
    path = paste0(INDIR, "/4 AB_RESULTS/", sep=""),
    pattern = paste0("^", bid, sep=""),
    full.names = TRUE
  )
#  print(fpath_list)
  if (length(fpath_list) == 0) {
    process_df <- process_df %>% filter(!(batch_id == bid))
  } else {
    for (fpath in fpath_list) {
  		up_df <- data.frame(
  			batch_type = "result",
  			batch_id = bid,
  			src_path = fpath,
			  file_name = basename(fpath)
  		)
  		process_df <- rbind(process_df, up_df)
    }
  }
}

if (nrow(process_df) > 0) {
  outfn <- paste0(UPDIR, "/update.batch_files.txt", sep="")
  write.table(process_df, file = outfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
	
	for (i in 1:nrow(process_df)) {
		src_file <- process_df$src_path[i]
		fname <- process_df$file_name[i]
		file.copy(from = src_file, to = paste0(UPDIR, "/batches/", fname, sep=""))
	}
}

# Need to send the result count back to bash for further processing.
ret_val <- as.numeric(nrow(process_df))
cat(ret_val)



