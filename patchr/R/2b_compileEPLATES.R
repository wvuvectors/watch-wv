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

#f <- file("stdin")
#open(f)
#while(length(line <- readLines(f, n = 1)) > 0) {
#  fpaths <- strsplit(line, " ")[[1]]
#}
#close(f)

fpaths <- c("../patchr/tmp/")
UPDIR <- fpaths[1]

input_df <- as.data.frame(
	read.table(
		paste0(UPDIR, "/update.batch_files.txt", sep=""), 
		quote="", 
		sep="\t", 
		header=TRUE, 
		check.names=FALSE)
)

update_df <- input_df %>% filter(batch_type == "extraction")

main_df <- data.frame(
  extraction_id = character(),
  sample_id = character(),
  extraction_well = character(),
  extraction_batch_id = character(),
  extraction_comment = character(),
  extraction_storage_location = character()
)

batch_df <- data.frame(
  extraction_batch_id = character(),
  extraction_batch_type = character(),
  extraction_batch_record_version = character(),
  extraction_date = character(),
  extraction_run_by = character(),
  extraction_machine = character(),
  extraction_batch_comment = character(),
  extraction_input_ul = character(),
  extraction_output_ul = character(),
  extraction_method = character(),
  extraction_method_lot_id = character(),
  extraction_eluant = character()
)


for (i in nrow(update_df)) {
  this_fn <- update_df$file_name[i]
  this_fpath <- paste0(UPDIR, "/batches/", this_fn, sep="")

	platef_in <- excel2df(this_fpath)
	plate_df <- platef_in$Plate_Map %>% column_to_rownames(var = "...1")
	plate_df <- plate_df %>% select(where(~!all(is.na(.))))
	
	storage_df <- platef_in$Storage_Map %>% column_to_rownames(var = "...1")
	storage_df <- storage_df %>% select(where(~!all(is.na(.))))
	comment_df <- platef_in$Comment_Map %>% column_to_rownames(var = "...1")
	comment_df <- comment_df %>% select(where(~!all(is.na(.))))
	
	metadata_df <- as.data.frame(t(platef_in$Metadata))
	metadata_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(metadata_df) <- as.character(metadata_df[1, ])
	metadata_df <- metadata_df[-c(1:1), ]
	metadata_df <- metadata_df %>% slice(1:1)
	metadata_df <- metadata_df %>% select(where(~!all(is.na(.))))
	
	batch_keys <- colnames(metadata_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- tolower(batch_keys)
	batch_keys <- str_replace_all(batch_keys, "^", "extraction_")
	batch_keys <- str_replace_all(batch_keys, "^extraction_extraction_", "extraction_")
	colnames(metadata_df) <- batch_keys

	metadata_df$extraction_date <- convertToDateTime(as.numeric(metadata_df$extraction_date))

	batch_df <- rbind(batch_df, metadata_df)

		# Init the main update table.
	mainup_df <- data.frame(
	  extraction_id = character(),
	  sample_id = character(),
	  extraction_well = character(),
	  extraction_batch_id = character(),
	  extraction_comment = character(),
	  extraction_storage_location = character()
	)
	
#	batchup_df <- as.data.frame(
#	  setNames(
#	    replicate(length(batch_keys), character(0), simplify = FALSE),
#	    batch_keys
#	  )
#	)

	batch_id <- metadata_df$extraction_batch_id[1]
	date_str <- str_replace_all(metadata_df$extraction_date[1], "-", "")

	for (i in 1:nrow(plate_df)) {
		well_r <- toupper(rownames(plate_df)[i])
		for (j in 1:ncol(plate_df)) {
			well_c <- colnames(plate_df)[j]
			if (as.numeric(well_c) < 10) {
				well_c <- paste0(0, well_c, sep="")
			}
			sample_id <- plate_df[i, j]
			if (!is.na(sample_id) & sample_id != "") {
  			mainup_df <- add_row(
  			  mainup_df,
  			  extraction_id = paste0(sample_id, ".", date_str, sep=""),
  			  sample_id = sample_id,
  			  extraction_well = paste0(well_r, well_c, sep=""),
  				extraction_batch_id = batch_id,
  				extraction_comment = comment_df[i, j],
  				extraction_storage_location = storage_df[i, j]
  			)
			}
		}
	}
	
	main_df <- rbind(main_df, mainup_df)
}

if (nrow(main_df) > 0) {
	dfn <- paste0(UPDIR, "update.extraction.txt", sep="")
	write.table(main_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

	bfn <- paste0(UPDIR, "update.ebatch.txt", sep="")
	write.table(batch_df, file = bfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
}

ret_val <- nrow(main_df)
cat(ret_val)

