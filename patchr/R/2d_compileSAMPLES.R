#! /usr/bin/env Rscript --vanilla

library(tidyverse)
library(dplyr)
library(data.table)
library(DT)
library(zoo)
library(rlang)
library(glue)
library(scales)
library(lubridate)


# f <- file("stdin")
# open(f)
# while(length(line <- readLines(f, n = 1)) > 0) {
#  fpaths <- strsplit(line, " ")[[1]]
# }
# close(f)

fpaths <- c("../patchr/tmp/")
UPDIR <- fpaths[1]

input_df <- as.data.frame(
	read.csv(
		paste0(UPDIR, "/batches/Samples.csv", sep=""), 
		quote="", 
		header=TRUE, 
		check.names=FALSE)
)

input_df <- input_df %>% mutate(across(where(is.character), ~ na_if(.,"")))

cleaned_df <- input_df %>% filter(
  !is.na(`Sample Collection Start`) & 
  !is.na(`Sample Collection End`)
)

allsample_df <- cleaned_df %>% select(
	sample_id = `Asset Tag ID`, 
	location_id = `Location`, 
	sample_collection_method = `Sample Collection Method`, 
	sample_collection_by = `Sample Collection By`, 
	sample_event = `Event Type`,
	sample_qc = `Sample QC Check`,
	sample_collection_start_datetime = `Sample Collection Start`,
	sample_collection_end_datetime = `Sample Collection End`,
	sample_recovered_datetime = `Sample Retrieved Date/Time`,
	sample_flow = `Sample Flow (MGD)`,
	sample_received_by = `Sample Received By`,
	sample_received_date = `Sample Received Date`,
	sample_comment = `Comments`)

allsample_df$sample_event[is.na(allsample_df$sample_event)] <- "Routine Surveillance"

newsample_ids <- c()
id_files <- c("update.concentration.txt", "update.extraction.txt", "update.assay.txt")
for (fn in id_files) {
	fpath <- paste0(UPDIR, "/", fn, sep="")
	if (file.exists(fpath)) {
		id_df <- as.data.frame(
			read.table(fpath, sep = "\t", quote="", header=TRUE, check.names=FALSE)
		)
		newsample_ids <- c(newsample_ids, id_df$sample_id)
	}
}
newsample_ids <- unique(newsample_ids)
newsample_df <- allsample_df %>% filter(sample_id %in% newsample_ids)

if (nrow(newsample_df) > 0) {
	sfn <- paste0(UPDIR, "update.sample.txt", sep="")
	write.table(newsample_df, file = sfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
}

ret_val <- nrow(newsample_df)
cat(ret_val)

