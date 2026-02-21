#! /usr/bin/env Rscript --vanilla

source("addins/baselib.R")

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

# Loop through CPLATE files in this run.
# These are passed on STDIN as a space-delimited set of args.
# The first arg is always to destination folder.
#f <- file("stdin")
#open(f)
#while(length(line <- readLines(f, n = 1)) > 0) {
#  fpaths <- strsplit(line, " ")[[1]]
#}
#close(f)

fpaths <- c("../patchr/tmp/", "../patchr/tmp/batches/CB218.xlsx")
outdir <- fpaths[1]

main_df <- data.frame(
  concentration_id = character(),
  concentration_well = character(),
  concentration_batch_id = character(),
  concentration_comment = character(),
  concentration_spike = character()
)

batch_df <- data.frame(
  concentration_batch_id = character(),
  concentration_batch_type = character(),
  concentration_batch_record_version = character(),
  concentration_date = character(),
  concentration_run_by = character(),
  concentration_machine = character(),
  concentration_batch_comment = character(),
  concentration_input_ml = character(),
  concentration_output_ml = character(),
  concentration_method = character(),
  concentration_method_lot_id = character()
)

for (i in 2:length(fpaths)) {
  this_fpath <- fpaths[i]

	platef_in <- excel2df(this_fpath)
	plate_df <- platef_in$Plate_Map %>% column_to_rownames(var = "...1")
	plate_df <- plate_df %>% select(where(~!all(is.na(.))))
	
	spike_df <- platef_in$Spike_Map %>% column_to_rownames(var = "...1")
	comment_df <- platef_in$Comment_Map %>% column_to_rownames(var = "...1")

	metadata_df <- as.data.frame(t(platef_in$Metadata))
	metadata_df <- rownames_to_column(metadata_df, var = "V0")
	colnames(metadata_df) <- as.character(metadata_df[1, ])
	metadata_df <- metadata_df[-c(1:1), ]
	metadata_df <- metadata_df %>% select(where(~!all(is.na(.))))
	
	batch_keys <- colnames(metadata_df)
	batch_keys <- str_replace_all(batch_keys, " ", "_")
	batch_keys <- str_replace_all(batch_keys, "^", "concentration_")
	batch_keys <- tolower(batch_keys)
	colnames(metadata_df) <- batch_keys

	batch_df <- rbind(batch_df, metadata_df)

		# Init the main update table.
	mainup_df <- data.frame(
	  concentration_id = character(),
	  concentration_well = character(),
	  concentration_batch_id = character(),
	  concentration_comment = character(),
	  concentration_spike = character()
	)
	
#	batchup_df <- as.data.frame(
#	  setNames(
#	    replicate(length(batch_keys), character(0), simplify = FALSE),
#	    batch_keys
#	  )
#	)

	
	batch_id <- metadata_df$concentration_batch_id[1]
	for (i in 1:nrow(plate_df)) {
		well_r <- rownames(plate_df)[i]
		for (j in 1:ncol(plate_df)) {
			well_c <- colnames(plate_df)[j]
			if (as.numeric(well_c) < 10) {
				well_c <- paste0(0, well_c, sep="")
			}
			this_id <- plate_df[i, j]
			
			mainup_df <- add_row(
			  mainup_df,
			  concentration_id = this_id,
			  concentration_well = paste0(well_r,well_c,sep=""),
				concentration_batch_id = batch_id,
				concentration_comment = comment_df[i, j],
				concentration_spike = spike_df[i, j]
			)
		}
	}
	
	main_df <- rbind(main_df, mainup_df)
}

dfn <- paste0(outdir, "update.concentration.txt", sep="")
write.table(main_df, file = dfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

bfn <- paste0(outdir, "update.cbatch.txt", sep="")
write.table(batch_df, file = bfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)


