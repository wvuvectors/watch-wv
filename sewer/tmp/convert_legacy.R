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

library(rstatix)
library(openxlsx)

excel2df <- function(fname) { 
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x, col_names = FALSE))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # return the data frame
  return(data_frame)
} 

abatch_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.abatch.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

assay_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.assay.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

qbatch_df <- abatch_df %>% select(
  quantification_batch_id = assay_batch_id,
  quantification_batch_date = assay_date,
  quantification_batch_machine = assay_machine,
  quantification_batch_amplification_method = assay_amplification_method,
  quantification_batch_amplification_method_lot_id = assay_amplification_method_lot_id,
  quantification_batch_quantification_type = assay_quantification_type,
  quantification_batch_record_version = assay_batch_record_version,
  quantification_batch_analysis_software_version = assay_analysis_software_version,
  quantification_batch_run_by = assay_run_by,
  quantification_batch_comment = assay_batch_comment
)
qbatch_df$batch_type <- "Quantification"
qbatch_df <- qbatch_df %>% mutate(
	run_id = quantification_batch_id
)


qblock_df <- abatch_df %>% select(
  assay_block_reaction_vol_ul = assay_reaction_ul, 
  assay_block_method = assay_method, 
  assay_block_method_lot_id = assay_method_lot_id, 
  quantification_batch_id = assay_batch_id, 
  assay_block_id = assay_method
)
qblock_df$assay_block_comment <- "Legacy import"
qblock_df$assay_block_samples_in_block <- -1

qblock_df <- qblock_df %>% mutate(
  assay_block_id = quantification_batch_id, 
  assay_block_sample_input_vol_ul = 5
)

#qblock_df <- left_join(
#  qblock_df, 
#  assay_df %>% select(quantification_batch_id = assay_batch_id, assay_block_sample_input_vol_ul = assay_input_ul),
#  by = join_by(quantification_batch_id)
#)


qassay_df <- assay_df %>% select(
  assay_id,
  assay_type, 
  assay_class,
  sample_id, 
  assay_well = assay_location_in_batch, 
  quantification_batch_id = assay_batch_id, 
  assay_target_name = assay_target, 
  assay_target_locus = assay_target_genetic_locus, 
  assay_target_fluorophore, 
  assay_target_copies_per_ul_reaction, 
  assay_accepted_partitions = assay_accepted_droplets
)
qassay_df$validation_key <- 1
qassay_df$assay_positive_partitions <- NA

qassay_df <- qassay_df %>% mutate(
  assay_target_name = paste0(assay_target_name, " (", assay_target_locus, ")", sep=""),
  assay_block_id = quantification_batch_id,
  assay_target_abbreviation = case_when(
    str_detect(assay_target_name, "Human") ~ "Human", 
    str_detect(assay_target_name, "Influenza Virus A") ~ "FluA", 
    str_detect(assay_target_name, "Influenza Virus B") ~ "FluB", 
    str_detect(assay_target_name, "N1") ~ "N1", 
    str_detect(assay_target_name, "N2") ~ "N2", 
    str_detect(assay_target_name, "SC2") ~ "SC2", 
    str_detect(assay_target_name, "Syncitial Virus, Human") ~ "RSV", 
    str_detect(assay_target_name, "Syncitial Virus, Bovine") ~ "BRSV", 
    str_detect(assay_target_name, "Norovirus GII") ~ "Noro II", 
    str_detect(assay_target_name, "Norovirus GI") ~ "Noro I", 
  )
)

qassay_df <- qassay_df %>% select(-assay_target_locus)

qcontrol_df <- qassay_df %>% filter(assay_class == "control")
qassay_df <- qassay_df %>% filter(assay_class != "control")

qcontrol_df <- qcontrol_df %>% mutate(
  assay_type = case_when(
    str_detect(assay_type, "PCR Positive") ~ "pcr_positive_control", 
    str_detect(assay_type, "PCR Negative") ~ "pcr_negative_control", 
    .default = "unknown_control"
  ),
  sample_id = case_when(
    str_detect(assay_type, "pcr_positive_control") ~ paste0("ppc-", as.numeric(Sys.time()), sep=""), 
    str_detect(assay_type, "PCR Negative") ~ paste0("pnc-", as.numeric(Sys.time()), sep=""), 
    .default = paste0("unc-", as.numeric(Sys.time()), sep="")
  )
)
qcontrol_df <- qcontrol_df %>% select(-assay_class)



############################################



cbatch_OLD_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.cbatch.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

conc_OLD_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.concentration.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

cbatch_df <- cbatch_OLD_df %>% select(
  concentration_batch_record_version,
  concentration_batch_machine = concentration_machine, 
  concentration_batch_matrix = concentration_method, 
  concentration_batch_matrix_lot_id = concentration_method_lot_id, 
  concentration_batch_sample_input_vol_ml = concentration_input_ml, 
  concentration_batch_elution_output_vol_ul = concentration_output_ml, 
  concentration_batch_date = concentration_date, 
  concentration_batch_run_by = concentration_run_by, 
  concentration_batch_id,
  concentration_batch_comment
)

cbatch_df <- cbatch_df %>% mutate(
  batch_type = "Concentration", 
  concentration_batch_matrix_vol_ul = 150, 
  concentration_batch_enhancement_reagent = "ER1", 
  concentration_batch_enhancement_reagent_lot_id = NA, 
  concentration_batch_enhancement_reagent_vol_ul = 100, 
  concentration_batch_spike_material = "Respiratory Syncytial Virus, Bovine (BRSV)", 
  concentration_batch_spike_vol_ul = 15, 
  run_id = concentration_batch_id,
  validation_key = -1, 
  concentration_batch_elution_output_vol_ul = concentration_batch_elution_output_vol_ul * 1000
)

conc_df <- conc_OLD_df %>% select(
  concentration_id, 
  sample_id, 
  concentration_well = concentration_location_in_batch, 
  concentration_batch_id
)
conc_df$validation_key <- -1



############################################



ebatch_OLD_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.ebatch.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

extraction_OLD_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.extraction.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

ebatch_df <- ebatch_OLD_df %>% select(
  extraction_batch_record_version, 
  extraction_batch_id, 
  extraction_batch_machine = extraction_machine, 
  extraction_batch_method = extraction_method, 
  extraction_batch_method_lot_id = extraction_method_lot_id, 
  extraction_batch_sample_input_vol_ul = extraction_input_ul, 
  extraction_batch_elution_output_vol_ul = extraction_output_ul, 
  extraction_batch_eluant_material = extraction_eluant, 
  extraction_batch_date = extraction_date, 
  extraction_batch_run_by = extraction_run_by, 
  extraction_batch_comment
)

ebatch_df <- ebatch_df %>% mutate(
  batch_type = "Extraction", 
  run_id = extraction_batch_id,
  validation_key = -1
)


extraction_df <- extraction_OLD_df %>% select(
  extraction_id, 
  extraction_well = extraction_location_in_batch, 
  extraction_batch_id
)

extraction_df <- extraction_df %>% mutate(
  sample_id = sub("^(.+?)\\..+$", "\\1", extraction_id), 
  validation_key = -1
)



############################################



sample_OLD_df <- as.data.frame(
  read.table(
    paste0("../data/v1_format/watchdb.sample.txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)

sample_df <- sample_OLD_df %>% select(
  sample_id,
  location_id,
  sample_collection_by,
  sample_event,
  sample_qc,
  sample_collection_start_datetime,
  sample_collection_end_datetime,
  sample_recovered_datetime,
  sample_flow,
  sample_received_by,
  sample_received_date,
  sample_comment
)

# We need info from the resources Excel file for the new sample table.
#
resource_tables <- c("location", "wwtp", "county", "target", "lab", "loci", "methods")
platef_in <- excel2df("resources/watchdb.all_tables.xlsx")
resource_list <- list()
for (tname in resource_tables) {
  r_df <- as.data.frame(platef_in[[tname]])
  colnames(r_df) <- as.character(r_df[1, ])
  r_df <- r_df[-1, ]
  resource_list[[tname]] <- r_df
}

sample_df <- left_join(
  sample_df, 
  resource_list[["location"]] %>% select(
    location_id, 
    location_type, 
    county_id=location_counties_served
  ), 
  by=join_by(location_id)
)

sample_df <- sample_df %>% mutate(
  validation_key = -1,
  sample_collection_method = "Composite"
)


write.table(cbatch_df, file = "../data/latest/watchdb.cbatch.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(conc_df, file = "../data/latest/watchdb.concentration.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(ebatch_df, file = "../data/latest/watchdb.ebatch.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(extraction_df, file = "../data/latest/watchdb.extraction.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(qassay_df, file = "../data/latest/watchdb.qassay.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(qbatch_df, file = "../data/latest/watchdb.qbatch.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(qblock_df, file = "../data/latest/watchdb.qblock.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(qcontrol_df, file = "../data/latest/watchdb.qcontrol.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.table(sample_df, file = "../data/latest/watchdb.sample.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)


