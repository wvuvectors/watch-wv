#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 argument on STDIN:
#   INCDIR is the directory that contains the incremental update files.
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
INCDIR <- args[1]
}
#INCDIR <- "../data/latest"
#RUNID <- basename(INCDIR)


# We need info from the resources Excel file.
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

# Read in the data tables.
#
data_tables <- c(
	"cbatch", "concentration", 
	"ebatch", "extraction",
	"qbatch", "qblock", "qcontrol", "qassay", 
	"sample"
)

data_list <- list()
for (dname in data_tables) {
	d_df <- as.data.frame(
	read.table(
    paste0(INCDIR, "/watchdb." , dname, ".txt", sep=""), 
    quote="", 
    sep="\t", 
    header=TRUE, 
		fill=TRUE,
		comment.char="", 
    check.names=FALSE)
)
	data_list[[dname]] <- d_df
}

cjoin_df <- left_join(
	data_list[["concentration"]],
	data_list[["cbatch"]],
	by=join_by(concentration_batch_id)
)

ejoin_df <- left_join(
	data_list[["extraction"]],
	data_list[["ebatch"]],
	by=join_by(extraction_batch_id)
)

qjoin_df <- left_join(
	data_list[["qassay"]],
	data_list[["qblock"]],
	by=join_by(assay_block_id, quantification_batch_id)
)


calc_df <- qjoin_df %>% select(
	assay_id, 
	sample_id, 
	target=assay_target_name, 
	assay_target_copies_per_ul_reaction, 
	assay_input_ul=assay_block_sample_input_vol_ul, 
	assay_reaction_ul=assay_block_reaction_vol_ul
)

calc_df <- left_join(
	calc_df, 
	data_list[["sample"]] %>% select(
		sample_id, 
		location_id, 
		event_type=sample_event, 
		collection_start_datetime=sample_collection_start_datetime, 
		collection_end_datetime=sample_collection_end_datetime, 
		sample_flow, 
		sample_qc
	), 
	by=join_by(sample_id)
)

calc_df$collection_start_datetime <- mdy_hm(calc_df$collection_start_datetime)
calc_df$collection_end_datetime <- mdy_hm(calc_df$collection_end_datetime)


calc_df <- left_join(
	calc_df, 
	cjoin_df %>% select(
		sample_id, 
		concentration_batch_matrix_vol_ul, 
		concentration_batch_enhancement_reagent_vol_ul, 
		conc_input_ml=concentration_batch_sample_input_vol_ml,
		conc_output_ul=concentration_batch_elution_output_vol_ul
	), 
	by=join_by(sample_id)
)

calc_df <- left_join(
	calc_df, 
	ejoin_df %>% select(
		sample_id, 
		extraction_input_ul=extraction_batch_sample_input_vol_ul,
		extraction_output_ul=extraction_batch_elution_output_vol_ul
	), 
	by=join_by(sample_id)
)

calc_df <- left_join(
	calc_df, 
	resource_list[["location"]] %>% select(
		location_id, 
		location_population_served,
		location_counties_served,
		lab_id=location_primary_lab,
		location_category
	), 
	by=join_by(location_id)
)


# Now do the calculations, in calc_df.
#
calc_df <- calc_df %>% mutate(
	collection_window_hours = as.numeric(difftime(as.POSIXct(collection_end_datetime), as.POSIXct(collection_start_datetime), units="hours")),
	across(c(assay_target_copies_per_ul_reaction, 
					 assay_reaction_ul, 
					 assay_input_ul, 
					 extraction_output_ul, 
					 extraction_input_ul, 
					 conc_output_ul, 
					 conc_input_ml, 
					 location_population_served, 
					 sample_flow), as.numeric)
)

# Copies of target per mL of wastewater is the basic result, so we calculate that first:
# 1000 uL/mL * copies/uL_rxn * (rxn_total_uL/rxn_input_uL) * (extraction_output_uL/extraction_input_uL) * (concentration_output_mL/concentration_input_mL)
#
calc_df <- calc_df %>% mutate(
	target_copies_per_ml = 1000 * assay_target_copies_per_ul_reaction * (assay_reaction_ul/assay_input_ul) * (extraction_output_ul/extraction_input_ul) * ((conc_output_ul/1000)/conc_input_ml), 
	target_copies_per_ml_per_day = target_copies_per_ml / (collection_window_hours/24)
)

# Normalize by daily flow, where flow is available, and by population, where population is available:
# copies/mL * 3785.41 mL/gallon * flow million_gallons/day * 10^6 gallons/million_gallons * collection_time_hrs * 1 day/24 hours
#
calc_df <- calc_df %>% mutate(
	target_copies_per_ml_per_day_per_capita = target_copies_per_ml_per_day / as.numeric(location_population_served),
	target_copies_flow_normalized = target_copies_per_ml * 3785.41 * sample_flow * 10^6 * (collection_window_hours/24),
	target_copies_flow_normalized_per_capita = target_copies_flow_normalized / as.numeric(location_population_served)
)


# Validate the results against required columns.
# 	1 = all data present and properly formed.
# 	0 = all required data present and properly formed, but at least one optional data field is missing or malformed.
# 	-1 = at least one required data value missing or malformed.
# We return all the results and let the downstream applications decide how to filter.
	
result_df <- calc_df %>% mutate(
	validation_key = case_when(
		is.na(sample_qc) | str_detect(sample_qc, regex("Pass", ignore_case=TRUE), negate=TRUE) ~ -1, 
		is.na(location_id) | location_id == "" ~ -1, 
		is.na(location_category) | location_category == "" ~ 0, 
		is.na(target_copies_per_ml) | as.numeric(target_copies_per_ml) < 0 ~ -1, 
		.default = 1
	)
)

result_df <- result_df %>% select(
	assay_id,
  sample_id,
	target,
	lab_id,
	location_id,
	location_category,
	location_population_served,
	location_counties_served, 
	event_type,
	collection_start_datetime, 
	collection_end_datetime,
	collection_window_hours, 
	sample_flow, 
	sample_qc, 
	target_copies_per_ml,
	target_copies_per_ml_per_day, 
	target_copies_per_ml_per_day_per_capita, 
	target_copies_flow_normalized, 
	target_copies_flow_normalized_per_capita
)

# Print the result df.
rfn <- paste0(INCDIR, "/watchdb.result.txt", sep="")
write.table(result_df, file = rfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)


ret_val <- nrow(result_df)
cat(ret_val)
