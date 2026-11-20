#! /usr/bin/env Rscript --vanilla

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 4 argument on STDIN:
#   INDIR is the directory that contains the latest watch data tables.
#		OUTDIR is the directory that will contain the output of this script.
#		MUF is the full path to the Marshall U data file, a csv in NWSS format.
#		DTSTAMP is a unique id to tag the archive file name (usually a datetime stamp).
#
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) > 0) {
#	INDIR <- args[1]
#	OUTDIR <- args[2]
#	MUF <- args[3]
# DTSTAMP <- args[4]
#}
INDIR <- "../data/latest"
OUTDIR <- "../data/xtend/nwss"
MUF <- "/Users/tpd0001/Library/CloudStorage/GoogleDrive-wvuvectors@gmail.com/My Drive/WaTCH-WV/WaTCH-WV SHARED/DATA FOR NWSS/READY/mu_nwss.LATEST.csv"

if (!exists("DTSTAMP")) {
	DTSTAMP <- format(Sys.time(), "%Y%m%d%H%M")
}


# Read in the MU data file. This should be a csv file in NWSS format already.
# Eventually we'll merge this with the watch data; to start, we'll keep it separate.
#
mu_df <- fread(MUF, quote="\"")


# Load in the watch data tables.
#
watch_tables <- c(
  "cbatch", "concentration", 
  "ebatch", "extraction",
  "qassay", "qbatch", "qblock", "qcontrol", 
  "sample", "result"
)

watch_list <- list()
for (dname in watch_tables) {
  d_df <- as.data.frame(
    read.table(
      paste0(INDIR, "/watchdb." , dname, ".txt", sep=""), 
      quote="", 
      sep="\t", 
      header=TRUE, 
      fill=TRUE,
      comment.char="", 
      check.names=FALSE)
  )
  d_df <- d_df %>% 
    mutate(across(everything(), as.character)) %>% 
    mutate(across(where(is.character), tolower))
  watch_list[[dname]] <- d_df
}
result_df <- watch_list[["result"]]


# We need info from the resources Excel file.
#
resource_tables <- c("location", "wwtp", "county", "target", "lab", "loci", "methods", "misc_nwss")
platef_in <- excel2df("resources/watchdb.all_tables.xlsx")
resource_list <- list()
for (tname in resource_tables) {
	r_df <- as.data.frame(platef_in[[tname]])
	colnames(r_df) <- as.character(r_df[1, ])
	r_df <- r_df[-1, ]
	r_df <- r_df %>% 
	  mutate(across(everything(), as.character)) %>% 
	  mutate(across(where(is.character), tolower))
	resource_list[[tname]] <- r_df
}
target_df <- resource_list[["target"]]
location_df <- resource_list[["location"]] %>% filter(location_status == "active")
location_df <- left_join(location_df, resource_list[["wwtp"]], by = join_by(location_primary_wwtp_id == wwtp_id))
location_df <- left_join(location_df, resource_list[["county"]], by = join_by(location_counties_served == county_name))

concentration_df <- left_join(watch_list[["concentration"]], watch_list[["cbatch"]], by="concentration_batch_id")
concentration_df <- left_join(concentration_df, resource_list[["misc_nwss"]], by = join_by(concentration_batch_matrix == method_watch_value))

extraction_df <- left_join(watch_list[["extraction"]], watch_list[["ebatch"]], by="extraction_batch_id")
extraction_df <- left_join(extraction_df, resource_list[["misc_nwss"]], by = join_by(extraction_batch_method == method_watch_value))

qcontrol_min_df <- watch_list[["qcontrol"]]
qcontrol_min_df <- qcontrol_min_df %>% 
  filter(assay_type == "pcr_negative_control") %>% 
  select(
    quantification_batch_id, 
    assay_target_name, 
    control_positive_partitions = assay_positive_partitions, 
    control_target_copies_per_ul_reaction = assay_target_copies_per_ul_reaction)

quant_df <- left_join(watch_list[["qassay"]], watch_list[["qbatch"]], by="quantification_batch_id")
quant_df <- left_join(quant_df, qcontrol_min_df, by=c("quantification_batch_id", "assay_target_name"))
quant_df <- left_join(quant_df, resource_list[["misc_nwss"]], by = join_by(quantification_batch_quantification_type == method_watch_value))

# Restrict our results to just those targets that are found in our target dataframe.
# This gets rid of BRSV, human RP, etc.
# Also restrict to active locations only.
#
result_df <- result_df %>% 
  filter(target %in% target_df$target_id) %>% 
  filter(location_id %in% location_df$location_id)


# Init a new dataframe with the NWSS fields. This will hold the WaTCH data.
# Start with the fields pulled directly from the result table, and constant values.
#
sewer_df <- data.frame(
  flow_rate = result_df$sample_flow, 
  lab_id = result_df$lab_id, 
  pcr_target_avg_conc = as.numeric(result_df$target_copies_per_ml), 
	population_served = result_df$location_population_served, 
  sample_collect_date = mdy_hm(result_df$collection_end_datetime), 
  sample_collect_time = mdy_hm(result_df$collection_end_datetime), 
  sample_id = result_df$sample_id, 
  result_id = result_df$result_id, 
  TMP_target = result_df$target, 
  TMP_location = result_df$location_id, 
  inhibition_detect = "not tested", 
  inhibition_method = "none", 
  lod_ref = "https://doi.org/10.1039/D0EW00946F", 
  lod_sewage = "2855", 
  num_no_target_control = "1", 
  pasteurized = "no", 
  pcr_target_units = "copies/L wastewater", 
  pretreatment = "no", 
  quant_stan_type = "dna", 
  rec_eff_percent = "-1", 
  reporting_jurisdiction = "WV", 
  stan_ref = "https://doi.org/10.1038/s41598-020-80715-1", 
  major_lab_method = "1", 
  major_lab_method_desc = "Nanotrap concentration of raw influent and quantification by ddPCR",
  inhibition_adjust = "", 
  other_jurisdiction = "", 
  rec_eff_spike_conc = "", 
  rec_eff_spike_matrix = "", 
  rec_eff_target_name = "", 
  dashboard_ignore = "no"
)
sewer_df$sample_collect_date = format(sewer_df$sample_collect_date, "%m/%d/%Y")
sewer_df$sample_collect_time = format(sewer_df$sample_collect_time, "%H:%M")

# Add target data to output df.
#
sewer_df <- sewer_df %>% 
  left_join(
    target_df %>% select(
      TMP_target = target_id, 
      pcr_gene_target = nwss_gene_target, 
      pcr_target = nwss_target, 
      pcr_gene_target_ref = nwss_gene_target_ref
    ), by="TMP_target")

# Add location data to output df.
#
sewer_df <- sewer_df %>% 
  left_join(
    location_df %>% select(
      TMP_location = location_id, 
      sample_matrix = location_collection_matrix, 
      sample_type = location_collection_type, 
      institution_type = nwss_category, 
      sample_location = location_type, 
      site_id = nwss_site_id, 
      wwtp_name = location_primary_wwtp_id, 
      county_names = county_fips, 
      zipcode = location_zipcode, 
      capacity_mgd = wwtp_capacity_mgd, 
      epaid = wwtp_epaid_id, 
      sample_location_specify = location_common_name, 
      wwtp_jurisdiction = location_primary_wwtp_id
    ), by="TMP_location")

sewer_df <- sewer_df %>% mutate(
  sample_location_specify = if_else(tolower(sample_location) == "wwtp", "", sample_location_specify),
  wwtp_jurisdiction = if_else(tolower(sample_location) == "wwtp", "wv", wwtp_jurisdiction)
)

# Add concentration data to output df.
#
sewer_df <- sewer_df %>% 
  left_join(
    concentration_df %>% select(
      sample_id, 
      concentration_method = method_nwss_value
    ), 
    by="sample_id", relationship = "many-to-many")

# Add extraction data to output df.
#
sewer_df <- sewer_df %>% 
  left_join(
    extraction_df %>% select(
      sample_id, 
      extraction_method = method_nwss_value
    ), 
    by="sample_id", relationship = "many-to-many")


# Add quantification data to output df.
#
sewer_df <- sewer_df %>% 
  left_join(
    quant_df %>% select(
      result_id = assay_id, 
      test_result_date = quantification_batch_date, 
      pcr_type = method_nwss_value, 
      control_positive_partitions, 
      control_target_copies_per_ul_reaction
    ), 
    by="result_id", relationship = "many-to-many")


# Remove duplicate rows. Stil not quite sure how this is happening, but it is.
#
sewer_df <- sewer_df %>% distinct()


sewer_df <- sewer_df %>% 
  mutate(ntc_amplify = if_else(
      as.numeric(control_positive_partitions) > 3 | as.numeric(control_target_copies_per_ul_reaction) > 1, 
      "yes", 
      "no", 
      missing = "no"), 
    sample_id = result_id
  )


sewer_df <- sewer_df %>% select(!c(control_positive_partitions, control_target_copies_per_ul_reaction, result_id, TMP_location, TMP_target))

# rbind SEWER and MU tables.
#
nwss_df <- rbind(sewer_df, mu_df)

# Replace NA values with empty string.
#
nwss_df <- nwss_df %>% mutate(across(everything(), ~ replace_na(as.character(.), "")))


# Remove rows with missing required fields.
#
nwss_df <- nwss_df %>% 
  filter(sample_collect_date != "") %>%
  filter(flow_rate != "")
  


# Remove targets not to send to NWSS yet.
#
nwss_df <- nwss_df %>% filter(pcr_target != "mev_wt")


# Write output table to files.
#
mfn <- paste0(OUTDIR, "/merged_nwss.LATEST.csv", sep="")
#write.table(nwss_df, file = mfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.csv(nwss_df, file = mfn, row.names = FALSE)

afn <- paste0(OUTDIR, "/archives/merged_nwss.", DTSTAMP, ".csv", sep="")
#write.table(nwss_df, file = afn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.csv(nwss_df, file = afn, row.names = FALSE)

sfn <- paste0(OUTDIR, "/archives/sewer_nwss.", DTSTAMP, ".csv", sep="")
#write.table(sewer_df, file = sfn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.csv(sewer_df, file = sfn, row.names = FALSE)

mufn <- paste0(OUTDIR, "/archives/mu_nwss.", DTSTAMP, ".csv", sep="")
#write.table(mu_df, file = mufn, sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
write.csv(mu_df, file = mufn, row.names = FALSE)


# p1 <- ggplot(df_all %>% filter(target == TARGETS[1]), aes()) + 
#   geom_histogram(aes(x=z, y=ifelse(after_stat(count) > 0, after_stat(count), NA)), color="darkgreen", fill="darkgray", binwidth=0.1) + 
#   scale_x_continuous(limits=c(-5, 5)) + 
#   labs(y="count", title=paste0(TARGETS[1], " at binwidth 0.1", sep="")) 
# ggsave(paste0("z_", TARGETS[1], ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = "analyses/")



ret_val <- nrow(nwss_df)
cat(ret_val)

