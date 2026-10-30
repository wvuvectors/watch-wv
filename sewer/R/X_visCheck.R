#! /usr/bin/env Rscript

source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

# Accepts 1 argument on STDIN:
#   RFILE is the file that contains the result table.
#

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) > 0) {
#	RFILE <- args[1]
#}
RFILE <- "../data/latest/watchdb.result.txt"

# Load data files.
result_df <- as.data.frame(read.table(RFILE, sep="\t", header=TRUE, check.names=FALSE))

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
target_df <- resource_list[["target"]]


#
# Dates, dates, dates.
#

# Convert date strings into Date objects.
result_df$collection_start_datetime <- mdy_hms(result_df$collection_start_datetime)
result_df$collection_end_datetime <- mdy_hms(result_df$collection_end_datetime)

# Add epi week info.
result_df$year <- lubridate::year(result_df$collection_end_datetime)
result_df$month <- lubridate::month(result_df$collection_end_datetime)
result_df$epi_week <- lubridate::epiweek(result_df$collection_end_datetime)
result_df <- result_df %>%
  mutate(epi_year = case_when(
    result_df$month == 12 & result_df$epi_week == 1 ~ result_df$year + 1,
    TRUE ~ result_df$year
  ))

result_df <- result_df %>% 
	mutate(epi_date = as.integer(epi_year) + as.integer(epi_week)/100)

#
# Data cleaning.
#

now <- Sys.time()
this_year <- lubridate::year(now)
this_month <- lubridate::month(now)
this_epi_week <- lubridate::epiweek(now)
this_epi_year = case_when(
	this_month == 12 & this_epi_week == 1 ~ this_year + 1,
	TRUE ~ this_year
)
this_epi_date = as.integer(this_epi_year) + as.integer(this_epi_week)/100


# Restrict results to those that pass sample QC.
cleaned_df <- result_df %>% 
	filter(
	  tolower(sample_qc) == "pass", 
  	this_epi_date - epi_date <= 0.8,
	  !is.na(target_copies_flow_normalized_per_capita)
	)


for (i in 1:length(target_df$target_id)) {
	targ <- target_df$target_id[i]
	plot_df <- cleaned_df %>% filter(target == targ)
	
	if (nrow(plot_df) > 0) {
		p1 <- ggplot(plot_df, aes()) + 
			geom_point(aes(x=as.factor(epi_date), y=target_copies_flow_normalized_per_capita, color=as.factor(location_id), fill=as.factor(location_id))) + 
			#scale_x_continuous(limits=c(-5, 5)) + 
			labs(y = "Target copies/mL", title=paste0(targ, sep=""))

		ggsave(paste0(str_replace_all(targ, "\\/", "_"), ".png", sep=""), plot = p1, width = 3306, height = 1892, units = c("px"), path = paste0("analyses/", this_epi_date, sep=""))
	}
}



