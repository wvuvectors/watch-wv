source("addins/base_vars.R")
source("addins/base_functions.R")

source("addins/dash_version.R")
source("addins/dash_libs.R")
source("addins/dash_vars.R")
source("addins/dash_functions.R")
source("addins/dash_sources.R")


# color palette is ggthemes$calc
# 1 Chart 1  #004586	dark blue
# 2 Chart 2  #ff420e	red orange
# 3 Chart 3  #ffd320	yellow
# 4 Chart 4  #579d1c	green
# 5 Chart 5  #7e0021	brown
# 6 Chart 6  #83caff	light blue
# 7 Chart 7  #314004	dark green
# 8 Chart 8  #aecf00	light green
# 9 Chart 9  #4b1f6f	purple
#10 Chart 10 #ff950e	orange
#11 Chart 11 #c5000b	dark red
#12 Chart 12 #0084d1	Carolina blue


# Set some variables that are constant only for this particular instance.
# Other variables that are constant across all instances can be found in addins/dash_vars.
#
Sys.setenv(TZ="America/New_York")
TODAY <- Sys.Date()
#testd <- as.Date("2025-12-20")

# When the epi week is 53 and the current month is January, we have to use the previous 
# year when calculating the start date of the current epi week. Epi weeks begin on Sunday.
# The function get_date_from_epi_week can be found in addins/base_functions.
#
TODAY_EPIWEEK <- lubridate::epiweek(TODAY)
TODAY_EPIYEAR <- year(TODAY)
if (TODAY_EPIWEEK == 53 & month(TODAY) == 1) {
  TODAY_EPIYEAR <- TODAY_EPIYEAR - 1
}
TODAY_WEEKSTART <- get_date_from_epi_week(TODAY_EPIYEAR, TODAY_EPIWEEK)


# Load the data.
#
rollbyregion_df <- as.data.frame(
	read.table(
		REGION_DATAF, 
		quote="",
		sep="\t", 
		comment.char="",
		header=TRUE, 
		fill=TRUE,
		check.names=FALSE)
)


# Load the resource tables, primarily to retrieve metadata for active locations.
# Also need the active facilities so we can report the response rate for each county.
# We simplify the county and facility id columns to unnamed vectors for certain applications.
#
resource_tables <- c("target", "county")
platef_in <- excel2df(WATCH_RESOURCEF)
resource_list <- list()
for (tname in resource_tables) {
	r_df <- as.data.frame(platef_in[[tname]])
	colnames(r_df) <- as.character(r_df[1, ])
	r_df <- r_df[-1, ]
	r_df <- r_df %>% mutate(across(everything(), as.character))
	resource_list[[tname]] <- r_df
}
target_df <- resource_list[["target"]]
county_df <- resource_list[["county"]]

COUNTIES <- unname(unique(county_df$county_id))
#FACILITIES <- unname(unique(rollbyregion_df$location_id))

# Summary table comes in with the year and epi week stored in separate columns.
# We extract the Sunday of that week to use as a primary date, for plotting, etc.
#
rollbyregion_df <- rollbyregion_df %>% mutate(
	primary_date = as.Date(lubridate::ymd(get_date_from_epi_week(epi_year, epi_week)))
)

# Identify any stale data and flag in a new column called stale.
# This is set by the global var STALE_THRESHOLD_DAYS and compares the current and 
# row epi week (so its a bit crude). This is equivalent to the function isStale, 
# but that function is not vectorized to work within mutate. Yet.
#
rollbyregion_df <- rollbyregion_df %>% mutate(
	is_stale = ifelse(
		as.numeric(difftime(TODAY_WEEKSTART, primary_date, units = "days")) > STALE_THRESHOLD_DAYS, 
		TRUE, 
		FALSE)
)

# To generate the maps, we need the county polygons which are available in an SPDF file.
# Only need the NAME and geometry columns for now.
county_spdf <- read_sf(COUNTY_SHAPEF)

# Populate the map colors with appropriate values.
# map_county_spdf is the primary data frame used to build the map layers.
#
map_county_spdf <- county_spdf %>% 
  #filter(NAME %in% COUNTIES) %>%
  select(county_id = NAME, geometry)

colors_df <- rollbyregion_df %>% 
	filter(location_id %in% COUNTIES) %>% 
	group_by(location_id, target) %>% 
	arrange(epi_date) %>% 
	summarize(
		this_epi_date = max(epi_date, na.rm = TRUE),
		is_stale = last(is_stale), 
		this_alevel = ifelse(is_stale == TRUE, 1, last(abundance_level)), 
		this_tlevel = ifelse(is_stale == TRUE, 1, last(trend_level)), 
		this_acolor = ifelse(is_stale == TRUE, ALEVEL_COLORS[1], ALEVEL_COLORS[last(abundance_level)]), 
		this_tcolor = ifelse(is_stale == TRUE, TLEVEL_COLORS[1], TLEVEL_COLORS[last(trend_level)]))
		

map_county_spdf <- map_county_spdf %>% 
	left_join(colors_df, by = c("county_id" = "location_id"))



#################################################
#
# Processing for the seqr data. This is soon-to-be updated!
#
#df_seqr <- as.data.frame(read.table(WVD_SEQR_F, sep="\t", header=TRUE, check.names=FALSE))
#df_seqr <- df_seqr %>% rename(location_id = location)
#
# df_seqr$sample_collection_start_datetime <- mdy_hm(df_seqr$sample_collection_start_datetime)
# df_seqr$sample_collection_end_datetime <- mdy_hm(df_seqr$sample_collection_end_datetime)
# df_seqr$date_primary <- as.Date(df_seqr$sample_collection_end_datetime)
# df_seqr$epi_week <- lubridate::epiweek(df_seqr$date_primary)
# df_seqr$week_starting <- floor_date(df_seqr$date_primary, "week", week_start = 7)
# df_seqr$week_ending <- df_seqr$week_starting+6
#df_seqr$week_ending <- ceiling_date(df_seqr$date_primary, "week", week_start = 7)


# Assign overarching lineage groups to our seqr data.
# https://covid.cdc.gov/covid-data-tracker/#variant-summary
#
#df_seqr$percent <- as.numeric(df_seqr$variant_proportion) * 100

# df_seqr <- df_seqr %>% mutate(
# 	color_group = case_when(
# 		str_detect(variant, "^B\\.") ~ "Alpha-Kappa B*",
# 		str_detect(variant, "^B\\.1\\.1\\.529") ~ "Omicron B*",
# 		str_detect(variant, "^BA\\.2\\.86") ~ "Pirola BA.2.86",
# 		str_detect(variant, "^BA") ~ "Omicron BA*",
# 		str_detect(variant, "^XBB") ~ "Omicron XBB",
# 		str_detect(variant, "^EG") ~ "Eris EG*",
# 		str_detect(variant, "^KP\\.") ~ "Pirola KP*",
# 		str_detect(variant, "^LB\\.") ~ "Pirola LB*",
# 		str_detect(variant, "^JN\\.") ~ "Pirola JN*",
# 		str_detect(variant, "^XEC") ~ "Pirola XEC"
# 	)
# )
# df_seqr$color_group <- replace_na(df_seqr$color_group, "Other Omicron")
# 


