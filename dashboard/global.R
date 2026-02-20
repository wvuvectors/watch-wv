source("addins/baselib.R")
source("addins/basevar.R")
source("addins/basefun.R")

source("addins/version.R")
source("addins/dbsources.R")



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

Sys.setenv(TZ="America/New_York")
today <- Sys.Date()
#testd <- as.Date("2025-12-20")

# This might be a kludge, but when the epi week is 53 and the current month is January, 
# we have to use the previous year when calculating the start date of the current epi week.
this_epiweek <- lubridate::epiweek(today)
this_epiyear <- year(today)
if (this_epiweek == 53 & month(today) == 1) {
  this_epiyear <- this_epiyear - 1
}
this_week <- get_date_from_epi_week(this_epiyear, this_epiweek)



# Load input files:
# First, all the abundance data.
df_all <- as.data.frame(read.table(WVD_RSS_F, sep="\t", header=TRUE, check.names=FALSE))

# Load the resource tables, primarily to retrieve metadata for active locations.
# Also need the active facilities so we can report the response rate for each county.
# We simplify the county and facility id columns to unnamed vectors for certain applications.
resources <- excel2df(WVD_RESOURCE_F)
df_counties <- resources$county %>% filter(county_id %in% df_all$location_id)
COUNTIES <- unname(unique(df_counties$county_id))
df_facilities <- resources$location %>% filter(tolower(location_status) == "active" & tolower(location_category) == "wwtp")
FACILITIES <- unname(unique(df_facilities$location_id))

# We use all the abundance data to identify counties; however, we don't need the facility data.
# So here we simplify to just the county and state rows.
# This is the primary data frame used across the dashboard.
df_rss <- df_all %>% filter(location_id %in% COUNTIES | location_id == "WV")

# Summary table comes in with the year and epi week stored in separate columns.
# We extract the Sunday of that week to use as a primary date, for plotting, etc.
df_rss <- df_rss %>% mutate(primary_date = as.Date(lubridate::ymd(get_date_from_epi_week(epi_year, epi_week))))

# Identify any stale data in df_rss and flag in a new column called stale.
# This is set by the global var STALE_THRESHOLD_DAYS and compares the current and 
# row epi week (so its a bit crude). This is equivalent to the function isStale, 
# but that function is not vectorized to work within mutate. Yet.
df_rss <- df_rss %>% mutate(stale = ifelse(as.numeric(difftime(this_week, primary_date, units = "days")) > STALE_THRESHOLD_DAYS, TRUE, FALSE))
#df_rss <- df_rss %>% mutate(stale = ifelse(isStale(primary_date), TRUE, FALSE))

# To generate the maps, we need the county polygons which are available in an SPDF file.
# Only need the NAME and geometry columns for now.
county_spdf <- read_sf(WVD_COUNTY_F)

# Populate the map colors with appropriate values.
# map_county_spdf is the primary data frame used to build the map layers.
map_county_spdf <- county_spdf %>% 
  filter(NAME %in% COUNTIES) %>%
  select(county_id = NAME, geometry)
  
df_colors <- df_rss %>% 
	filter(location_id %in% COUNTIES) %>% 
	group_by(location_id, target) %>% 
	arrange(epi_date) %>% 
	summarize(
		this_epi_date = max(epi_date, na.rm = TRUE),
		stale = last(stale), 
		this_alevel = ifelse(stale == TRUE, 1, last(abundance_level)), 
		this_tlevel = ifelse(stale == TRUE, 1, last(trend_level)), 
		this_acolor = ifelse(stale == TRUE, ALEVEL_COLORS[1], ALEVEL_COLORS[last(abundance_level)]), 
		this_tcolor = ifelse(stale == TRUE, TLEVEL_COLORS[1], TLEVEL_COLORS[last(trend_level)]))
		

map_county_spdf <- left_join(map_county_spdf, df_colors, by = c("county_id" = "location_id"))



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


