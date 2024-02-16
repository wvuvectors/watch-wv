source("addins/baselib.R")
source("addins/basevar.R")
source("addins/basefun.R")

source("addins/version.R")
source("addins/dbsources.R")


#rsconnect::deployApp('path/to/your/app')

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
#today <- as.Date("2022-07-12")

# Load shape files
county_spdf <- readOGR( 
  dsn= paste0("shapefiles/wv_counties/") , 
  layer="WV_Counties",
  verbose=FALSE
)


# Load data files
df_result <- as.data.frame(read.table(DB_RESULTS, sep="\t", header=TRUE, check.names=FALSE))
df_sample <- as.data.frame(read.table(DB_SAMPLES, sep="\t", header=TRUE, check.names=FALSE))
df_seqr <- as.data.frame(read.table(DB_SEQR, sep="\t", header=TRUE, check.names=FALSE))

resources <- excel2df(RES_ALL)
df_active_loc <- resources$location %>% filter(tolower(location_status) == "active")
df_active_wwtp <- resources$wwtp %>% filter(wwtp_id %in% df_active_loc$location_primary_wwtp_id)


# Add colors to county layer.
#
county_spdf$color_group <- case_when(
  county_spdf$NAME %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "#EAAA00", 
  county_spdf$NAME %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "#00B140", 
  .default = "#eeeeee")


# Restrict results to those that pass sample QC
#
df_result <- df_result %>% filter(tolower(sample_qc) == "pass")
#df_result <- df_result %>% filter(tolower(sample_qc) == "pass" & location_id %in% df_active_loc$location_id)


# Convert date strings into Date objects.
#
df_result$collection_start_datetime <- mdy_hm(df_result$collection_start_datetime)
df_result$collection_end_datetime <- mdy_hm(df_result$collection_end_datetime)

df_sample$sample_collection_start_datetime <- mdy_hm(df_sample$sample_collection_start_datetime)
df_sample$sample_collection_end_datetime <- mdy_hm(df_sample$sample_collection_end_datetime)
df_sample$sample_recovered_datetime <- mdy_hm(df_sample$sample_recovered_datetime)
df_sample$sample_received_date <- mdy(df_sample$sample_received_date)

df_seqr$sample_collection_start_datetime <- mdy_hm(df_seqr$sample_collection_start_datetime)
df_seqr$sample_collection_end_datetime <- mdy_hm(df_seqr$sample_collection_end_datetime)


# Create a date to use as the primary reference and make sure all tables use the same date.
#
df_result$date_primary <- as.Date(df_result$collection_end_datetime)
df_seqr$date_primary <- as.Date(df_seqr$sample_collection_end_datetime)


# Add some date objects, for use in summary stats.
#
df_result$mmr_year <- lubridate::year(df_result$date_primary)
df_result$mmr_week <- lubridate::week(df_result$date_primary)
df_result$week_starting <- floor_date(df_result$date_primary, "week", week_start = 1)
df_result$week_ending <- ceiling_date(df_result$date_primary, "week", week_start = 1)

df_seqr$mmr_year <- lubridate::year(df_seqr$date_primary)
df_seqr$mmr_week <- lubridate::week(df_seqr$date_primary)
df_seqr$week_starting <- floor_date(df_seqr$date_primary, "week", week_start = 1)
df_seqr$week_ending <- ceiling_date(df_seqr$date_primary, "week", week_start = 1)


# Assign overarching lineage groups to our seqr data.
#
df_seqr$percent <- as.numeric(df_seqr$variant_proportion) * 100

df_seqr <- df_seqr %>% mutate(
	color_group = case_when(
		str_detect(lineage_group, "^B\\.") ~ "B*",
		str_detect(lineage_group, "^BA") ~ "BA*",
		str_detect(lineage_group, "^BA\\.2\\.86") ~ "BA.2.86",
		str_detect(lineage_group, "^XBB") ~ "XBB*",
		str_detect(lineage_group, "^EG") ~ "EG*",
		str_detect(lineage_group, "^HV") ~ "HV*",
		str_detect(lineage_group, "^JN") ~ "JN*"
	)
)
df_seqr$color_group <- replace_na(df_seqr$color_group, "Other XBB")


#################################################

# Split the full data set for use in the different environments (tabs).
#

df_rs <- df_result %>% filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow))

df_rs_meta <- df_active_loc %>% filter(tolower(location_category) == "wwtp") %>% 
	select(location_common_name, location_group, location_counties_served, location_population_served, location_primary_lab) %>% 
	rename(Site = location_common_name, Group = location_group, County = location_counties_served, Pop = location_population_served, Lab = location_primary_lab)


# Assign a default date column for plotting data.
#
df_rs$date_to_plot <- df_rs$week_ending
df_seqr$date_to_plot <- df_seqr$week_ending



