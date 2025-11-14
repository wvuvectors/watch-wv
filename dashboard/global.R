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
#today <- as.Date("2022-07-12")

this_epiweek <- lubridate::epiweek(today)


county_spdf <- read_sf(paste0(RES_BASE, "/shapefiles/wv_counties/WV_Counties.shp"))

# Load data files
df_pcr_wvu <- as.data.frame(read.table(DB_RESULTS_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr_mu <- as.data.frame(read.table(DB_RESULTS_MU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr <- rbind(df_pcr_wvu, df_pcr_mu)

df_s_wvu <- as.data.frame(read.table(DB_SAMPLES_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_s_mu <- as.data.frame(read.table(DB_SAMPLES_MU, sep="\t", header=TRUE, check.names=FALSE))
df_sample <- rbind(df_s_wvu, df_s_mu)

df_seqr <- as.data.frame(read.table(DB_SEQR, sep="\t", header=TRUE, check.names=FALSE))
df_seqr <- df_seqr %>% rename(location_id = location)

df_alerts <- as.data.frame(read.table(DB_ALERTS, sep="\t", header=TRUE, check.names=FALSE))


resources <- excel2df(RES_ALL)
df_active_loc <- resources$location %>% filter(tolower(location_status) == "active")
df_active_wwtp <- resources$wwtp %>% filter(wwtp_id %in% df_active_loc$location_primary_wwtp_id)
df_active_county <- resources$county %>% filter(county_id %in% df_active_loc$location_counties_served)
df_active_loci <- resources$loci %>% filter(target_id %in% TARGETS)


df_active_loc$dotsize <- case_when(
  df_active_loc$location_population_served < 26000 ~ 3250, 
  .default = df_active_loc$location_population_served/8)
  
# Restrict results to those that pass sample QC and have NTC below threshold
#
df_pcr <- df_pcr %>% filter(location_id %in% df_active_loc$location_id)
df_pcr <- df_pcr %>% filter(tolower(sample_qc) == "pass")
df_pcr <- df_pcr %>% filter(tolower(target_result_validated) != "ntc above threshold")
#df_pcr <- df_pcr %>% filter(location_id != "CheatLakeWWTP-01")

df_seqr <- df_seqr %>% filter(tolower(sample_id) != "ntc")

# Convert date strings into Date objects.
#
df_pcr$collection_start_datetime <- mdy_hm(df_pcr$collection_start_datetime)
df_pcr$collection_end_datetime <- mdy_hm(df_pcr$collection_end_datetime)

df_sample$sample_collection_start_datetime <- mdy_hm(df_sample$sample_collection_start_datetime)
df_sample$sample_collection_end_datetime <- mdy_hm(df_sample$sample_collection_end_datetime)
df_sample$sample_recovered_datetime <- mdy_hm(df_sample$sample_recovered_datetime)
df_sample$sample_received_date <- mdy(df_sample$sample_received_date)

df_seqr$sample_collection_start_datetime <- mdy_hm(df_seqr$sample_collection_start_datetime)
df_seqr$sample_collection_end_datetime <- mdy_hm(df_seqr$sample_collection_end_datetime)


# Create a date to use as the primary reference and make sure all tables use the same date.
#
df_pcr$date_primary <- as.Date(df_pcr$collection_end_datetime)
df_seqr$date_primary <- as.Date(df_seqr$sample_collection_end_datetime)


# Add some date objects, for use in summary stats.
#

df_pcr$epi_week <- lubridate::epiweek(df_pcr$date_primary)

df_pcr$week_starting <- floor_date(df_pcr$date_primary, "week", week_start = 7)
df_pcr$week_ending <- df_pcr$week_starting+6
#df_pcr$week_ending <- ceiling_date(df_pcr$date_primary, "week", week_start = 7)

df_seqr$epi_week <- lubridate::epiweek(df_seqr$date_primary)

df_seqr$week_starting <- floor_date(df_seqr$date_primary, "week", week_start = 7)
df_seqr$week_ending <- df_seqr$week_starting+6
#df_seqr$week_ending <- ceiling_date(df_seqr$date_primary, "week", week_start = 7)


# Assign overarching lineage groups to our seqr data.
# https://covid.cdc.gov/covid-data-tracker/#variant-summary
#
#df_seqr$percent <- as.numeric(df_seqr$variant_proportion) * 100

df_seqr <- df_seqr %>% mutate(
	color_group = case_when(
		str_detect(variant, "^B\\.") ~ "Alpha-Kappa B*",
		str_detect(variant, "^B\\.1\\.1\\.529") ~ "Omicron B*",
		str_detect(variant, "^BA\\.2\\.86") ~ "Pirola BA.2.86",
		str_detect(variant, "^BA") ~ "Omicron BA*",
		str_detect(variant, "^XBB") ~ "Omicron XBB",
		str_detect(variant, "^EG") ~ "Eris EG*",
		str_detect(variant, "^KP\\.") ~ "Pirola KP*",
		str_detect(variant, "^LB\\.") ~ "Pirola LB*",
		str_detect(variant, "^JN\\.") ~ "Pirola JN*",
		str_detect(variant, "^XEC") ~ "Pirola XEC"
	)
)
df_seqr$color_group <- replace_na(df_seqr$color_group, "Other Omicron")



#################################################

df_pcr$target_copies_fn_per_cap <- df_pcr$target_copies_fn_per_cap/df_pcr$target_per_capita_basis
df_pcr$target_copies_per_ldcap <- df_pcr$target_copies_per_ldcap/df_pcr$target_per_capita_basis
df_pcr$target_per_capita_basis <- 1

df_rs <- df_pcr %>% filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow))

# Assign a default date column for plotting data.
#
df_rs$date_to_plot <- df_rs$week_ending
df_seqr$date_to_plot <- df_seqr$week_ending
this_week <- floor_date(today, "week", week_start = 7) + 6

# Use this block instead if you want to refer everything to the start of the week
#
#df_rs$date_to_plot <- df_rs$week_starting
#df_seqr$date_to_plot <- df_seqr$week_starting
#this_week <- floor_date(today, "week", week_start = 1)
#


# Split the full data set for use in the different environments (tabs).
#
dflist_rs <- list()
for (i in 1:length(TARGETS)) {
#	dflist_rs[[i]] <- df_rs %>% filter(target == TARGETS[i] & target_genetic_locus == GENLOCI[i])
	dflist_rs[[i]] <- df_rs %>% filter(target == TARGETS[i])
}

dflist_alerts <- list()
for (i in 1:length(TARGETS)) {
	dflist_alerts[[i]] <- df_alerts %>% filter(target == TARGETS[i])
	dflist_alerts[[i]]$region_geolevel <- case_when(
		dflist_alerts[[i]]$region_name %in% df_active_county$county_id ~ "county", 
		dflist_alerts[[i]]$region_name %in% df_active_loc$location_id ~ "facility", 
		.default = "state")
		
	dflist_alerts[[i]]$abundance_color <- case_when(
		dflist_alerts[[i]]$abundance_pct_change < ALERT_LEVEL_THRESHOLDS[1] ~ ALERT_LEVEL_COLORS[1], 
		dflist_alerts[[i]]$abundance_pct_change >= ALERT_LEVEL_THRESHOLDS[1] & dflist_alerts[[i]]$abundance_pct_change < ALERT_LEVEL_THRESHOLDS[2] ~ ALERT_LEVEL_COLORS[2], 
		dflist_alerts[[i]]$abundance_pct_change >= ALERT_LEVEL_THRESHOLDS[2] & dflist_alerts[[i]]$abundance_pct_change < ALERT_LEVEL_THRESHOLDS[3] ~ ALERT_LEVEL_COLORS[3], 
		dflist_alerts[[i]]$abundance_pct_change >= ALERT_LEVEL_THRESHOLDS[3] ~ ALERT_LEVEL_COLORS[4], 
		.default = ALERT_LEVEL_COLORS[5])

	dflist_alerts[[i]]$abundance_level <- case_when(
		dflist_alerts[[i]]$abundance_pct_change < ALERT_LEVEL_THRESHOLDS[1] ~ ALERT_LEVEL_STRINGS[1], 
		dflist_alerts[[i]]$abundance_pct_change >= ALERT_LEVEL_THRESHOLDS[1] & dflist_alerts[[i]]$abundance_pct_change < ALERT_LEVEL_THRESHOLDS[2] ~ ALERT_LEVEL_STRINGS[2], 
		dflist_alerts[[i]]$abundance_pct_change >= ALERT_LEVEL_THRESHOLDS[2] & dflist_alerts[[i]]$abundance_pct_change < ALERT_LEVEL_THRESHOLDS[3] ~ ALERT_LEVEL_STRINGS[3], 
		dflist_alerts[[i]]$abundance_pct_change >= ALERT_LEVEL_THRESHOLDS[3] ~ ALERT_LEVEL_STRINGS[4], 
		.default = ALERT_LEVEL_STRINGS[5])

	dflist_alerts[[i]]$trend_color <- case_when(
		dflist_alerts[[i]]$trend == TREND_STRINGS[1] ~ ALERT_LEVEL_COLORS[1], 
		dflist_alerts[[i]]$trend == TREND_STRINGS[2] ~ ALERT_LEVEL_COLORS[2], 
		dflist_alerts[[i]]$trend == TREND_STRINGS[3] ~ ALERT_LEVEL_COLORS[3], 
		dflist_alerts[[i]]$trend == TREND_STRINGS[4] ~ ALERT_LEVEL_COLORS[4], 
		dflist_alerts[[i]]$trend == TREND_STRINGS[5] ~ ALERT_LEVEL_COLORS[5], 
		dflist_alerts[[i]]$trend == TREND_STRINGS[6] ~ ALERT_LEVEL_COLORS[6], 
		dflist_alerts[[i]]$trend == TREND_STRINGS[7] ~ ALERT_LEVEL_COLORS[7], 
		.default = "#f3f3e1")

}


dflist_map_f <- list()
dflist_map_c <- list()
df_mappable_c <- county_spdf %>% 
	filter(NAME %in% df_active_loc$location_counties_served) %>%
	select(NAME,OBJECTID,AREA_,PERIMETER,DEP_24K_,DEP_24K_ID,STATE,FIPS,Shape_Leng,County,SHAPE_Le_1,SHAPE_Area,geometry)

for (i in 1:length(TARGETS)) {
	dflist_map_f[[i]] <- merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[i]], by.x="location_id", by.y="region_name", all.x = TRUE)

	dflist_map_f[[i]]$abundance_color <- replace_na(dflist_map_f[[i]]$abundance_color, ALERT_LEVEL_COLORS[5])
	dflist_map_f[[i]]$abundance_level <- replace_na(dflist_map_f[[i]]$abundance_level, ALERT_LEVEL_STRINGS[5])
	dflist_map_f[[i]]$trend_color <- replace_na(dflist_map_f[[i]]$trend_color, ALERT_LEVEL_COLORS[5])
	dflist_map_f[[i]]$trend <- replace_na(dflist_map_f[[i]]$trend, TREND_STRINGS[5])


	dflist_map_c[[i]] <- merge(df_mappable_c, dflist_alerts[[i]], by.x="NAME", by.y="region_name", all.x = TRUE)

	dflist_map_c[[i]]$abundance_color <- replace_na(dflist_map_c[[i]]$abundance_color, ALERT_LEVEL_COLORS[5])
	dflist_map_c[[i]]$abundance_level <- replace_na(dflist_map_c[[i]]$abundance_level, ALERT_LEVEL_STRINGS[5])
	dflist_map_c[[i]]$trend_color <- replace_na(dflist_map_c[[i]]$trend_color, ALERT_LEVEL_COLORS[5])
	dflist_map_c[[i]]$trend <- replace_na(dflist_map_c[[i]]$trend, TREND_STRINGS[5])
}


df_rs_meta <- df_active_loc %>% filter(tolower(location_category) == "wwtp") %>% 
	select(location_common_name, location_group, location_counties_served, location_population_served, location_primary_lab) %>% 
	rename(Site = location_common_name, Group = location_group, County = location_counties_served, Pop = location_population_served, Lab = location_primary_lab)


anames <- c("region_name", "region_geolevel")
df_regions <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, anames)), stringsAsFactors = FALSE)

for (county in df_active_county$county_id) {
	#print(county)
	this_rowc <- list("region_name" = c(county), "region_geolevel" = c("county"))
	#this_rowf <- list("region_name" = c(county), "region_geolevel" = c("county"))
	df_regions <- rbind(df_regions, as.data.frame(this_rowc))
	
	locations <- (df_active_loc %>% filter(tolower(location_category) == "wwtp" & location_counties_served == county))$location_id
	
	for (loc_id in locations) {
		this_rowl <- list("region_name" = c(loc_id), "region_geolevel" = c("facility"))
		df_regions <- rbind(df_regions, as.data.frame(this_rowl))
	}
	
}
colnames(df_regions) <- anames
df_regions <- df_regions %>% arrange(region_name)

df_regions$region_lab <- case_when(
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_id ~ "zoowvu", 
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_id ~ "muidsl", 
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "zoowvu", 
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "muidsl", 
  .default = "other")
df_regions$region_lab <- replace_na(df_regions$region_lab, "other")



# Pre-calculate risk, abundance, trend, and variant for latest week (all regions).
# Enables map coloring and informs the county table.
#
# This info appears in dflist_alerts.
#

# dflist_alerts <- list()
# for (i in 1:length(TARGETS)) {
# #	dflist_rs[[i]] <- df_rs %>% filter(target == TARGETS[i] & target_genetic_locus == GENLOCI[i])
# 	dflist_alerts[[i]] <- df_regions
# 	latest_date <- max(dflist_rs[[i]]$date_to_plot)
# 	df_latest <- dflist_rs[[i]] %>% 
# 		filter(date_to_plot >= latest_date-31) %>% 
# 		select("date_to_plot", "target_copies_fn_per_cap", "location_id")
# 	
# 	for (this_region in df_regions$region_name) {
# 		if ((df_regions %>% filter(region_name == this_region))$region_geolevel == "county") {
# 			# roll up county active facilities
# 			loc_ids <- (df_active_loc %>% filter(location_counties_served == this_region & location_category == "wwtp"))$location_id
# #		} else {
# 			# just the single facility, but need it as a vector
# #			loc_ids <- (df_active_loc %>% filter(location_id == this_region))$location_id
# #		}
# 		
# 		suppressMessages(
# 		status_df <- df_latest %>% 
# 								 filter(location_id %in% loc_ids) %>%
# 								 group_by(date_to_plot) %>% 
# 								 arrange(date_to_plot) %>%
# 								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE),
# 													 date_to_plot := date_to_plot)
# 		)
# 	
# 		vec_risk <- getRiskLevel(status_df)
# 		dflist_alerts[[i]]$risk_text <- vec_risk[1]
# 		dflist_alerts[[i]]$risk_color <- vec_risk[2]
# 
# 		vec_abundance <- getAbundance(status_df)
# 		dflist_alerts[[i]]$abundance_text <- vec_abundance[1]
# 		dflist_alerts[[i]]$abundance_color <- vec_abundance[2]
# 
# 		vec_trend <- getTrend(status_df)
# 		vec_variant <- getDominantVariant(status_df)
# 		dflist_alerts[[i]]$variant_text <- vec_variant[1]
# 		dflist_alerts[[i]]$variant_color <- vec_variant[2]
# 		} else {
# 		}
# 	}
# }

