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
  dsn= paste0(RES_BASE, "/shapefiles/wv_counties/") , 
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
df_active_county <- resources$county %>% filter(county_id %in% df_active_loc$location_counties_served)

# Add colors to county and location dataframes.
#
county_spdf$colorby <- case_when(
  county_spdf$NAME %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "#EAAA00", 
  county_spdf$NAME %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "#00B140", 
  .default = "#eeeeee")

df_active_county$colorby <- case_when(
  df_active_county$county_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "#EAAA00", 
  df_active_county$county_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "#00B140", 
  .default = "#eeeeee")


df_active_loc$colorby <- case_when(
  tolower(df_active_loc$location_primary_lab) == "zoowvu" ~ "#EAAA00", 
  tolower(df_active_loc$location_primary_lab) == "muidsl" ~ "#00B140", 
  .default = "#cccccc")


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
df_seqr$color_group <- replace_na(df_seqr$color_group, "Other")



#################################################

# Split the full data set for use in the different environments (tabs).
#
df_result$target_copies_fn_per_cap <- df_result$target_copies_fn_per_cap/df_result$target_per_capita_basis
df_result$target_copies_per_ldcap <- df_result$target_copies_per_ldcap/df_result$target_per_capita_basis
df_result$target_per_capita_basis <- 1

df_rs <- df_result %>% filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow))

# Assign a default date column for plotting data.
#
df_rs$date_to_plot <- df_rs$week_ending
df_seqr$date_to_plot <- df_seqr$week_ending

dflist_rs <- list()
for (i in 1:length(TARGETS_RS)) {
	dflist_rs[[i]] <- df_rs %>% filter(target == TARGETS_RS[i] & target_genetic_locus == GENLOCI_RS[i])
}


df_rs_meta <- df_active_loc %>% filter(tolower(location_category) == "wwtp") %>% 
	select(location_common_name, location_group, location_counties_served, location_population_served, location_primary_lab) %>% 
	rename(Site = location_common_name, Group = location_group, County = location_counties_served, Pop = location_population_served, Lab = location_primary_lab)


# Pre-calculate freshness & deltas for each county.
# Add this info to df_active_county.
#

anames <- c("County", DISEASE_RS)
df_alert_county <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, anames)), stringsAsFactors = FALSE)

for (county in df_active_county$county_id) {
	#print(county)
	this_row <- list("County" = c(county))
	
	locations <- (df_active_loc %>% filter(location_counties_served == county))$location_id
	
	for (i in 1:length(TARGETS_RS)) {
		this_disease <- DISEASE_RS[[i]]
		df_this <- dflist_rs[[i]] %>% filter(location_id %in% locations)
		#print(df_this)
		delta <- calcDelta(df_this, 12)
		this_row <- append(this_row, list(this_disease = c(delta)))
	}

	df_alert_county <- rbind(df_alert_county, as.data.frame(this_row))
}
colnames(df_alert_county) <- anames

for (colname in colnames(df_alert_county)) {
	if (colname != "County") {
		df_alert_county[[colname]] <- as.numeric(df_alert_county[[colname]])
		#df_alert_county[[colname]] <- replace_na(df_alert_county[[colname]], "-")
	}
}

df_alert_county$avg <- rowMeans(df_alert_county[,DISEASE_RS], na.rm = TRUE)
df_alert_county <- df_alert_county %>% arrange(desc(avg))




# for (i in 1:length(TARGETS_RS)) {
# 	col_f <- paste0("freshness_", i, sep="")
# 	col_d <- paste0("delta_", i, sep="")
# 	#df_active_loc <- df_active_county %>% add_column(!!col_f)
# 	#df_active_loc <- df_active_county %>% add_column(!!col_d)
# 	for (loc in df_active_loc$location_id) {
# 		df_this <- df_rs %>% filter(location_id == loc & target == TARGETS_RS[i] & target_genetic_locus == GENLOCI_RS[i])
# 		if (length(df_this$date_to_plot) == 0) {
# 			freshness <- "NA"
# 			delta <- "NA"
# 		} else {
# 			last_date <- max(df_this$date_to_plot, na.rm = TRUE)
# 			freshness <- ymd(today)-ymd(last_date)
# 			delta <- 100 * 
# 				mean((df_this %>% filter(ymd(date_to_plot) == ymd(last_date)))$target_copies_fn_per_cap, na.rm = TRUE) / 
# 				mean((df_this %>% filter(date_to_plot > (today %m-% months(12))))$target_copies_fn_per_cap, na.rm = TRUE)
# 		}
# 		df_active_loc[[col_f]] <- freshness
# 		df_active_loc[[col_d]] <- delta
# 	}
# }
# 
# 
