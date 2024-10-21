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

county_spdf <- read_sf(paste0(RES_BASE, "/shapefiles/wv_counties/WV_Counties.shp"))

# Load data files
df_pcr_wvu <- as.data.frame(read.table(DB_RESULTS_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr_mu <- as.data.frame(read.table(DB_RESULTS_MU, sep="\t", header=TRUE, check.names=FALSE))
df_pcr <- rbind(df_pcr_wvu, df_pcr_mu)

df_s_wvu <- as.data.frame(read.table(DB_SAMPLES_WVU, sep="\t", header=TRUE, check.names=FALSE))
df_s_mu <- as.data.frame(read.table(DB_SAMPLES_MU, sep="\t", header=TRUE, check.names=FALSE))
df_sample <- rbind(df_s_wvu, df_s_mu)

df_seqr <- as.data.frame(read.table(DB_SEQR, sep="\t", header=TRUE, check.names=FALSE))

resources <- excel2df(RES_ALL)
df_active_loc <- resources$location %>% filter(tolower(location_status) == "active")
df_active_wwtp <- resources$wwtp %>% filter(wwtp_id %in% df_active_loc$location_primary_wwtp_id)
df_active_county <- resources$county %>% filter(county_id %in% df_active_loc$location_counties_served)
df_active_loci <- resources$loci %>% filter(target_id %in% TARGETS)



# Add colors to county and location dataframes.
#
# KLUDGE for now; working on a dynamic color palette approach.
#
county_spdf$colorby <- case_when(
  county_spdf$NAME %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "#EAAA00", 
  county_spdf$NAME %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "#00B140", 
  .default = "#EEEEEE")

df_active_county$colorby <- case_when(
  df_active_county$county_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "#EAAA00", 
  df_active_county$county_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "#00B140", 
  .default = "#EEEEEE")

df_active_loc$colorby <- case_when(
  tolower(df_active_loc$location_primary_lab) == "zoowvu" ~ "#EAAA00", 
  tolower(df_active_loc$location_primary_lab) == "muidsl" ~ "#00B140", 
  .default = "#CCCCCC")


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
df_pcr$mmr_year <- lubridate::year(df_pcr$date_primary)
df_pcr$mmr_week <- lubridate::week(df_pcr$date_primary)
df_pcr$week_starting <- floor_date(df_pcr$date_primary, "week", week_start = 1)
df_pcr$week_ending <- ceiling_date(df_pcr$date_primary, "week", week_start = 1)

df_seqr$mmr_year <- lubridate::year(df_seqr$date_primary)
df_seqr$mmr_week <- lubridate::week(df_seqr$date_primary)
df_seqr$week_starting <- floor_date(df_seqr$date_primary, "week", week_start = 1)
df_seqr$week_ending <- ceiling_date(df_seqr$date_primary, "week", week_start = 1)


# Assign overarching lineage groups to our seqr data.
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

# Split the full data set for use in the different environments (tabs).
#
df_pcr$target_copies_fn_per_cap <- df_pcr$target_copies_fn_per_cap/df_pcr$target_per_capita_basis
df_pcr$target_copies_per_ldcap <- df_pcr$target_copies_per_ldcap/df_pcr$target_per_capita_basis
df_pcr$target_per_capita_basis <- 1

df_rs <- df_pcr %>% filter(tolower(event_type) == "routine surveillance" & !is.na(sample_flow))

# Assign a default date column for plotting data.
#
df_rs$date_to_plot <- df_rs$week_starting
df_seqr$date_to_plot <- df_seqr$week_starting
#df_rs$date_to_plot <- df_rs$week_ending
#df_seqr$date_to_plot <- df_seqr$week_ending
#df_seqr$date_to_plot <- df_seqr$date_primary

# https://covid.cdc.gov/covid-data-tracker/#variant-summary


dflist_rs <- list()
for (i in 1:length(TARGETS)) {
#	dflist_rs[[i]] <- df_rs %>% filter(target == TARGETS[i] & target_genetic_locus == GENLOCI[i])
	dflist_rs[[i]] <- df_rs %>% filter(target == TARGETS[i])
}


df_rs_meta <- df_active_loc %>% filter(tolower(location_category) == "wwtp") %>% 
	select(location_common_name, location_group, location_counties_served, location_population_served, location_primary_lab) %>% 
	rename(Site = location_common_name, Group = location_group, County = location_counties_served, Pop = location_population_served, Lab = location_primary_lab)


# Pre-calculate freshness & deltas for each county.
# Add this info to df_active_county.
#
anames <- c("region_name", "region_geolevel", DISEASES)
df_regions <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, anames)), stringsAsFactors = FALSE)
df_fresh <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, anames)), stringsAsFactors = FALSE)

for (county in df_active_county$county_id) {
	#print(county)
	this_rowc <- list("region_name" = c(county), "region_geolevel" = c("county"))
	this_rowf <- list("region_name" = c(county), "region_geolevel" = c("county"))
	
	locations <- (df_active_loc %>% filter(tolower(location_category) == "wwtp" & location_counties_served == county))$location_id
	
	for (i in 1:length(TARGETS)) {
		this_disease <- DISEASES[[i]]
		df_this <- dflist_rs[[i]] %>% filter(location_id %in% locations)
		#print(df_this)
		delta <- calcDelta(df_this, VIEW_RANGE_PRIMARY)
		fresh <- calcFresh(df_this)
		this_rowc <- append(this_rowc, list(this_disease = c(delta)))
		this_rowf <- append(this_rowf, list(this_disease = c(fresh)))
	}
	df_regions <- rbind(df_regions, as.data.frame(this_rowc))
	df_fresh <- rbind(df_fresh, as.data.frame(this_rowf))

	for (loc_id in locations) {
		this_rowl <- list("region_name" = c(loc_id), "region_geolevel" = c("facility"))
		this_rowg <- list("region_name" = c(loc_id), "region_geolevel" = c("facility"))
		for (i in 1:length(TARGETS)) {
			this_disease <- DISEASES[[i]]
			df_this <- dflist_rs[[i]] %>% filter(location_id == loc_id)
			#print(df_this)
			delta <- calcDelta(df_this, VIEW_RANGE_PRIMARY)
			fresh <- calcFresh(df_this)
			this_rowl <- append(this_rowl, list(this_disease = c(delta)))
			this_rowg <- append(this_rowg, list(this_disease = c(fresh)))
		}
		df_regions <- rbind(df_regions, as.data.frame(this_rowl))
		df_fresh <- rbind(df_fresh, as.data.frame(this_rowg))
	}
	
}
colnames(df_regions) <- anames

for (colname in colnames(df_regions)) {
	if (colname != "region_name" & colname != "region_geolevel") {
		df_regions[[colname]] <- as.numeric(df_regions[[colname]])
		#df_regions[[colname]] <- replace_na(df_regions[[colname]], "-")
	}
}
df_regions$avg <- rowMeans(df_regions[,DISEASE_RS], na.rm = TRUE)
#df_regions <- df_regions %>% arrange(desc(avg))
df_regions <- df_regions %>% arrange(region_name)


df_regions$region_lab <- case_when(
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_id ~ "zoowvu", 
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_id ~ "muidsl", 
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "zoowvu"))$location_counties_served ~ "zoowvu", 
  df_regions$region_name %in% (df_active_loc %>% filter(tolower(location_category) == "wwtp" & tolower(location_primary_lab) == "muidsl"))$location_counties_served ~ "muidsl", 
  .default = "Other")

df_regions$region_lab <- replace_na(df_regions$region_lab, "Other")


colnames(df_fresh) <- anames

for (colname in colnames(df_fresh)) {
	if (colname != "region_name" & colname != "region_geolevel") {
		df_fresh[[colname]] <- as.numeric(df_fresh[[colname]])
	}
}


# getColorPal <- function(basis) {
# 	pal <- case_when(
# 		basis == "lab" ~ labPal()
# 	)
# 	return(pal)
# }
# 
# labPal <- colorFactor(
# 	palette = c("#EAAA00", "#00B140", "#EEEEEE"),
# 	domain = df_regions$lab
# )
# 
# alertPal <- colorFactor(
# 	palette = c("#EAAA00", "#00B140", "#EEEEEE"),
# 	domain = df_regions$lab
# )



# for (i in 1:length(TARGETS)) {
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
