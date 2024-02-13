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

# Convert this palette to match alert level!
mypalette <- colorBin(palette="plasma", domain=county_spdf@data$POPCH_PCT, na.color="transparent")


# Load data files
df_result <- as.data.frame(read.table(DB_RESULTS, sep="\t", header=TRUE, check.names=FALSE))
df_sample <- as.data.frame(read.table(DB_SAMPLES, sep="\t", header=TRUE, check.names=FALSE))

resources <- excel2df(RES_ALL)

# Clean up the input a bit:
# Restrict it to results from routine surveillance at WWTP sites.
# Also make sure they all have valid flow values.
#
df_wwtp <- resources$location %>% filter(location_category == "wwtp" & location_status == "active")
df_result <- df_result %>% filter(event_type == "Routine Surveillance" & 
                                    sample_qc == "Pass" & 
                                    location_id %in% df_wwtp$location_id & 
                                    !is.na(sample_flow))

#TARGET_CLASS = unique(resources$target %>% filter(target_id == TARGET_PRIMARY))$target_class
#DISEASE_PRIMARY = unique(resources$target %>% filter(target_id == TARGET_PRIMARY))$target_disease


#df_hosp1 <- as.data.frame(read.table("data/hospitalizations.csv", sep=",", header=TRUE, check.names=TRUE))
#colnames(df_hosp1)[1] <- "i"
#df_hosp1$time_value <- ymd(df_hosp1$time_value)
#df_hosp2 <- df_hosp1 %>% group_by(time_value, location_name) %>% summarize(rolling_weekly_mean = sum(value)*7)
#df_hosp2$mmr_year <- lubridate::year(df_hosp2$time_value)
#df_hosp2$mmr_week <- lubridate::week(df_hosp2$time_value)
#df_hospital <- df_hosp2 %>% group_by(mmr_year, mmr_week, location_name) %>% summarize(weekly_sum = round(sum(rolling_weekly_mean), digits = 0))


# Convert date strings into Date objects.
#
df_result$collection_start_datetime <- mdy_hm(df_result$collection_start_datetime)
df_result$collection_end_datetime <- mdy_hm(df_result$collection_end_datetime)
df_sample$sample_collection_start_datetime <- mdy_hm(df_sample$sample_collection_start_datetime)
df_sample$sample_collection_end_datetime <- mdy_hm(df_sample$sample_collection_end_datetime)
df_sample$sample_recovered_datetime <- mdy_hm(df_sample$sample_recovered_datetime)
df_sample$sample_received_date <- mdy(df_sample$sample_received_date)


# Create a date to use for the plot axis.
#
df_result$date_primary <- as.Date(df_result$collection_end_datetime)

# Add some date objects for use in summary stats.
#
df_result$mmr_year <- lubridate::year(df_result$date_primary)
df_result$mmr_week <- lubridate::week(df_result$date_primary)
df_result$week_starting <- floor_date(df_result$date_primary, "week", week_start = 1)
df_result$week_ending <- ceiling_date(df_result$date_primary, "week", week_start = 1)

df_result$date_to_plot <- df_result$week_ending

