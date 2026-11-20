
source("addins/sewer_version.R")
source("addins/sewer_libs.R")
source("addins/base_vars.R")
source("addins/base_functions.R")
source("addins/sewer_sources.R")

INDIR <- "../data/latest"

# Load in the sewer result table.
#
this_df <- as.data.frame(
  read.table(
    file.path(INDIR, "watchdb.sample.txt"), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)
this_df <- this_df %>% mutate(across(everything(), as.character))

# Convert date strings into Date objects.
#
#this_df$collection_start_datetime <- mdy_hm(this_df$collection_start_datetime)
#this_df$collection_end_datetime <- mdy_hm(this_df$collection_end_datetime)



newd_df <- as.data.frame(
  read.table(
    file.path(INDIR, "new_dates.txt"), 
    quote="", 
    sep="\t", 
    header=TRUE, 
    fill=TRUE,
    comment.char="", 
    check.names=FALSE)
)
newd_df <- newd_df %>% mutate(across(everything(), as.character))

# Convert date strings into Date objects.
#
#newd_df$collection_start_datetime <- mdy_hm(newd_df$collection_start_datetime)
#newd_df$collection_end_datetime <- mdy_hm(newd_df$collection_end_datetime)


merge_df <- this_df %>%
  rows_update(newd_df, by = "sample_id")


merge_df$coll_start <- mdy_hm(merge_df$sample_collection_start_datetime)
merge_df$coll_end <- mdy_hm(merge_df$sample_collection_end_datetime)

bad_dates_df <- merge_df %>% filter(difftime(coll_end, coll_start, units = "hours") < 3) %>%
    distinct(sample_id, .keep_all = TRUE)

write.table(bad_dates_df, file = "bad_dates.r2.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)

merge_df <- merge_df %>% select(!c(coll_start, coll_end))
#merge_df$sample_collection_start_datetime <- format(merge_df$sample_collection_start_datetime, "%m/%d/%Y %H:%M")
#merge_df$sample_collection_end_datetime <- format(merge_df$sample_collection_end_datetime, "%m/%d/%Y %H:%M")


write.table(merge_df, file = "NEW_watchdb.sample.txt", sep = "\t", row.names = FALSE, quote = FALSE, append = FALSE)
