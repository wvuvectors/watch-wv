
source("addins/baselib.R")

excel2df <- function(fname) { 
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
  
  # assigning names to data frames
  names(data_frame) <- sheets
  
  # return the data frame
  data_frame
} 

TEST_F <- "../patchr/data/updates/2026-02-13_16-38/batches/CB211.xlsx"

# Load the resource tables, primarily to retrieve metadata for active locations.
# Also need the active facilities so we can report the response rate for each county.
# We simplify the county and facility id columns to unnamed vectors for certain applications.
resources <- excel2df(TEST_F)
#df_samples <- resources$county %>% filter(county_id %in% df_all$location_id)


