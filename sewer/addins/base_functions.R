
ci90 <- function(x) {
	0.5 * qt(0.80, length(x) - 1) * (sd(x) / sqrt(length(x)))
#	m <- mean(x)
#	se <- sd(x) / sqrt(x.length)
#	ci = qt(1 - (0.05 / 2), x.length - 1) * se
#	lower.ci = m - ci
#	upper.ci = m + ci
}

ci95 <- function(x) {
	0.5 * qt(0.95, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

ci98 <- function(x) {
	0.5 * qt(0.98, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

ci99 <- function(x) {
	0.5 * qt(0.99, length(x) - 1) * (sd(x) / sqrt(length(x)))
}

se <- function(x) {
	0.5 * sd(x) / sqrt(length(x))
}

get_date_from_epi_week <- function(year, epi_week) {
  # 1. Find the date of the first day of the year
  jan1 <- ymd(paste0(year, "-01-01", sep=""))

  # 2. Determine the day of the week for Jan 1st (Sunday = 1, Saturday = 7)
  # wday() with week_start= "Sunday" (default) is used
  jan1_wday <- wday(jan1)

  # 3. Calculate the date of the first Saturday of the epi_year
  # The first epi week ends on the first Saturday of January, provided it falls >= 4 days into the month
  # This logic is complex; a more reliable way is to find the date corresponding
  # to the *start* of the specific epi_week directly.

  # The standard approach is to work with the ISO week system (which is similar but starts on Monday)
  # and adjust for the Sunday start.
  
  # A robust approach using base R functionality with lubridate helpers:
  # The first day of the epi_week (Sunday) can be calculated relative to Jan 1.
  
  # Days to add to Jan 1 to get to the first Sunday of the year (this needs careful handling of year start logic)

  # A simplified robust calculation is to find the date of a known day of the week in that specific epi_week.
  # Let's find the Sunday (start date) of the target epi_week.

  # Use the `epiweek()` function's underlying logic reference:
  # The first epi week of the year ends on the first Saturday of January if it falls at least four days into the month.
  # Days from Jan 1 to first Saturday:
  days_to_first_sat <- (7 - jan1_wday + 7) %% 7
  first_sat <- jan1 + days(days_to_first_sat)
  
# 	if (length(mday(first_sat)) > 1) {
# 		print(paste0(year, ":", epi_week, ". first_sat is ",first_sat, sep=""))
# 	}
	
  # Check if first Sat is week 1 according to the rule (must be Jan 4th or later)
  if (as.numeric(mday(first_sat)[1]) < 4) {
    # If not week 1, the first epiweek starts the following Sunday
    start_date_week1 <- first_sat + days(1)
  } else {
    # Otherwise, week 1 starts the Sunday 6 days before this Saturday
    start_date_week1 <- first_sat - days(6)
  }

  # Calculate the start date of the desired epi_week
  # The start date of a specific epi_week is the start date of week 1 plus (week_num - 1) weeks
	target_date <- start_date_week1 + weeks(epi_week - 1)

  return(target_date)
}


format_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
		true = paste(day, month, years, sep = "\n"), 
		false = if_else(is.na(lag(month)) | lag(month) != month, true = paste(day, month, sep = "\n"), false = day)
	)
}


printy_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	paste(month, ". ", day, ", ", years, sep = "")
}


excel2df <- function(fname) { 

	# getting info about all excel sheets
	sheets <- readxl::excel_sheets(fname)
	tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
	data_frame <- lapply(tibble, as.data.frame)

	# assigning names to data frames
	names(data_frame) <- sheets

	# return the data frame
	return(data_frame)
} 

get_LTcoefs <- function(window_data) {
  # window_data is a data frame slice
  model <- lm(mean_abundance ~ epi_date, data = as.data.frame(window_data))
  return(coef(model))
}

