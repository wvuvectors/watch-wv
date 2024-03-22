library(readxl)

plot_theme <- function () { 
	theme(axis.text = element_text(size = 8),
				axis.title = element_text(size = 9, color="#333333"),
				axis.line.x = element_line(color="#bbbbbb", linewidth=1),
				axis.line.y = element_line(color="#bbbbbb", linewidth=1),
				axis.ticks.length.y = unit(-0.5, "cm"), 
				strip.text = element_text(size = 8),
				panel.grid.major = element_line(color="#eeeeee", linewidth=1), 
				panel.grid.minor.x = element_line(color="#eeeeee", linewidth=1),
				panel.background = element_rect(fill="transparent"), 
				panel.border = element_rect(fill=NA, color="#bbbbbb", linewidth=1), 
				legend.position = "none",
				legend.justification = c("left", "top"),
				#legend.direction = "horizontal",
				legend.box.just = "center",
				#legend.margin = margin(6, 6, 6, 6),
				legend.title = element_blank(),
				legend.background = element_rect(fill="transparent"), 
				legend.text = element_text(size = 8, color = "#333333"),
				plot.background = element_rect(fill="transparent"), 
				plot.title = element_text(size = 10, color="#045a8d", hjust=0.5)
)}

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

format_dates <- function(x) {
	month <- strftime(x, format = "%b")           		# Abbreviated name of the month.
	day <- strftime(x, format = "%d")           			# Abbreviated name of the day.
	years <- lubridate::year(x)                       # Year as a 4-digit number.
	if_else(is.na(lag(years)) | lag(years) != years,  # Conditions for pasting.
		true = paste(day, month, years, sep = "\n"), 
		false = if_else(is.na(lag(month)) | lag(month) != month, true = paste(day, month, sep = "\n"), false = day)
	)
}


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
  

calcTrend <- function(df_this, mo_base) {
	
	if (length(df_this$date_primary) == 0) {
		return(NA)
	}
	
	most_recent_date <- max(df_this$date_primary, na.rm = TRUE)

	vec_all <- (df_this %>% filter(date_primary > (most_recent_date %m-% months(mo_base))))$val

	if (length(vec_all) == 0) {
		trend <- NA
	} else {
		trend <- mean(vec_all, na.rm = TRUE)
		trend <- as.numeric(trend)
	}
	
	return(trend)
}


calcDelta <- function(df_this, mo_base) {
	
	if (length(df_this$date_primary) == 0) {
		return(NA)
	}
	most_recent_date <- max(df_this$date_primary, na.rm = TRUE)
	
	vec_now <- (df_this %>% filter(ymd(date_primary) == ymd(most_recent_date)))$target_copies_fn_per_cap
	vec_all <- (df_this %>% filter(date_primary > (most_recent_date %m-% months(mo_base))))$target_copies_fn_per_cap

	if (length(vec_now) == 0 | length(vec_all) == 0) {
		delta <- NA
	} else {
		d_now <- mean(vec_now, na.rm = TRUE)
		d_base <- mean(vec_all, na.rm = TRUE)
		delta <- 100 * (d_now - d_base)/d_base
		delta <- formatC(as.numeric(delta), format="d")
	}
	
	return(delta)
}


calcFresh <- function(df_this) {
	
	if (length(df_this$date_primary) == 0) {
		return(NA)
	}
	
	most_recent_date <- max(df_this$date_primary, na.rm = TRUE)
	
	freshness <- ymd(today)-ymd(most_recent_date)
	
	return(freshness)
}


#
# Generate an alert color string based on the target levels at the given location(s).
#
getAlertColor <- function(freshness, delta) {
	
# 		df_targ <- df_rs %>% filter(
# 			target == inputTarget & 
# 			target_genetic_locus == inputLocus & 
# 			location_id %in% locations)
	
	if (is.na(freshness) | as.numeric(freshness) >= STALE_DATA_THRESHOLDS[3]) {
		this_color <- ALERT_LEVEL_COLORS[5]
	} else {
		delta <- as.numeric(delta)
		this_color <- case_when(
			delta <= ALERT_LEVEL_THRESHOLDS[1] ~ ALERT_LEVEL_COLORS[1],
			delta > ALERT_LEVEL_THRESHOLDS[1] & delta <= ALERT_LEVEL_THRESHOLDS[2] ~ ALERT_LEVEL_COLORS[2],
			delta > ALERT_LEVEL_THRESHOLDS[2] & delta <= ALERT_LEVEL_THRESHOLDS[3] ~ ALERT_LEVEL_COLORS[3],
			delta >= ALERT_LEVEL_THRESHOLDS[3] ~ ALERT_LEVEL_COLORS[4]
		)
	}
	
	return(this_color)
}


#
# Generate a color string based on the data freshness at the given location(s).
#
getFreshnessColor <- function(freshness) {
	
	this_color <- case_when(
		is.na(freshness) ~ STALE_DATA_COLORS[4],
		freshness >= 0 & freshness < STALE_DATA_THRESHOLDS[1] ~ STALE_DATA_COLORS[1],
		freshness >= STALE_DATA_THRESHOLDS[1] & freshness < STALE_DATA_THRESHOLDS[2] ~ STALE_DATA_COLORS[2],
		freshness >= STALE_DATA_THRESHOLDS[2] & freshness < STALE_DATA_THRESHOLDS[3] ~ STALE_DATA_COLORS[3],
		freshness >= STALE_DATA_THRESHOLDS[3] ~ STALE_DATA_COLORS[4]
	)
	
	return(this_color)
}



