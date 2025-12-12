
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
	data_frame
} 
  

calcTrend <- function(df_this, mo_base) {
	
	if (length(df_this$date_to_plot) == 0) {
		return(NA)
	}
	
	most_recent_date <- max(df_this$date_to_plot, na.rm = TRUE)

	vec_all <- (df_this %>% filter(date_to_plot > (most_recent_date %m-% months(mo_base))))$val

	if (length(vec_all) == 0) {
		trend <- NA
	} else {
		trend <- mean(vec_all, na.rm = TRUE)
		trend <- as.numeric(trend)
	}
	
	return(trend)
}



watchPal <- function(name) {
	if (tolower(name) == "abundance") {
		return(abundance_level_colors)
	} else if (tolower(name) == "trend") {
		return(trend_level_colors)
	} else if (tolower(name) == "lab") {
		return(lab_colors)
	} else {
		print(paste0("Problem with watchPal! name is ", name, sep=""))
		return(default_colors)
	}
}

# watch_colors <- list(
# 	Trend = trend_level_colors, 
# 	Abundance = abundance_level_colors
# )
# 
# 
# watch_palettes <- function(name, n, all_palettes = watch_colors, type = c("discrete", "continuous")) {
#   palette <- all_palettes[[name]]
#   if (missing(n)) {
#     n = length(palette)
#   }
#   type = match.arg(type)
#   out = switch(type,
#                continuous = grDevices::colorRampPalette(palette)(n),
#                discrete = palette[1:n]
#   )
#   structure(out, name = name, class = "palette")
# }
# 
# 
# scale_color_watch_d <- function(name) {
# 	ggplot2::scale_color_manual(values = watch_palettes(name, type = "discrete"))
# }
# scale_colour_watch_d = scale_color_watch_d
# 
# scale_fill_watch_d <- function(name) {
# 	ggplot2::scale_fill_manual(values = watch_palettes(name, type = "discrete"))
# }


