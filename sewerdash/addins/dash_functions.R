
plot_theme <- function () { 
	theme(axis.text = element_text(size = 8),
				axis.title = element_text(size = 9, color="#333333"),
				axis.line.x = element_line(color="#bbbbbb", linewidth=1),
				axis.line.y = element_line(color="#bbbbbb", linewidth=1),
				axis.ticks.length.y = unit(-0.5, "cm"), 
				strip.text = element_blank(),
				strip.background = element_rect(fill="#ffffff"),
#				strip.text = element_text(size = 8, color="#045a8d", hjust=0, vjust=0.5),
				panel.grid.major = element_line(color="#eeeeee", linewidth=1), 
#				panel.grid.minor.x = element_line(color="#eeeeee", linewidth=1),
				panel.grid.minor = element_line(color="#eeeeee", linewidth=0.7),
				panel.background = element_rect(fill="transparent"), 
				panel.border = element_rect(fill=NA, color="#bbbbbb", linewidth=1), 
				panel.spacing.y = unit(2, "lines"),
				legend.position = "none",
				legend.justification = c("left", "top"),
				#legend.direction = "horizontal",
				legend.box.just = "center",
				#legend.margin = margin(6, 6, 6, 6),
				legend.title = element_blank(),
				legend.background = element_rect(fill="transparent"), 
				legend.text = element_text(size = 8, color = "#333333"),
				plot.background = element_rect(fill="transparent"), 
				plot.title = element_text(size = 9, color="#045a8d", face="italic", hjust=0, vjust=0.5)
)}


# This is used in server.R, but NOT global.R so be careful.
isStale <- function(d) {

	date_diff <- as.numeric(difftime(this_week, d, units = "days"))
	#print(paste0(d, ": ", date_diff, sep=""))
	
	if (date_diff > STALE_THRESHOLD_DAYS) {
		return(TRUE)
	} else {
		return(FALSE)
	}
	
}


calcTrend <- function(df_this, mo_base) {
	
	if (length(df_this$primary_date) == 0) {
		return(NA)
	}
	
	most_recent_date <- max(df_this$primary_date, na.rm = TRUE)

	vec_all <- (df_this %>% filter(primary_date > (most_recent_date %m-% months(mo_base))))$mean_abundance

	if (length(vec_all) == 0) {
		trend <- NA
	} else {
		trend <- mean(vec_all, na.rm = TRUE)
		trend <- as.numeric(trend)
	}
	
	return(trend)
}


getAlertDetail <- function(disease, region, alevel, tlevel) {
	if (alevel$level == 1) {
		txt <- paste0(
		"Abundance data for ", disease, " in ", region, " is ", alevel$detail, " ", tlevel$detail, sep="")
	} else {
		txt <- paste0(
		"The level of ", disease, " in ", region, " wastewater is ", alevel$detail, " ", tlevel$detail, sep="")
	}
		
	return(txt)
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

