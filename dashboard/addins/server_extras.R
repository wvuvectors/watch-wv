	# Accepts a menu option.
	#
	# Change the way the active map is colored.
	#
	changeMapColor <- function(clicked) {
    #print("##### changeMapColor called!")

		mapProxy <- getMapProxy(controlRV$mapIndex)
		
		# Update some reactive elements
		controlRV$activeMapColor[controlRV$mapIndex] <- clicked
		
		colorByIndex <- getColorByIndex(clicked)
		
		mapProxy %>% clearMarkers() %>% clearShapes() %>% 
			addPolygons( 
				data = dflist_map_c[[colorByIndex]], 
				layerId = ~NAME,
				fill = TRUE,
				fillColor = ~abundance_color, 
				fillOpacity = 0.4, 
				stroke = TRUE, 
				color = ~trend_color, 
				weight = 2, 
				group="county",
				label = ~as.character(paste0(NAME, " county (", abundance_level, " & ", trend,")")), 
				highlightOptions = highlightOptions(
					weight = 1,
					color = "#00F900",
					#dashArray = "",
					fillOpacity = 0.5))
	}
	



	# Accepts a dataframe and an integer representing # of months to plot.
	#
	# Generate a ggplotly object of the variant data within the date window.
	#
	plotVariants <- function(df_variants, months_to_plot, target_index) {
		#print("##### plotVariants called!")

		dlab <- case_when(
			months_to_plot == 1 ~ DATE_LABELS[1],
			months_to_plot == 3 ~ DATE_LABELS[2],
			months_to_plot == 6 ~ DATE_LABELS[3],
			months_to_plot == 12 ~ DATE_LABELS[4],
			months_to_plot == 24 ~ DATE_LABELS[5]
		)

		dbrk <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS[1],
			months_to_plot == 3 ~ DATE_BREAKS[2],
			months_to_plot == 6 ~ DATE_BREAKS[3],
			months_to_plot == 12 ~ DATE_BREAKS[4],
			months_to_plot == 24 ~ DATE_BREAKS[5]
		)
				
		dbrk_minor <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS_MINOR[1],
			months_to_plot == 3 ~ DATE_BREAKS_MINOR[2],
			months_to_plot == 6 ~ DATE_BREAKS_MINOR[3],
			months_to_plot == 12 ~ DATE_BREAKS_MINOR[4],
			months_to_plot == 24 ~ DATE_BREAKS_MINOR[5]
		)

		#end_date <- max(df_variants$date_to_plot, na.rm=TRUE)
    end_date <- this_week
		date_limits <- c(end_date %m-% months(months_to_plot), end_date)
		df_plot <- df_variants %>% filter(date_to_plot >= date_limits[1] & date_to_plot <= date_limits[2])
    
    if (nrow(df_plot) > 0) {
			gplot <- ggplot(df_plot, aes(fill=color_group, y=total_pct, x=date_to_plot)) + labs(y = "", x = "") + 
	#							geom_bar(position="stack", stat="identity", aes(fill=factor(color_group), text=paste0("Week that starts ", printy_dates(date_to_plot), "\nVariant family: ", color_group, "\nProportion: ", prettyNum(total_pct, digits=2), "%", sep=""))) + 
								geom_bar(position="stack", stat="identity", aes(fill=factor(color_group), text=paste0("Week of ", printy_dates(date_to_plot-7), " - ", printy_dates(date_to_plot), "\nVariant family: ", color_group, "\nProportion: ", prettyNum(total_pct, digits=2), "%", sep=""))) + 
								scale_fill_brewer(type="div", palette = "RdYlBu", direction = -1, na.value = "#a8a8a8") + 
								labs(x="", y="", fill=NULL) + 
								scale_x_date(date_breaks = dbrk, date_minor_breaks = dbrk_minor, date_labels = dlab, limits = date_limits) + 
	#							scale_y_continuous(name = NULL, limits = c(0, 110), breaks = c(0, 25, 50, 75, 100)) +  DOESN'T WORK FOR SOME REASON!?
								plot_theme() + 
								theme(legend.position = "right", legend.title=element_blank())
		} else {
			gplot <- ggplot()
			fireDataWarnings(c(0))
		}

		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}

	# Accepts a location id.
	#
	# Return a dataframe of variant data for the current bio target at the given location.
	#
	getVariantData <- function(loc_id) {
	  #print("##### getVariantData called!")
	
		if (loc_id == "WV") {
			# get all active facilities
			loc_ids <- (df_active_loc %>% filter(location_category == "wwtp"))$location_id
		} else {
			if ((df_regions %>% filter(region_name == loc_id))$region_geolevel == "county") {
				# roll up county active facilities
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			} else {
				# just the single facility, but need it as a vector
				loc_ids <- (df_active_loc %>% filter(location_id == loc_id))$location_id
			}
		}
		
		df_trans <- df_seqr %>% filter(location_id %in% loc_ids)

		t1 <- df_trans %>% group_by(date_to_plot, color_group) %>% tally(variant_proportion)	# n = sum of variant prop across all samples in a date
		t2 <- df_trans %>% group_by(date_to_plot) %>% mutate(location_count = n_distinct(location_id, na.rm = TRUE)) %>% select(date_to_plot, location_count) %>% distinct()
		t3 <- df_trans %>% group_by(date_to_plot) %>% mutate(collection_count = n_distinct(sample_collection_end_datetime, na.rm = TRUE)) %>% select(date_to_plot, collection_count) %>% distinct()
		t23 <- merge(t2, t3, by = "date_to_plot")
		df_this <- merge(t1, t23, by = "date_to_plot")
		df_this$total_prop <- (df_this$n / df_this$location_count) / df_this$collection_count
		df_this$total_pct <- as.numeric(formatC(100*df_this$total_prop, format="f", digits=2))

		return(df_this)

	}


	updateVariantFreshness <- function() {
	  #print("##### updateDataFreshness called!")
	  
	  loc_id <- controlRV$mapClick[controlRV$mapIndex]	# This is a county name, or WV

		# Calculate the freshness and completeness of the variant data (COVID only).
		#
		df_fresh_sq <- df_seqr %>% filter(location_id %in% loc_ids)

		if (nrow(df_fresh_sq) > 0) {
			most_recent_sample_date_sq <- max(df_fresh_sq$date_to_plot, na.rm=TRUE)
			most_recent_sample_win_sq <- c(
				floor_date(most_recent_sample_date_sq, "week", week_start = 7), 
				floor_date(most_recent_sample_date_sq, "week", week_start = 7) + 6
			)
			epi_week_sq <- lubridate::epiweek(most_recent_sample_date_sq)

	
			sample_freshness_text_sq <- paste0(
				"week ", epi_week_sq, 
				" (", 
				printy_dates(most_recent_sample_win_sq[1]), " - ", printy_dates(most_recent_sample_win_sq[2]), 
				")",
				sep="")

			most_recent_contributors_sq <- df_fresh_sq %>% filter(date_to_plot == most_recent_sample_date_sq)
			most_recent_contributor_count_sq <- length(unique(most_recent_contributors_sq$location_id))
			if (most_recent_contributor_count_sq > 1) {
				site_suffix_sq <- "sites"
			} else {
				site_suffix_sq <- "site"
			}
			sample_completeness_text_sq <- paste0(
				most_recent_contributor_count_sq, " reporting ", site_suffix_sq, ".", 
				sep="")

		} else {
			sample_freshness_text_sq <- paste0("not reported for this region.", sep="")
			sample_completeness_text_sq <- ""
		}
		
				fr_text_sq <- paste0(
					"Latest variant data for this region is from ", 
					sample_freshness_text_sq, " and includes ", sample_completeness_text_sq, ".", 
					sep=""
				)
				cp_text_sq <- paste0("It includes samples from ", sample_completeness_text_sq, sep="")

	}
