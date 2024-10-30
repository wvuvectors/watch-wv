#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


shinyServer(function(input, output, session) {
	
	#
	# Assign Leaflet proxies
	#
	leafletProxyCovid <- leafletProxy(mapId="map_covid", session)
	leafletProxyFlu <- leafletProxy(mapId="map_flu", session)
	leafletProxyRsv <- leafletProxy(mapId="map_rsv", session)
	
	getMapProxy <- function(indx) {
		if (indx == 1) {
			return(leafletProxyCovid)
		}
		if (indx == 2) {
			return(leafletProxyFlu)
		}
		if (indx == 3) {
			return(leafletProxyRsv)
		}
		return(leafletProxyCovid)
	}
	
	
	#
	# Initialize reactive values with some defaults.
	#
	controlRV <- reactiveValues(
		activeGeoLevel = c(GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT), 
		mapClick = c("WV", "WV", "WV"), 
		trendLines = c(TRUE, TRUE),
		viewMonths = c(VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY), 
		mapClickLat = c(0, 0, 0),
		mapClickLng = c(0, 0, 0)
	)
	
	
	togglePanels <- function(on, off) {
		if (!missing(on)) {
			for (targ in on) {
				shinyjs::show(targ)
			}
		}
		if (!missing(off)) {
			for (targ in off) {
				shinyjs::hide(targ)
			}
		}
	}

	
	# Generate a dataframe for plotting target data.
	#
	getPlotData <- function(mapIndex, loc_id, date_win) {
#	  print("##### getPlotData called!")
	
		if (loc_id == "WV") {
			df_this <- dflist_rs[[mapIndex]] %>% 
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE),
								 					date_primary := date_primary)
		} else {
			#print(loc_name)
			if (controlRV$activeGeoLevel[mapIndex] == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_id == loc_id))$location_id
			}
	
			df_this <- dflist_rs[[mapIndex]] %>% 
								 filter(location_id %in% loc_ids) %>%
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE),
								 					date_primary := date_primary)
		}
	
		# Calculate the mean signal for each target over the last month and last year.
	#		if (controlRV$trendLines[1] == TRUE) {
			trendline_short <- calcTrend(df_this, 1)
	#		} else {
	#			trendline_short <- 0
	#		}
	#		if (controlRV$trendLines[2] == TRUE) {
			trendline_long <- calcTrend(df_this, VIEW_RANGE_PRIMARY)
	#		} else {
	#			trendline_long <- 0
	#		}

		end_date <- max(df_this$date_to_plot, na.rm=TRUE)
		dates <- c(end_date %m-% months(date_win), end_date)
	
		df_return <- df_this %>% filter(date_to_plot >= dates[1] & date_to_plot <= dates[2])
	
		if (length(df_return$date_to_plot) == 0) {
			return(list())
		} else {
			return(list(df_return, trendline_short, trendline_long))
		}
	}


	getSeqrData <- function(mapIndex, loc_id, date_win) {
#	  print("##### getSeqrData called!")
	
		if (loc_id == "WV") {
# 			df_this <- df_seqr %>% 
# 								 group_by(date_to_plot, variant) %>% 
# 								 summarize(val := mean(percent, na.rm = TRUE)) %>% 
# 								 arrange(date_to_plot, variant)
			df_trans <- df_seqr
		} else {
	
			if (controlRV$activeGeoLevel[mapIndex] == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_id == loc_id))$location_id
			}
			
			df_trans <- df_seqr %>% filter(location %in% loc_ids)
# 			df_this <- df_seqr %>% 
# 								 filter(location %in% loc_ids) %>%
# 								 group_by(date_to_plot, color_group) %>% 
# 								 summarize(val := sum(percent, na.rm = TRUE)) %>% 
# 								 arrange(date_to_plot, color_group)
		}

		t1 <- df_trans %>% group_by(date_to_plot, color_group) %>% tally(variant_proportion)	# n = sum of variant prop across all samples in a date
		t2 <- df_trans %>% group_by(date_to_plot) %>% mutate(location_count = n_distinct(location, na.rm = TRUE)) %>% select(date_to_plot, location_count) %>% distinct()
		t3 <- df_trans %>% group_by(date_to_plot) %>% mutate(collection_count = n_distinct(sample_collection_end_datetime, na.rm = TRUE)) %>% select(date_to_plot, collection_count) %>% distinct()
		t23 <- merge(t2, t3, by = "date_to_plot")
		df_this <- merge(t1, t23, by = "date_to_plot")
		df_this$total_prop <- (df_this$n / df_this$location_count) / df_this$collection_count
		df_this$total_pct <- as.numeric(formatC(100*df_this$total_prop, format="f", digits=2))

		end_date <- max(df_this$date_to_plot, na.rm=TRUE)
		dates <- c(end_date %m-% months(date_win), end_date)
	
		df_return <- df_this %>% filter(date_to_plot >= dates[1] & date_to_plot <= dates[2])

		#View(df_return)

		if (length(df_return$date_to_plot) == 0) {
			return(list())
		} else {
			return(list(df_return, -1, -1))
		}
	}


	# Generate an abundance plot on the given data frame.
	#
	plotAbundance <- function(data_in, date_win) {
#	  print("##### plotAbundance called!")

    df_plot <- as.data.frame(data_in[1])
    tl_mo <- as.numeric(data_in[2])
    tl_yr <- as.numeric(data_in[3])
    #View(df_plot)

		dlab <- case_when(
			date_win == 1 ~ DATE_LABELS[1],
			date_win == 3 ~ DATE_LABELS[2],
			date_win == 6 ~ DATE_LABELS[3],
			date_win == 12 ~ DATE_LABELS[4],
			date_win == 24 ~ DATE_LABELS[5],
		)

		dbrk <- case_when(
			date_win == 1 ~ DATE_BREAKS[1],
			date_win == 3 ~ DATE_BREAKS[2],
			date_win == 6 ~ DATE_BREAKS[3],
			date_win == 12 ~ DATE_BREAKS[4],
			date_win == 24 ~ DATE_BREAKS[5],
		)
		
		end_date <- max(df_plot$date_to_plot, na.rm=TRUE)
		
		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
											scale_x_date(date_breaks = dbrk, date_labels = dlab, limits = c(end_date %m-% months(date_win), end_date)) + 
											#scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme() + 
											labs(x = NULL, y = NULL, color = NULL) + 
											geom_hline(aes(yintercept=tl_mo), color=TRENDL_MO_COLOR, linetype="dotted", alpha=0.6, linewidth=0.5) + 
											geom_hline(aes(yintercept=tl_yr), color=TRENDL_YR_COLOR, linetype="solid", alpha=0.8, linewidth=1) + 
#											geom_point(aes(x = date_to_plot, y = val), shape = 1, size = 2, alpha=0.9) + 
#											geom_point(aes(x = date_to_plot, y = val, color = target_genetic_locus, text=paste0(date_to_plot, ": ", prettyNum(val, digits=1, big.mark=","), sep="")), shape = 1, size = 2, alpha=0.9) + 
											geom_point(aes(x = date_to_plot, y = val, text=date_to_plot), shape = 1, size = 2, alpha=0.9) + 
											geom_line(aes(x = date_to_plot, y = val), alpha=0.4, na.rm = TRUE)

		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	# Generate a seqr plot on the given data frame.
	#
	plotSeqr <- function(data_in, date_win) {
#		print("##### plotSeqr called!")

    df_plot <- as.data.frame(data_in[1])
    #tl_mo <- as.numeric(data_in[2])
    #tl_yr <- as.numeric(data_in[3])

		dlab <- case_when(
			date_win == 1 ~ DATE_LABELS[1],
			date_win == 3 ~ DATE_LABELS[2],
			date_win == 6 ~ DATE_LABELS[3],
			date_win == 12 ~ DATE_LABELS[4],
			date_win == 24 ~ DATE_LABELS[5],
		)

		dbrk <- case_when(
			date_win == 1 ~ DATE_BREAKS[1],
			date_win == 3 ~ DATE_BREAKS[2],
			date_win == 6 ~ DATE_BREAKS[3],
			date_win == 12 ~ DATE_BREAKS[4],
			date_win == 24 ~ DATE_BREAKS[5],
		)
    print(dbrk)
		
		#end_date <- max(df_plot$date_to_plot, na.rm=TRUE)
    #print(end_date)
    
		gplot <- ggplot(df_plot, aes(fill=color_group, y=total_pct, x=date_to_plot)) + labs(y = "", x = "") + 
							geom_bar(position="stack", stat="identity", aes(fill=factor(color_group))) + 
							scale_fill_brewer(type="div", palette = "RdYlBu", direction = -1, na.value = "#a8a8a8") + 
							labs(x="", y="") + 
		          scale_x_date(date_breaks = dbrk, date_labels = dlab) + 
		          #scale_x_date(date_breaks = dbrk, date_labels = dlab, limits = c(end_date %m-% months(date_win), end_date)) + 
#							scale_y_continuous(name = NULL, limits = c(0, 110), breaks = c(0, 25, 50, 75, 100)) +  DOESN'T WORK FOR SOME REASON!?
							plot_theme() +
							theme(legend.position = "right", legend.title=element_blank())

		ggplotly(gplot) %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}


  updateAllPlots <- function(mapIndex) {
#	  print("##### updateAllPlots called!")

# Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in dplyr 1.1.0.
# Please use `reframe()` instead.
# When switching from `summarise()` to `reframe()`, remember that `reframe()` always returns an ungrouped data
#  frame and adjust accordingly.

    pdlist <- getPlotData(mapIndex, controlRV$mapClick[mapIndex], controlRV$viewMonths[mapIndex])    
		
		loc_id <- controlRV$mapClick[mapIndex]
		
		if (loc_id == "WV") {
			plot_location <- "West Virginia"
		} else {
			if (controlRV$activeGeoLevel[mapIndex] == "County") {
				#this_county <- (df_active_loc %>% filter(location_id == loc_id))$location_counties_served
				this_state <- "WV"
				plot_location <- paste0(loc_id, " county, ", this_state, sep="")
			} else {
				if (controlRV$activeGeoLevel[mapIndex] == "Facility") {
					plot_location <- unique((df_active_loc %>% filter(location_id == loc_id))$location_common_name)
				}
			}
		}
		
		if (mapIndex == 1) {
			output$plot_title_covid <- renderText(paste0(TARGETS[mapIndex], " abundance in ", plot_location, " wastewater"))
			output$plotsq_title_covid <- renderText(paste0(TARGETS[mapIndex], " variants in ", plot_location, " wastewater"))

			if (length(pdlist) > 0) {
				output$plot_covid <- renderPlotly({
					plotAbundance(pdlist, controlRV$viewMonths[mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
				})
			} else {
				output$plot_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
				print(paste0("No data for ", TARGETS[mapIndex], "!"))
			}

	    sdlist <- getSeqrData(mapIndex, controlRV$mapClick[mapIndex], controlRV$viewMonths[mapIndex])    
			if (length(sdlist) > 0) {
				output$plotsq_covid <- renderPlotly({
					plotSeqr(sdlist, controlRV$viewMonths[mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
				})
			} else {
				output$plotsq_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
				print(paste0("No variant data for ", TARGETS[mapIndex], "!"))
			}
		} else {
			if (mapIndex == 2) {
				output$plot_title_flu <- renderText(paste0(TARGETS[mapIndex], " Abundance in ", plot_location))
				#output$plotsq_title_flu <- renderText(paste0(TARGETS[mapIndex], " Variant Proportions"))
			} else {
				if (mapIndex == 3) {
					output$plot_title_rsv <- renderText(paste0(TARGETS[mapIndex], " Abundance in ", plot_location))
					#output$plotsq_title_rsv <- renderText(paste0(TARGETS[mapIndex], " Variant Proportions"))
				}
			}
		}
		
  }

  
	#
	# Update the status row (on reaction to map click or site selection).
	#
	updateStatus <- function(mapIndex) {
#	  print("##### updateStatus called!")
	  
	  # Generate trend text
	  # Generate abundance level
	  # Generate data age
		geolevel <- controlRV$activeGeoLevel[mapIndex]
		loc_id <- controlRV$mapClick[mapIndex]
				
		if (loc_id == "WV") {
			
			#loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			#df_this_site <- df_rs
			title_text <- "State of West Virginia"

		} else {
		
			if (geolevel == "Facility") {
			
				loc_name <- unique((df_active_loc %>% filter(location_id == loc_id))$location_common_name)
				title_text <- paste0(loc_name, " Facility", sep="")

			} else {
				# county click!
				#loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
				#df_this_site <- df_rs %>% filter(location_id %in% loc_ids)
				title_text <- paste0(loc_id, " county, WV", sep="")
			}
		}
		
		if (loc_id == "WV") {
			delta <- mean((df_regions %>% filter(region_geolevel == "county"))[[DISEASES[mapIndex]]], na.rm = TRUE)
			fresh <- mean((df_fresh %>% filter(region_geolevel == "county"))[[DISEASES[mapIndex]]], na.rm = TRUE)
		} else {
			delta <- (df_regions %>% filter(region_name == loc_id))[[DISEASES[mapIndex]]]
			fresh <- (df_fresh %>% filter(region_name == loc_id))[[DISEASES[mapIndex]]]
		}
#			delta <- calcDelta(df_this, 12)
#			fresh <- calcFresh(df_this)
			
		if (is.na(fresh)) {
			output$fresh_covid <- renderText("-")
		} else {
			fresh <- formatC(as.numeric(fresh), format="d")
			output$fresh_covid <- renderText(fresh)
		}
		color <- getFreshnessColor(fresh)
#			frunner <- paste0("document.getElementById('rs_fresh_", i, "').style.color = '", color, "';")
		frunner <- paste0("document.getElementById('fresh_covid').style.color = '", color, "';")
		runjs(frunner)
	
		if (is.na(delta)) {
			output$level_covid <- renderText("-")
		} else {
			delta <- formatC(as.numeric(delta), format="d")
			output$level_covid <- renderText(paste0(delta, "%"))
		}

#			drunner <- paste0("document.getElementById('rs_delta_", i, "').style.color = '", color, "';")
#			drunner <- paste0("document.getElementById('delta_covid').style.color = '", color, "';")
#			runjs(drunner)

#			bgcolor <- getAlertColor(fresh, delta)
#			arunner <- paste0("document.getElementById('rs_CSL_", i, "').style.backgroundColor = '", bgcolor, "';")
#			arunner <- paste0("document.getElementById('rs_CSL_", i, "').style.backgroundColor = '", bgcolor, "';")
#			runjs(arunner)	
	}
	
	
	#
	# Write the selection info block (on reaction to map click or site selection).
	#
	updateSelectionInfo <- function(mapIndex) {
#	  print("##### updateSelectionInfo called!")
	  
		geolevel = controlRV$activeGeoLevel[mapIndex]
		loc_id = controlRV$mapClick[mapIndex]
		date_win = controlRV$viewMonths[mapIndex]
		
		if (loc_id == "WV") {
			
			loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			df_this <- df_rs
			
			total_popserved <- sum(distinct(df_active_loc %>% filter(location_category == "wwtp"), location_id, location_population_served)$location_population_served)
			if (total_popserved == -1) {
				total_popserved = "Unknown"
			}

			total_cap <- sum(distinct(df_active_wwtp, wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
			num_facilities <- n_distinct(df_this$location_id)
			num_counties <- n_distinct((df_active_loc %>% filter(location_category == "wwtp"))$location_counties_served)
			
			selection_text <- "West Virginia"
			risk_level_text <- "TBD"
			title_text <- "State of West Virginia"
			selection_details_text <- paste0(
					"The State of West Virginia represents ", num_facilities, " active treatment facilities serving ", 
					prettyNum(total_popserved, big.mark=","), " residents across ", 
					num_counties, " counties.", sep="")
		} else {
			
		  if (geolevel == "Facility") {
				
				loc_ids <- unique((df_active_loc %>% filter(location_id == loc_id))$location_id)
				this_wwtp_id <- unique((df_active_loc %>% filter(location_id == loc_id))$location_primary_wwtp_id)
				loc_name <- unique((df_active_loc %>% filter(location_id == loc_id))$location_common_name)
				
				df_this <- df_rs %>% filter(location_id == loc_ids)

				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id == loc_ids), location_id, location_population_served)$location_population_served)
				if (total_popserved == -1) {
					total_popserved = "Unknown"
				}

				total_cap <- sum(distinct(df_active_wwtp %>% filter(wwtp_id == this_wwtp_id), wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
				num_facilities <- 1
				this_county <- (df_active_loc %>% filter(location_id == loc_ids))$location_counties_served
			
				selection_text <- loc_name
				risk_level_text <- "TBD"
				title_text <- loc_name
				selection_details_text <- paste0("The ", loc_name, " facility serves ", 
					 prettyNum(total_popserved, big.mark=","), " residents in ", 
					 this_county, " county, WV.",
					 sep="")

			} else {
				# County click!
				
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
				num_facilities <- n_distinct(loc_ids)
				facility_post <- "facilities"
				if (num_facilities == 1) {
					facility_post = "facility"
				}
				
				df_this <- df_rs %>% filter(location_id %in% loc_ids)
				
				total_cap <- sum(distinct(df_active_wwtp %>% filter(wwtp_id %in% df_this), wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id %in% loc_ids), location_id, location_population_served)$location_population_served)
				pct_served <- 100 * total_popserved / (resources$county %>% filter(county_id == loc_id))$county_population

				selection_text <- paste0(loc_id, " county")
				risk_level_text <- "TBD"
				title_text <- paste0(loc_id, " county, WV", sep="")
				selection_details_text <- paste0(loc_id, " county in WV supports ", 
					 num_facilities, " WaTCH ", facility_post, ", serving a total of ", 
					 prettyNum(total_popserved, big.mark=","), " residents (", 
					 prettyNum(pct_served, digits=1), "% of the county).",
					 sep="")
			}
			
		}
		
		# Calculate the freshness of abundance data.
		#		
		most_recent_sample_date <- max(df_this$date_to_plot, na.rm=TRUE)
		most_recent_contributors <- df_this %>% filter(date_to_plot == most_recent_sample_date)
		most_recent_contributor_count <- length(unique(most_recent_contributors$location_id))
		if (most_recent_contributor_count > 1) {
			site_text <- "sites"
		} else {
			site_text <- "site"
		}
		sample_date_text <- paste0(
			"Most recent data is from the week of ", printy_dates(most_recent_sample_date), " and includes ", 
			most_recent_contributor_count, " reporting ", site_text, " (", 
			prettyNum(100*(most_recent_contributor_count/num_facilities), digits=1), "%).", 
			sep="")


		# Calculate the freshness of variant data (if COVID).
		#	
# 		most_recent_sample_date_sq <- max(df_this$date_to_plot, na.rm=TRUE)
# 		most_recent_contributors <- df_this %>% filter(date_to_plot == most_recent_sample_date)
# 		most_recent_contributor_count <- length(unique(most_recent_contributors$location_id))
# 		if (most_recent_contributor_count > 1) {
# 			site_text <- "sites"
# 		} else {
# 			site_text <- "site"
# 		}
# 		sample_date_text_sq <- paste0(
# 			"Most recent data is from the week of ", printy_dates(most_recent_sample_date), " and includes ", 
# 			most_recent_contributor_count, " reporting ", site_text, " (", 
# 			prettyNum(100*(most_recent_contributor_count/num_facilities), digits=1), "%).", 
# 			sep="")


		# Print the title, selection, and last_sample strings to the UI
		if (mapIndex == 1) {
			output$selection_covid <- renderText(selection_text)
			output$risk_level_covid <- renderText(risk_level_text)
			output$selection_details_covid <- renderText(selection_details_text)
			output$selection_title_covid <- renderText(title_text)
			output$selection_samples_covid <- renderText(sample_date_text)
#			output$selectionsq_samples_covid <- renderText(sample_date_text_sq)
		} else {
			if (mapIndex == 2) {
				output$selection_flu <- renderText(selection_text)
				output$risk_level_flu <- renderText(risk_level_text)
				output$selection_details_flu <- renderText(selection_details_text)
				output$selection_title_flu <- renderText(title_text)
				output$selection_samples_flu <- renderText(sample_date_text)
			} else {
				if (mapIndex == 3) {
					output$selection_rsv <- renderText(selection_text)
					output$risk_level_rsv <- renderText(risk_level_text)
					output$selection_title_rsv <- renderText(title_text)
					output$selection_details_rsv <- renderText(selection_details_text)
					output$selection_samples_rsv <- renderText(sample_date_text)
				}
			}
		}


	}

	
	getDownloadTableRS <- function(inputTarget, inputLocus, loc_name, date_win) {
		return(df_rs)
	}
		

	output$downloadData <- downloadHandler(
		filename = "test.csv",
		content = function(file) {
			write.csv(getDownloadTableRS(), file, row.names = FALSE)
		}
	)


	#
	#	Reset a map
	#
	resetMap <- function(mapIndex) {
		#print("resetMap called!")
		mapProxy <- getMapProxy(mapIndex)
		mapProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
	}
	
	#
	# Change the geolevel
	#
	changeGeolevel <- function(clicked, mapIndex) {

		mapProxy <- getMapProxy(mapIndex)

		# Update some reactive elements
		controlRV$activeGeoLevel[mapIndex] <- clicked
		
		if (clicked == "County") {
			mapProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
										 layerId = ~location_id, 
										 lat = ~location_lat, 
										 lng = ~location_lng, 
										 radius = ~dotsize, 
#										 radius = 5, 
										 stroke = FALSE,
										 weight = 4, 
										 opacity = 0.5,
#										 color = ~alertPal(current_fold_change_smoothed), 
										 fill = TRUE,
#										 fillColor = ~alertPal(current_fold_change_smoothed), 
										 fillColor = ~colorby,
										 group = "facility", 
										 label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
										 fillOpacity = 0.6) %>%
					addPolygons( 
						data = county_spdf, 
						layerId = ~NAME, 
						fillColor = ~colorby, 
						stroke=TRUE, 
						fillOpacity = 0.7, 
						color="black", 
						weight=0.5, 
						group="county",
						label = ~NAME, 
						highlightOptions = highlightOptions(
							weight = 1,
							color = "#fff",
							fillOpacity = 0.9,
							bringToFront = TRUE)
					)
		} else {
			if (clicked == "Facility") {
				mapProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = county_spdf, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(colorby), 
							stroke=TRUE,
							fillOpacity = 0, 
							color="#666666", 
							weight=1, 
							group="county"
						) %>% 
						addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
											layerId = ~location_id, 
											lat = ~location_lat, 
											lng = ~location_lng, 
											radius = ~dotsize, 
	#										radius = 5, 
											stroke = TRUE,
											weight = 2, 
											opacity = 0.9,
											color = "#000000", 
											fill = TRUE,
	#										fillColor = ~alertPal(current_fold_change_smoothed), 
											fillColor = ~colorby,
											group = "facility", 
											label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
											highlightOptions = highlightOptions(
												weight = 3,
												color = "#fff",
												fillOpacity = 0.9,
												bringToFront = TRUE),
											fillOpacity = 0.6)
			}
		}

		# Update the plots only if necessary
		if (controlRV$mapClick[mapIndex] != "WV") {
			controlRV$mapClick[mapIndex] <- "WV"
			updateAllPlots(mapIndex)
			updateSelectionInfo(mapIndex)
			updateStatus(mapIndex)
		}
		

	}

	#
	# Click on a map marker
	#
	clickMapMarker <- function(clicked, mapIndex) {
    if (length(clicked) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- clicked$id
    }
		
		loc_id <- unique((df_active_loc %>% filter(location_common_name == clickedLocation))$location_id)
    
		#print(loc_id)
		
		if (loc_id %in% df_rs$location_id | clickedLocation == "WV") {
		
			# Update the reactive element
			controlRV$mapClick[mapIndex] <- clickedLocation
		
		  updateAllPlots(mapIndex)
			updateSelectionInfo(mapIndex)
			updateStatus(mapIndex)

		} else {
			print(paste0("No data for ", clickedLocation))
		}
	}

	#
	# Click on a map shape
	#
	clickMapShape <- function(clicked, mapIndex) {
#    print("##### clickMapShape called!")
		if (length(clicked) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- clicked$id
			controlRV$mapClickLat[mapIndex] <- clicked$lat
			controlRV$mapClickLng[mapIndex] <- clicked$lng
    }

    
		if (clickedLocation %in% df_active_loc$location_id | clickedLocation %in% df_active_loc$location_counties_served | clickedLocation == "WV") {
			controlRV$mapClick[mapIndex] <- clickedLocation
		} else {
			clickedLocation <- "WV"
			controlRV$mapClick[mapIndex] <- clickedLocation
		}
		
    updateAllPlots(mapIndex)
		updateSelectionInfo(mapIndex)
		updateStatus(mapIndex)

	}
	
	#
	# Off-target click on a map
	#
	clickMapOffMarker <- function(clicked, mapIndex) {
#    print("##### clickMapOffMarker called!")
		# only respond if this click is in a new position on the map
		if (clicked$lat != controlRV$mapClickLat[mapIndex] | clicked$lng != controlRV$mapClickLng[mapIndex]) {
			#print("New position detected!")
			controlRV$mapClickLat[mapIndex] <- 0
			controlRV$mapClickLng[mapIndex] <- 0

			# Update the reactive element
			controlRV$mapClick[mapIndex] <- "WV"
		
			updateAllPlots(mapIndex)
			updateSelectionInfo(mapIndex)
			updateStatus(mapIndex)
		}
	}
	
	#
	# Render the map
	#
	output$map_covid <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
												 layerId = ~location_id, 
												 lat = ~location_lat, 
												 lng = ~location_lng, 
												 radius = ~dotsize, 
		#										 radius = 5, 
												 stroke = FALSE,
												 weight = 4, 
												 opacity = 0.5,
		#										 color = ~alertPal(current_fold_change_smoothed), 
												 color = "#000000",
												 fill = TRUE,
		#										 fillColor = ~alertPal(current_fold_change_smoothed), 
												 fillColor = ~colorby,
												 group = "facility", 
												 label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
												 fillOpacity = 0.6) %>%
			addPolygons( 
				data = county_spdf, 
				layerId = ~NAME, 
				fillColor = ~colorby, 
				stroke=TRUE, 
				fillOpacity = 0.7, 
				color="#000000", 
				weight=0.5, 
				group="county",
				label = ~NAME, 
				highlightOptions = highlightOptions(
					weight = 1,
					color = "#000000",
					#dashArray = "",
					fillOpacity = 0.9,
					bringToFront = TRUE)
			)		
	})


	
	
	###########################
	#
	# OBSERVER FUNCTIONS - COVID TAB
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_covid_marker_click, { 
#    print("##### Map MARKER click top")
    clickMapMarker(input$map_covid_marker_click, 1)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_covid_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_covid_shape_click, 1)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	#
  observeEvent(input$map_covid_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_covid_click, 1)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	#
	# Reset the map
	#
	observeEvent(input$map_reset_covid, {
		resetMap(1)
	}, ignoreInit = TRUE)

	#
	# Change the map active geolayer
	#
  observeEvent(input$geo_level_covid, {
		changeGeolevel(input$geo_level_covid, 1)
	}, ignoreInit = TRUE)

	#
	# Change the date range to view
	#
  observeEvent(input$view_range_covid, {

		# Update some reactive elements
		controlRV$viewMonths[1] <- as.numeric(input$view_range_covid)
		
		# Update the plots
    updateAllPlots(1)

	}, ignoreInit = TRUE)

	#
	# Open the color key popup
	#
	onevent("click", "risk_level_key_popup_covid", togglePanels(on=c("risk_level_key_popup_covid")))

	#
	# Close the color key popup
	#
	observeEvent(input$risk_level_key_popup_close_covid,{
    togglePanels(off=c("risk_level_key_popup_covid"))
	}, ignoreInit = TRUE)
	

	# 
	# React to plot click
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_rs_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)


})







#