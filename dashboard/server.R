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
	getPlotData <- function(index, loc_id, date_win) {
		#print("getPlotData called!")
		#print(index)
		#print(loc_name)
		#print(date_win)
	
		#df_targ <- df_rs %>% filter(target == inputTarget & target_genetic_locus == inputLocus)
		if (loc_id == "WV") {
			df_this <- dflist_rs[[index]] %>% 
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE),
								 					date_primary := date_primary)
		} else {
			#print(loc_name)
			if (controlRV$activeGeoLevel == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_id == loc_id))$location_id
			}
	
			df_this <- dflist_rs[[index]] %>% 
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


	getSeqrDF <- function(inputTarget, loc_id, date_win) {
		#print("getSeqrDF called!")
	
		if (loc_id == "WV") {
# 			df_this <- df_seqr %>% 
# 								 group_by(date_to_plot, variant) %>% 
# 								 summarize(val := mean(percent, na.rm = TRUE)) %>% 
# 								 arrange(date_to_plot, variant)
			df_trans <- df_seqr
		} else {
	
			if (controlRV$activeGeoLevel[1] == "County") {
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

		return(df_return)
	}


	# Generate an abundance plot on the given data frame.
	#
	plotAbundance <- function(data_in, date_win) {
		#print("plotBasic called!")
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
	plotSeqr <- function(df_plot, date_win) {
		print("plotSeqr called!")
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

# Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in dplyr 1.1.0.
# Please use `reframe()` instead.
# When switching from `summarise()` to `reframe()`, remember that `reframe()` always returns an ungrouped data
#  frame and adjust accordingly.

    pdlist <- getPlotData(mapIndex, controlRV$mapClick[mapIndex], controlRV$viewMonths[mapIndex])    
		
		if (mapIndex == 1) {
			if (length(pdlist) > 0) {
				output$plot_covid <- renderPlotly({
					plotAbundance(pdlist, controlRV$viewMonths[mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
				})
			} else {
				output$plot_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
				print(paste0("Empty data frame for SARS-CoV-2!"))
			}
	
	
			output$plotsq_covid <- renderPlotly({
				df_plot <- getSeqrDF("SARS-CoV-2", controlRV$mapClick[mapIndex], controlRV$viewMonths[mapIndex])
				plotSeqr(df_plot, controlRV$viewMonths[mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
			})	
	
			# Update the plot titles
			#output$plot_title_covid <- renderText(paste0("SARS-CoV-2 Abundance for ", controlRV$mapClick[1], sep=""))
			output$plotsq_title_covid <- renderText(paste0("SARS-CoV-2 Variant Proportions", sep=""))
# 		for (i in 1:length(TARGETS_RS)) {
# 			eval(parse(text = paste0("output$plot", i, "_rs_title <- renderText(TARGETS_RS[", i, "])")))
# 		}
		}
		
  }

  
	#
	# Update the alert level data (reaction to map click or site selection).
	#
	updateAlertLevel <- function(geolevel, loc_id) {
	  #print("##### updateAlertLevelRS reached!")
	  
		if (missing(geolevel)) { geolevel = controlRV$activeGeoLevel }
		if (missing(loc_id)) { loc_id = controlRV$mapClick }
				
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
					
		# Print the title
		output$site_rs_title <- renderText(title_text)


		for (i in 1:length(TARGETS_RS)) {
			eval(parse(text = paste0("output$rs_CSL_", i, " <- renderText(DISEASE_RS[", i, "])")))
			
			#df_this <- df_this_site %>% filter(target == TARGETS_RS[i] & target_genetic_locus == GENLOCI_RS[i])
			#df_this_d <- df_regions %>% filter(target == TARGETS_RS[i] & target_genetic_locus == GENLOCI_RS[i])
			
			if (loc_id == "WV") {
			  delta <- mean((df_regions %>% filter(region_geolevel == "county"))[[DISEASE_RS[i]]], na.rm = TRUE)
				fresh <- mean((df_fresh %>% filter(region_geolevel == "county"))[[DISEASE_RS[i]]], na.rm = TRUE)
			} else {
			  delta <- (df_regions %>% filter(region_name == loc_id))[[DISEASE_RS[i]]]
				fresh <- (df_fresh %>% filter(region_name == loc_id))[[DISEASE_RS[i]]]
			}
#			delta <- calcDelta(df_this, 12)
#			fresh <- calcFresh(df_this)
					
			if (is.na(fresh)) {
			  eval(parse(text = paste0("output$rs_fresh_", i, " <- renderText('-')")))
			} else {
			  fresh <- formatC(as.numeric(fresh), format="d")
				eval(parse(text = paste0("output$rs_fresh_", i, " <- renderText('", fresh, "')")))
			}
			
			if (is.na(delta)) {
			  eval(parse(text = paste0("output$rs_delta_", i, " <- renderText('-')")))
			} else {
			  delta <- formatC(as.numeric(delta), format="d")
				eval(parse(text = paste0("output$rs_delta_", i, " <- renderText('", delta, " %')")))
			}

			color <- getFreshnessColor(fresh)
			frunner <- paste0("document.getElementById('rs_fresh_", i, "').style.color = '", color, "';")
			runjs(frunner)
		
			drunner <- paste0("document.getElementById('rs_delta_", i, "').style.color = '", color, "';")
			runjs(drunner)

			bgcolor <- getAlertColor(fresh, delta)
			arunner <- paste0("document.getElementById('rs_CSL_", i, "').style.backgroundColor = '", bgcolor, "';")
			runjs(arunner)
		}
	
	}
	
	
	#
	# Write the selection info block (reaction to map click or site selection).
	#
	updateSelectionInfo <- function(mapIndex) {
	  #print("##### updateSelectionInfo reached!")
	  
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
			
			title_text <- "State of West Virginia"
			selection_text <- paste0(
					"The State of West Virginia represents ", num_facilities, " active treatment facilities serving ", 
					prettyNum(total_popserved, big.mark=","), " residents across ", 
					num_counties, " counties.", sep="")
			plot_title <- "West Virginia Wastewater"
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
			
				title_text <- loc_name
				selection_text <- paste0("The ", loc_name, " facility serves ", 
					 prettyNum(total_popserved, big.mark=","), " residents in ", 
					 this_county, " county, WV.",
					 sep="")
				plot_title <- title_text

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

				title_text <- paste0(loc_id, " county, WV", sep="")
				selection_text <- paste0(loc_id, " county in WV supports ", 
					 num_facilities, " WaTCH ", facility_post, ", serving a total of ", 
					 prettyNum(total_popserved, big.mark=","), " residents (", 
					 prettyNum(pct_served, digits=1), "% of the county).",
					 sep="")
				plot_title <- title_text
			}
			
		}
		
		# Calculate the mean flow.
		#		
		end_date <- max(df_this$date_primary, na.rm=TRUE)
		dates <- c(end_date %m-% months(VIEW_RANGE_PRIMARY), end_date)
		mean_flow <- mean((df_this %>% filter(date_primary >= dates[1] & date_primary <= dates[2]))$sample_flow)
		
		flow_text <- paste0("Mean flow for this time period is ", prettyNum(mean_flow, digits = 2), 
				" MGD with a total reported capacity of ", prettyNum(total_cap, digits = 2),
				" MGD.", sep="")

		# Print the title, selection, and flow strings to the UI
		if (mapIndex == 1) {
			output$selection_title_covid <- renderText(title_text)
			output$selection_covid <- renderText(selection_text)
			output$selection_flow_covid <- renderText(flow_text)
			output$plot_title_covid <- renderText(paste0(TARGETS[mapIndex], " in ", plot_title))
			output$plotsq_title_covid <- renderText(paste0(TARGETS[mapIndex], " Variant Proportions"))
		} else {
			if (mapIndex == 2) {
				output$selection_title_flu <- renderText(title_text)
				output$selection_flu <- renderText(selection_text)
				output$selection_flow_flu <- renderText(flow_text)
				output$plot_title_flu <- renderText(title_text)
				output$plot_title_flu <- renderText(paste0(TARGETS[mapIndex], " in ", plot_title))
			} else {
				if (mapIndex == 3) {
					output$selection_title_rsv <- renderText(title_text)
					output$selection_rsv <- renderText(selection_text)
					output$selection_flow_rsv <- renderText(flow_text)
					output$plot_title_rsv <- renderText(title_text)
					output$plot_title_rsv <- renderText(paste0(TARGETS[mapIndex], " in ", plot_title))
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
	# Render the maps
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
		
		if (clicked == "county") {
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
			if (clicked == "facility") {
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
			updateAllPlots()
			updateSelectionInfo(mapIndex)
			updateAlertLevel(controlRV$activeGeoLevel[mapIndex], controlRV$mapClick[mapIndex])
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
			updateAlertLevel(controlRV$activeGeoLevel[mapIndex], controlRV$mapClick[mapIndex])

		} else {
			print(paste0("No data for ", clickedLocation))
		}
	}


	#
	# Click on a map shape
	#
	clickMapShape <- function(clicked, mapIndex) {
		if (length(input$map_covid_shape_click) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- input$map_covid_shape_click$id
			controlRV$mapClickLat[mapIndex] <- input$map_covid_shape_click$lat
			controlRV$mapClickLng[mapIndex] <- input$map_covid_shape_click$lng
    }

    
		if (clickedLocation %in% df_active_loc$location_id | clickedLocation %in% df_active_loc$location_counties_served | clickedLocation == "WV") {
			controlRV$mapClick[mapIndex] <- clickedLocation
		} else {
			clickedLocation <- "WV"
			controlRV$mapClick[mapIndex] <- clickedLocation
		}
		
    updateAllPlots(mapIndex)
		updateSelectionInfo(mapIndex)
		updateAlertLevel(controlRV$activeGeoLevel[mapIndex], controlRV$mapClick[mapIndex])

	}
	
	#
	# Off-target click on a map
	#
	clickMapOffMarker <- function(clicked, mapIndex) {
		# only respond if this click is in a new position on the map
		if (clicked$lat != controlRV$mapClickLat[mapIndex] | clicked$lng != controlRV$mapClickLng[mapIndex]) {
			#print("New position detected!")
			controlRV$mapClickLat[mapIndex] <- 0
			controlRV$mapClickLng[mapIndex] <- 0

			# Update the reactive element
			controlRV$mapClick[mapIndex] <- "WV"
		
			updateAllPlots(mapIndex)
			updateSelectionInfo(mapIndex)
			updateAlertLevel(controlRV$activeGeoLevel[mapIndex], controlRV$mapClick[mapIndex])
		}
	}
	

	###########################
	#
	# OBSERVER FUNCTIONS
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_covid_marker_click, { 
    #print("##### Map MARKER click top")
    clickMapMarker(input$map_covid_marker_click, 1)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_covid_shape_click, {
    #print("##### Map Shape Click!")
    clickMapShape(input$map_covid_shape_click, 1)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	#
  observeEvent(input$map_covid_click, { 
		#print("##### Map off-target click top")
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
		changeGeolevel(tolower(input$geo_level_covid), 1)
	}, ignoreInit = TRUE)

	#
	# Change the date range to view
	#
  observeEvent(input$view_range_covid, {

		# Update some reactive elements
		controlRV$viewMonths <- as.numeric(input$view_range_covid)
		
		# Update the plots
    updateAllPlots()

	}, ignoreInit = TRUE)


	#
	# Respond to 1 month trend line checkbox click.
	#
  observeEvent(input$trendline_1mo_rs, {

		# Update some reactive elements
		if (controlRV$trendLines[1] == FALSE) {
			controlRV$trendLines[1] <- TRUE
		} else {
			controlRV$trendLines[1] <- FALSE
		}
		
		# Update the plots
    updateAllPlots()

	}, ignoreInit = TRUE)

	#
	# Respond to 1 year trend line checkbox click.
	#
  observeEvent(input$trendline_1yr_rs, {

		# Update some reactive elements
		if (controlRV$trendLines[2] == FALSE) {
			controlRV$trendLines[2] <- TRUE
		} else {
			controlRV$trendLines[2] <- FALSE
		}
		
		# Update the plots
    updateAllPlots()

	}, ignoreInit = TRUE)


	# 
	# React to plot click
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_rs_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)


	onevent("click", "alert_scale", togglePanels(on=c("alert_scale_info")))

	
	#
	# React to click on the site focus panel close button
	#
	observeEvent(input$alert_scale_info_close,{
    togglePanels(off=c("alert_scale_info"))
	}, ignoreInit = TRUE)
	

})







#