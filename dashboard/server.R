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
	rsLeafletProxy <- leafletProxy(mapId="map_rs", session)
	mgLeafletProxy <- leafletProxy(mapId="map_mg", session)
	stLeafletProxy <- leafletProxy(mapId="map_st", session)

	#
	# Initialize reactive values with some defaults.
	# Each element in a controlRV array corresponds to a UI tab.
	#
	controlRV <- reactiveValues(
								activeGeoLevel = c(GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT), 
								mapClick = c("WV", "WV", "WV"), 
								trendLines = c(TRUE, TRUE),
								viewMonths = c(VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY), 
								mapClickLat = c(0, 0, 0),
								mapClickLng = c(0, 0, 0),
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

	
	# Generate a dataframe for a basic RS plot.
	#
	getDataRS <- function(index, loc_name, date_win) {
		#print("getDataRS called!")
		#print(index)
		#print(loc_name)
		#print(date_win)
	
		#df_targ <- df_rs %>% filter(target == inputTarget & target_genetic_locus == inputLocus)
		if (loc_name == "WV") {
			df_this <- dflist_rs[[index]] %>% 
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE),
								 					date_primary := date_primary)
		} else {
			#print(loc_name)
			if (controlRV$activeGeoLevel[1] == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_common_name == loc_name))$location_id
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


	getSeqrDF <- function(inputTarget, loc_name, date_win) {
		#print("getBasicDF called!")
	
		if (loc_name == "WV") {
			df_this <- df_seqr %>% 
								 group_by(date_to_plot, color_group) %>% 
								 summarize(val := mean(percent, na.rm = TRUE)) %>% 
								 arrange(date_to_plot, color_group)
		} else {
	
			if (controlRV$activeGeoLevel[1] == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_common_name == loc_name))$location_id
			}

			df_this <- df_seqr %>% 
								 filter(location %in% loc_ids) %>%
								 group_by(date_to_plot, color_group) %>% 
								 summarize(val := sum(percent, na.rm = TRUE)) %>% 
								 arrange(date_to_plot, color_group)
		}
		#View(df_plot)

		end_date <- max(df_this$date_to_plot, na.rm=TRUE)
		dates <- c(end_date %m-% months(date_win), end_date)
	
		df_return <- df_this %>% filter(date_to_plot >= dates[1] & date_to_plot <= dates[2])

		return(df_return)
	}


	# Generate a basic plot on the given data frame.
	#
	plotBasic <- function(data_in, date_win) {
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
											geom_point(aes(x = date_to_plot, y = val), shape = 1, size = 2, alpha=0.9) + 
#											geom_point(aes(x = date_to_plot, y = val, color = target_genetic_locus, text=paste0(date_to_plot, ": ", prettyNum(val, digits=1, big.mark=","), sep="")), shape = 1, size = 2, alpha=0.9) + 
											geom_line(aes(x = date_to_plot, y = val), alpha=0.4, na.rm = TRUE)

		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	# Generate a seqr plot on the given data frame.
	#
	plotSeqr <- function(df_plot, date_win) {
		#print("plotSeqr called!")
    #View(df_plot)
    
		dlab <- case_when(
			controlRV$viewMonths[1] == 1 ~ DATE_LABELS[1],
			controlRV$viewMonths[1] == 3 ~ DATE_LABELS[2],
			controlRV$viewMonths[1] == 6 ~ DATE_LABELS[3],
			controlRV$viewMonths[1] == 12 ~ DATE_LABELS[4],
			controlRV$viewMonths[1] == 24 ~ DATE_LABELS[5],
		)

		dbrk <- case_when(
			controlRV$viewMonths[1] == 1 ~ DATE_BREAKS[1],
			controlRV$viewMonths[1] == 3 ~ DATE_BREAKS[2],
			controlRV$viewMonths[1] == 6 ~ DATE_BREAKS[3],
			controlRV$viewMonths[1] == 12 ~ DATE_BREAKS[4],
			controlRV$viewMonths[1] == 24 ~ DATE_BREAKS[5],
		)
		
		end_date <- max(df_plot$date_to_plot, na.rm=TRUE)

		gplot <- ggplot(df_plot, aes(fill=lineage_group, y=val, x=date_to_plot)) + labs(y = "", x = "") + 
							geom_bar(position="stack", stat="identity", aes(fill=factor(color_group))) + 
							scale_fill_brewer(type="qual", palette = "Dark2") + 
							labs(x="", y="") + 
							scale_x_date(date_breaks = dbrk, date_labels = dlab, limits = c(end_date %m-% months(date_win), end_date)) + 
							plot_theme() +
							theme(legend.position = "right", legend.title=element_blank())

		ggplotly(gplot) %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}


  updatePlotsRS <- function() {
    list1 <- getDataRS(1, controlRV$mapClick[1], controlRV$viewMonths[1])    
		if (length(list1) > 0) {
			output$plot1_rs <- renderPlotly({
				plotBasic(list1, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
			})
		} else {
			output$plot1_rs <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
			print(paste0("Empty data frame for ", TARGETS_RS[1], "!"))
		}

    list2 <- getDataRS(2, controlRV$mapClick[1], controlRV$viewMonths[1])
		if (length(list2) > 0) {
			output$plot2_rs <- renderPlotly({
				plotBasic(list2, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
			})
		} else {
			output$plot2_rs <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
			print(paste0("Empty data frame for ", TARGETS_RS[2], "!"))
		}
    
    list3 <- getDataRS(3, controlRV$mapClick[1], controlRV$viewMonths[1])
		if (length(list3) > 0) {
			output$plot3_rs <- renderPlotly({
				plotBasic(list3, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
			})
		} else {
			output$plot3_rs <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
			print(paste0("Empty data frame for ", TARGETS_RS[3], "!"))
		}
    
    list4 <- getDataRS(4, controlRV$mapClick[1], controlRV$viewMonths[1])
		if (length(list4) > 0) {
			output$plot4_rs <- renderPlotly({
				plotBasic(list4, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
			})
		} else {
			output$plot4_rs <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
			print(paste0("Empty data frame for ", TARGETS_RS[4], "!"))
		}

		output$plotsq_rs <- renderPlotly({
			df_plot <- getSeqrDF(TARGETS_RS[3], controlRV$mapClick[1], controlRV$viewMonths[1])
			plotSeqr(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		})	

    # Update the plot titles
		output$plotsq_rs_title <- renderText(paste0("SARS-CoV-2 Variant Proportions", sep=""))
		for (i in 1:length(TARGETS_RS)) {
			eval(parse(text = paste0("output$plot", i, "_rs_title <- renderText(TARGETS_RS[", i, "])")))
		}
  }

  
	#
	# Update the alert level data (reaction to map click or site selection).
	#
	updateAlertLevelRS <- function(geolevel, loc_name) {
		
		if (missing(geolevel)) { layer = controlRV$activeGeoLevel[1] }
		if (missing(loc_name)) { loc_name = controlRV$mapClick[1] }
		
		#print(loc_name)
		
		if (loc_name == "WV") {
			
			#loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			#df_this_site <- df_rs
			title_text <- "State of West Virginia"

		} else {
		
			if (geolevel == "Facility") {
			
				#loc_ids <- unique((df_active_loc %>% filter(location_common_name == loc_name))$location_id)
				#this_region_id <- unique((df_active_loc %>% filter(location_common_name == loc_name))$location_primary_wwtp_id)
			
				#df_this_site <- df_rs %>% filter(location_id == loc_ids)
				title_text <- loc_name

			} else {
				# county click!
				#loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
				#df_this_site <- df_rs %>% filter(location_id %in% loc_ids)
				title_text <- paste0(loc_name, " county, WV", sep="")
			}
		}
					
		# Print the title
		output$site_rs_title <- renderText(title_text)

		for (i in 1:length(TARGETS_RS)) {
			eval(parse(text = paste0("output$rs_CSL_", i, " <- renderText(DISEASE_RS[", i, "])")))
			
			#df_this <- df_this_site %>% filter(target == TARGETS_RS[i] & target_genetic_locus == GENLOCI_RS[i])
			#df_this_d <- df_regions %>% filter(target == TARGETS_RS[i] & target_genetic_locus == GENLOCI_RS[i])
			
			if (loc_name == "WV") {
				delta <- mean((df_regions %>% filter(region_geolevel == "county"))[[DISEASE_RS[i]]], na.rm = TRUE)
				fresh <- mean((df_fresh %>% filter(region_geolevel == "county"))[[DISEASE_RS[i]]], na.rm = TRUE)
			} else {
				delta <- (df_regions %>% filter(region_name == loc_name))[[DISEASE_RS[i]]]
				fresh <- (df_fresh %>% filter(region_name == loc_name))[[DISEASE_RS[i]]]
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
	# Write the site info block (reaction to map click or site selection).
	#
	updateSiteInfoRS <- function(geolevel, loc_name, date_win) {
		
		if (missing(geolevel)) { layer = controlRV$activeGeoLevel[1] }
		if (missing(loc_name)) { loc_name = controlRV$mapClick[1] }
		if (missing(date_win)) { date_win = controlRV$viewMonths[1] }
		
		#print(loc_name)

		if (loc_name == "WV") {
			
			loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			df_this <- df_rs
			
			total_popserved <- sum(distinct(df_active_loc %>% filter(location_category == "wwtp"), location_id, location_population_served)$location_population_served)
			if (total_popserved == -1) {
				total_popserved = "Unknown"
			}

			total_cap <- sum(distinct(df_active_wwtp, wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
			num_facilities <- n_distinct(df_this$location_id)
			num_counties <- n_distinct((df_active_loc %>% filter(location_category == "wwtp"))$location_counties_served)

			# Print the info
			#
			output$site_rs_info <- renderText(
				paste0("The State of West Virginia represents ", num_facilities, " active treatment facilities serving ", 
							 prettyNum(total_popserved, big.mark=","), " residents across ", 
							 num_counties, " counties.",
							 sep="")
			)

		} else {
			
			if (geolevel == "Facility") {
			
				loc_ids <- unique((df_active_loc %>% filter(location_common_name == loc_name))$location_id)
				this_wwtp_id <- unique((df_active_loc %>% filter(location_common_name == loc_name))$location_primary_wwtp_id)
			
				df_this <- df_rs %>% filter(location_id == loc_ids)

				title_text <- loc_name
				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id == loc_ids), location_id, location_population_served)$location_population_served)
				if (total_popserved == -1) {
					total_popserved = "Unknown"
				}

				total_cap <- sum(distinct(df_active_wwtp %>% filter(wwtp_id == this_wwtp_id), wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
				num_facilities <- 1
				this_county <- (df_active_loc %>% filter(location_id == loc_ids))$location_counties_served
			
				# Print the info
				#
				output$site_rs_info <- renderText(
					paste0("The ", loc_name, " facility serves ", 
								 prettyNum(total_popserved, big.mark=","), " residents in ", 
								 this_county, " county, WV.",
								 sep="")
				)
			} else {
				# county click!
				
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
				num_facilities <- n_distinct(loc_ids)
				facility_post <- "facilities"
				if (num_facilities == 1) {
					facility_post = "facility"
				}
				
				df_this <- df_rs %>% filter(location_id %in% loc_ids)
				
				title_text <- paste0(loc_name, " county, WV", sep="")

				total_cap <- sum(distinct(df_active_wwtp %>% filter(wwtp_id %in% df_this), wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id %in% loc_ids), location_id, location_population_served)$location_population_served)
				pct_served <- 100 * total_popserved / (resources$county %>% filter(county_id == loc_name))$county_population

				# Print the info
				#
				output$site_rs_info <- renderText(
					paste0(loc_name, " county in WV supports ", 
								 num_facilities, " WaTCH ", facility_post, ", serving a total of ", 
								 prettyNum(total_popserved, big.mark=","), " residents (", 
								 prettyNum(pct_served, digits=1), "% of the county).",
								 sep="")
				)
			}
			
		}
		
		# Calculate and print the site mean flow.
		#		
		end_date <- max(df_this$date_primary, na.rm=TRUE)
		dates <- c(end_date %m-% months(VIEW_RANGE_PRIMARY), end_date)
		
		mean_flow <- mean((df_this %>% filter(date_primary >= dates[1] & date_primary <= dates[2]))$sample_flow)
		output$site_rs_flow <- renderText(
			paste0("Mean daily flow is ", prettyNum(mean_flow, digits = 2), 
						 " MGD with a total reported capacity of ", prettyNum(total_cap, digits = 2),
						 " MGD.", sep="")
		)

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
	output$map_rs <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
												 layerId = ~location_common_name, 
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
												 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
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


	output$map_mg <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
												 layerId = ~location_common_name, 
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
												 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
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


	output$map_st <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
												 layerId = ~location_common_name, 
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
												 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
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
	# OBSERVER FUNCTIONS, ROUTINE SURVEILLANCE TAB
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_rs_marker_click, { 
		#print("Map MARKER click top")
		
 		if (length(input$map_rs_marker_click) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- input$map_rs_marker_click$id
    }
		
		loc_id <- unique((df_active_loc %>% filter(location_common_name == clickedLocation))$location_id)

		if (loc_id %in% df_rs$location_id | clickedLocation == "WV") {
		
			# Update the reactive element
			controlRV$mapClick[1] <- clickedLocation
		
			# Update the plots to include data for clicked marker
		  updatePlotsRS()
		  
		  # Update the site info panel
			updateSiteInfoRS(controlRV$activeGeoLevel[1], controlRV$mapClick[1], controlRV$viewMonths[1])
			updateAlertLevelRS(controlRV$activeGeoLevel[1], controlRV$mapClick[1])
		} else {
			#print(paste0("No data for ", clickedLocation))
		}
				
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_rs_shape_click, {
  	#print(input$map_rs_shape_click$id)
  	
		if (length(input$map_rs_shape_click) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- input$map_rs_shape_click$id
			controlRV$mapClickLat[1] <- input$map_rs_shape_click$lat
			controlRV$mapClickLng[1] <- input$map_rs_shape_click$lng
    }
		
		if (clickedLocation %in% df_active_loc$location_common_name | clickedLocation %in% df_active_loc$location_counties_served | clickedLocation == "WV") {
			controlRV$mapClick[1] <- clickedLocation
		} else {
			clickedLocation <- "WV"
			controlRV$mapClick[1] <- clickedLocation
		}
		
		updatePlotsRS()

		# Update the site info panel
		updateSiteInfoRS(controlRV$activeGeoLevel[1], clickedLocation, controlRV$viewMonths[1])
		updateAlertLevelRS(controlRV$activeGeoLevel[1], clickedLocation)

  }, ignoreNULL = FALSE, ignoreInit = FALSE)

	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	# We can tell 
	#
  observeEvent(input$map_rs_click, { 
		#print("Map click top")
		
		# only respond if this click is in a new position on the map
		if (input$map_rs_click$lat != controlRV$mapClickLat[1] | input$map_rs_click$lng != controlRV$mapClickLng[1]) {
			controlRV$mapClickLat[1] <- 0
			controlRV$mapClickLng[1] <- 0

			# Update the reactive element
			controlRV$mapClick[1] <- "WV"
		
			# Update the plots to show state-wide data
			updatePlotsRS()
			
			# Update the site info panel
			updateSiteInfoRS(controlRV$activeGeoLevel[1], "WV", controlRV$viewMonths[1])
			updateAlertLevelRS(controlRV$activeGeoLevel[1], "WV")
		}
				  				
  }, ignoreNULL = FALSE, ignoreInit = TRUE)



	#
	# Re-center the map
	#
	observeEvent(input$center_map_rs, {
		
		rsLeafletProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
#		rsLeafletProxy %>% flyTo(map_center$lng, map_center$lat, zoom = map_center$zoom, options = {animate = TRUE})
	}, ignoreInit = TRUE)

	#
	# Change the map active geolayer
	#
  observeEvent(input$geo_level_rs, {

		# Update some reactive elements
		controlRV$activeGeoLevel[1] <- input$geo_level_rs
		
		# Reconfigure the map markers and shapes
		if (tolower(input$geo_level_rs) == "county") {
			rsLeafletProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
										 layerId = ~location_common_name, 
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
										 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
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
			if (tolower(input$geo_level_rs) == "facility") {
				rsLeafletProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = county_spdf, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(colorby), 
							stroke=TRUE, 
							fillOpacity = 0, 
							color="black", 
							weight=1.0, 
							group="county"
						) %>% 
						addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
											 layerId = ~location_common_name, 
											 lat = ~location_lat, 
											 lng = ~location_lng, 
											 radius = ~dotsize, 
	#										 radius = 5, 
											 stroke = FALSE,
											 weight = 4, 
											 opacity = 0.9,
	#										 color = ~alertPal(current_fold_change_smoothed), 
											 fill = TRUE,
	#										 fillColor = ~alertPal(current_fold_change_smoothed), 
											 fillColor = ~colorby,
											 group = "facility", 
											 label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
											 fillOpacity = 0.6)
			}
		}
		
		# Update the plots only if necessary
		if (controlRV$mapClick[1] != "WV") {
			controlRV$mapClick[1] <- "WV"

			# Update the plots to reflect the mapClick
			updatePlotsRS()
		}
		
	}, ignoreInit = TRUE)

	#
	# Change the date range to view
	#
  observeEvent(input$view_range_rs, {

		# Update some reactive elements
		controlRV$viewMonths[1] <- as.numeric(input$view_range_rs)
		
		# Update the plots
    updatePlotsRS()

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
    updatePlotsRS()

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
    updatePlotsRS()

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







# WORKING!

# 	# Generate a dataframe for a plot of hospitalization data
# 	getHospPlotData <- function(date_win) {
# 		
# 		dates <- c(today %m-% months(date_win), today)
# 		#years <- c(lubridate::year(dates[1]), lubridate::year(dates[2]))
# 		#weeks <- c(lubridate::week(dates[1]), lubridate::week(dates[2]))
# 		
# 		df_plot <- df_hospital
# 		df_plot$date_primary <- as.Date(with(df_plot, paste(mmr_year, mmr_week, 1, sep="-")),"%Y-%U-%u")
# 		
# 		loc_name <- controlRV$mapClick[1]
# 		if (loc_name == "WV") {
# 			df_plot <- df_plot %>% 
# 								 filter(date_primary >= dates[1] & date_primary <= dates[2]) %>%
# 								 #filter(mmr_year >= years[1] & mmr_year <= years[2] & mmr_week >= weeks[1] & mmr_week <= weeks[2]) %>%
# 								 group_by(date_primary) %>% 
# 								 arrange(date_primary) %>%
# 								 summarize(val := mean(weekly_sum, na.rm = TRUE))
# 		} else {
# 			loc_id <- unique((df_location %>% filter(location_common_name == loc_name))$location_id)
# 		
# 			df_plot <- left_join(
# 					df_targ %>% filter(location_id == loc_id & date_primary >= dates[1] & date_primary <= dates[2]), 
# 					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
# 					by="location_id")
# 			df_plot <- df_plot %>% 
# 								 group_by(date_primary, target_genetic_locus) %>% 
# 								 arrange(date_primary, target_genetic_locus) %>%
# 								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
# 		}
# 
# 		df_plot$target <- 'SARS-CoV-2'
# 		return(df_plot)
# 	}
# 
#
# 	# Generate a dataframe for a class plot
# 	getClassPlotData <- function(date_win) {
# 		
# 		targ <- controlRV$visibleTarget
# 		tclass <- controlRV$visibleClass
# 		
# 		df_related <- resources$target %>% filter(target_class == tclass)
# 		#df_targ <- df_rs %>% filter(target %in% df_related$target_id & target != targ)
# 		df_targ <- df_rs %>% filter(target %in% df_related$target_id)
# 		#print(df_targ)
# 		dates <- c(today %m-% months(date_win), today)
# 		
# 		loc_name <- controlRV$mapClick[1]
# 		if (loc_name == "WV") {
# 			df_plot <- df_targ %>% 
# 								 filter(date_primary >= dates[1] & date_primary <= dates[2]) %>%
# 								 group_by(date_primary, target, target_genetic_locus) %>% 
# 								 arrange(date_primary, target, target_genetic_locus) %>%
# 								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
# 		} else {
# 			loc_id <- unique((resources$location %>% filter(location_common_name == loc_name))$location_id)
# 		
# 			df_plot <- left_join(
# 					df_targ %>% filter(location_id == loc_id & date_primary >= dates[1] & date_primary <= dates[2]), 
# 					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
# 					by="location_id")
# 			df_plot <- df_plot %>% 
# 								 group_by(date_primary, target, target_genetic_locus) %>% 
# 								 arrange(date_primary, target, target_genetic_locus) %>%
# 								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
# 		}
# 		return(df_plot)
# 	}
# 	
