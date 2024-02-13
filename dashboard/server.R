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
	# Assign Leaflet proxy
	#
	respLeafletProxy <- leafletProxy(mapId="map_resp", session)
	emergLeafletProxy <- leafletProxy(mapId="map_emerging", session)

	#
	# Initialize reactive values with some defaults
	#
	controlRV <- reactiveValues(
#								Dates = c(first_day, last_day),
								activeGeoLevel = c(GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT), 
								mapClick = c("WV", "WV"), 
								viewMonths = c(VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY), 
								clickLat=0, clickLng=0
	)
	
	
	# Generate a basic plot on the given data frame.
	#
	plotBasic <- function(df_plot, date_win) {
		#print("plotBasic called!")
    #View(df_plot)
    
		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(limits = c(today %m-% months(date_win), today)) + 
											#scale_color_manual(name = "Target", values = TARGET_COLORS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											#scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme() + 
											labs(x = NULL, y = NULL, color = NULL) + 
											geom_point(aes(x = date_to_plot, y = val), shape = 1, size = 2, alpha=0.9) + 
											geom_line(aes(x = date_to_plot, y = val), alpha=0.4, na.rm = TRUE)

		ggplotly(gplot) %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	# Generate a dataframe for a basic plot.
	#
	getBasicDF <- function(inputTarget, inputLocus, loc_name, date_win) {
	  #print("getBasicDF called!")
	  
		df_targ <- df_result %>% filter(target == inputTarget & target_genetic_locus == inputLocus)
		dates <- c(today %m-% months(date_win), today)
		
		if (loc_name == "WV") {
			df_plot <- df_targ %>% 
								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE))
		} else {
#			loc_id <- unique((resources$location %>% filter(location_common_name == loc_name))$location_id)
		
#			df_plot <- left_join(
#					df_targ %>% filter(location_id == loc_id & date_to_plot >= dates[1] & date_to_plot <= dates[2]), 
#					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
#					by="location_id")
			df_plot <- df_plot %>% 
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE))
		}
		#View(df_plot)
		return(df_plot)
	}
	

	# Generate a dataframe for a class plot
	getClassPlotData <- function(date_win) {
		
		targ <- controlRV$visibleTarget
		tclass <- controlRV$visibleClass
		
		df_related <- resources$target %>% filter(target_class == tclass)
		#df_targ <- df_result %>% filter(target %in% df_related$target_id & target != targ)
		df_targ <- df_result %>% filter(target %in% df_related$target_id)
		#print(df_targ)
		dates <- c(today %m-% months(date_win), today)
		
		loc_name <- controlRV$mapClick[1]
		if (loc_name == "WV") {
			df_plot <- df_targ %>% 
								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot, target, target_genetic_locus) %>% 
								 arrange(date_to_plot, target, target_genetic_locus) %>%
								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
		} else {
			loc_id <- unique((resources$location %>% filter(location_common_name == loc_name))$location_id)
		
			df_plot <- left_join(
					df_targ %>% filter(location_id == loc_id & date_to_plot >= dates[1] & date_to_plot <= dates[2]), 
					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
					by="location_id")
			df_plot <- df_plot %>% 
								 group_by(date_to_plot, target, target_genetic_locus) %>% 
								 arrange(date_to_plot, target, target_genetic_locus) %>%
								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
		}
		return(df_plot)
	}
	

	# Generate a dataframe for a plot of hospitalization data
	getHospPlotData <- function(date_win) {
		
		dates <- c(today %m-% months(date_win), today)
		#years <- c(lubridate::year(dates[1]), lubridate::year(dates[2]))
		#weeks <- c(lubridate::week(dates[1]), lubridate::week(dates[2]))
		
		df_plot <- df_hospital
		df_plot$date_to_plot <- as.Date(with(df_plot, paste(mmr_year, mmr_week, 1, sep="-")),"%Y-%U-%u")
		
		loc_name <- controlRV$mapClick[1]
		if (loc_name == "WV") {
			df_plot <- df_plot %>% 
								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 #filter(mmr_year >= years[1] & mmr_year <= years[2] & mmr_week >= weeks[1] & mmr_week <= weeks[2]) %>%
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(weekly_sum, na.rm = TRUE))
		} else {
			loc_id <- unique((df_location %>% filter(location_common_name == loc_name))$location_id)
		
			df_plot <- left_join(
					df_targ %>% filter(location_id == loc_id & date_to_plot >= dates[1] & date_to_plot <= dates[2]), 
					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
					by="location_id")
			df_plot <- df_plot %>% 
								 group_by(date_to_plot, target_genetic_locus) %>% 
								 arrange(date_to_plot, target_genetic_locus) %>%
								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
		}

		df_plot$target <- 'SARS-CoV-2'
		return(df_plot)
	}
	
  updatePlotsResp <- function() {
    df_plot1 <- getBasicDF(TARGETS_RESP[1], GENLOCI_RESP[1], controlRV$mapClick[1], controlRV$viewMonths[1])
    output$plot1_resp <- renderPlotly({plotBasic(df_plot1, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})
    
    df_plot2 <- getBasicDF(TARGETS_RESP[2], GENLOCI_RESP[2], controlRV$mapClick[1], controlRV$viewMonths[1])
    output$plot2_resp <- renderPlotly({plotBasic(df_plot2, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})
    
#    df_plot3 <- getBasicDF(TARGETS_RESP[3], GENLOCI_RESP[3], controlRV$mapClick[1], controlRV$viewMonths[1])
#    output$plot3_resp <- renderPlotly({plotBasic(df_plot3, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})
    
    df_plot4 <- getBasicDF(TARGETS_RESP[4], GENLOCI_RESP[4], controlRV$mapClick[1], controlRV$viewMonths[1])
    output$plot4_resp <- renderPlotly({plotBasic(df_plot4, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})

    # Update the plot titles
    output$plot1_resp_title = renderText(paste0(TARGETS_RESP[1], sep=""))
    output$plot2_resp_title = renderText(paste0(TARGETS_RESP[2], sep=""))
    output$plot3_resp_title = renderText(paste0(TARGETS_RESP[3], sep=""))
    output$plot4_resp_title = renderText(paste0(TARGETS_RESP[4], sep=""))
  }

  

	#
	# Render the WW table
	#
	output$tableWW <- renderDataTable(
		#df_result %>% filter(target == "SARS-CoV-2") %>% 
		#select(date_to_plot, location_id, target_genetic_locus, target_copies_flownorm_per_person) %>% 
		#rename(Date = date_to_plot, Location = location_id, Locus = target_genetic_locus, Copies = target_copies_flownorm_per_person) %>% 
	  resources$county %>% filter(county_lab_code != "") %>% select(county_name, county_population, county_lab_code) %>% rename(county = county_name, population = county_population, status = county_lab_code),
		options = list(paging = FALSE,   ## paginate the output
									 pageLength = 9,   ## number of rows to output for each page
									 scrollX = TRUE,   ## enable scrolling on X axis
									 scrollY = TRUE,   ## enable scrolling on Y axis
									 autoWidth = TRUE, ## use smart column width handling
									 server = FALSE,   ## use server- or client-side processing
									 dom = 't',
									 #buttons = c('csv', 'excel'),
									 columnDefs = list(list(targets = '_all', className = 'dt-center'))
		),
		#extensions = 'Buttons',
		selection = 'single',						## enable selection of a single row
		#filter = 'bottom',            	## include column filters at the bottom
		rownames = FALSE                ## don't show row numbers/names
	)	

	#
	# Render the plots
	#
	output$plot1_resp_title = renderText(paste0(TARGETS_RESP[1], sep=""))
	output$plot1_resp <- renderPlotly({
		df_plot <- getBasicDF(TARGETS_RESP[1], GENLOCI_RESP[1], controlRV$mapClick[1], controlRV$viewMonths[1])
		plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	

	output$plot2_resp_title = renderText(paste0(TARGETS_RESP[2], sep=""))
	output$plot2_resp <- renderPlotly({
	  df_plot <- getBasicDF(TARGETS_RESP[2], GENLOCI_RESP[2], controlRV$mapClick[1], controlRV$viewMonths[1])
	  plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	
	
	output$plot3_resp_title = renderText(paste0(TARGETS_RESP[3], sep=""))
#	output$plot3_resp <- renderPlotly({
#	  df_plot <- getBasicDF(TARGETS_RESP[3], GENLOCI_RESP[3], controlRV$mapClick[1], controlRV$viewMonths[1])
#	  plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
#	})	

	output$plot4_resp_title = renderText(paste0(TARGETS_RESP[4], sep=""))
	output$plot4_resp <- renderPlotly({
	  df_plot <- getBasicDF(TARGETS_RESP[4], GENLOCI_RESP[4], controlRV$mapClick[1], controlRV$viewMonths[1])
	  plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	

	#
	# Render the maps
	#
	output$map_resp <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircleMarkers(data = resources$location %>% filter(tolower(location_category) == "wwtp" & tolower(location_status) == "active"),
												 layerId = ~location_common_name, 
												 lat = ~location_lat, 
												 lng = ~location_lng, 
		#										 radius = 3500, 
												 radius = 5, 
												 stroke = FALSE,
												 weight = 4, 
												 opacity = 0.5,
		#										 color = ~alertPal(current_fold_change_smoothed), 
												 color = "black",
												 fill = TRUE,
		#										 fillColor = ~alertPal(current_fold_change_smoothed), 
												 fillColor = "gray",
												 group = "facility", 
												 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
												 fillOpacity = 0.6) %>%
			addPolygons( 
				data = county_spdf, 
				layerId = ~NAME, 
				fillColor = ~mypalette(POPCH_PCT), 
				stroke=TRUE, 
				fillOpacity = 0.7, 
				color="black", 
				weight=0.5, 
				group="county",
				label = ~NAME, 
				highlightOptions = highlightOptions(
					weight = 1,
					color = "#fff",
					#dashArray = "",
					fillOpacity = 0.9,
					bringToFront = TRUE)
			)		
	})


	output$map_emerg <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircleMarkers(data = resources$location %>% filter(tolower(location_category) == "wwtp" & tolower(location_status) == "active"),
												 layerId = ~location_common_name, 
												 lat = ~location_lat, 
												 lng = ~location_lng, 
		#										 radius = 3500, 
												 radius = 5, 
												 stroke = FALSE,
												 weight = 4, 
												 opacity = 0.5,
		#										 color = ~alertPal(current_fold_change_smoothed),
												 color = "black",  
												 fill = TRUE,
		#										 fillColor = ~alertPal(current_fold_change_smoothed), 
												 fillColor = "gray",
												 group = "facility", 
												 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
												 fillOpacity = 0.6) %>%
			addPolygons( 
				data = county_spdf, 
				layerId = ~NAME, 
				fillColor = ~mypalette(POPCH_PCT), 
				stroke=TRUE, 
				fillOpacity = 0.7, 
				color="black", 
				weight=0.5, 
				group="county",
				label = ~NAME, 
				highlightOptions = highlightOptions(
					weight = 1,
					color = "#fff",
					#dashArray = "",
					fillOpacity = 0.9,
					bringToFront = TRUE)
			)		
	})



	###########################
	#
	# OBSERVER FUNCTIONS, RESPIRATORY
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_resp_marker_click, { 
		#print("Map MARKER click top")
		
    clickedLocation <- input$map_resp_marker_click$id
		print(paste0("Facility:", clickedLocation))
		
		loc_id <- unique((resources$location %>% filter(location_common_name == clickedLocation))$location_id)

		if (loc_id %in% df_result$location_id) {
		
			# Update the reactive element
			controlRV$mapClick[1] <- clickedLocation
		
			# Update the map marker colors to indicate clicked marker
		  updatePlotsResp()
		  
		} else {
			print(paste0("No data for ", clickedLocation))
		}
				
  }, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_resp_shape_click, { 
    clickedLocation <- input$map_resp_shape_click$id
		print(paste0("County:", clickedLocation))
		
  }, ignoreInit = TRUE)

	#
	# Re-center the map
	#
	observeEvent(input$center_map_resp, {
		
		respLeafletProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
#		respLeafletProxy %>% flyTo(map_center$lng, map_center$lat, zoom = map_center$zoom, options = {animate = TRUE})
	}, ignoreInit = TRUE)

	#
	# Change the map active geolayer
	#
  observeEvent(input$geo_level_resp, {

		# Update some reactive elements
		controlRV$activeGeoLevel[1] <- input$geo_level_resp
		
		# Reconfigure the map markers and shapes
		if (tolower(input$geo_level_resp) == "county") {
			respLeafletProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircleMarkers(data = resources$location %>% filter(tolower(location_category) == "wwtp" & tolower(location_status) == "active"),
										 layerId = ~location_common_name, 
										 lat = ~location_lat, 
										 lng = ~location_lng, 
#										 radius = 3500, 
										 radius = 5, 
										 stroke = FALSE,
										 weight = 4, 
										 opacity = 0.5,
#										 color = ~alertPal(current_fold_change_smoothed), 
										 fill = TRUE,
#										 fillColor = ~alertPal(current_fold_change_smoothed), 
										 fillColor = "red",
										 group = "facility", 
										 #label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
										 fillOpacity = 0.6) %>%
					addPolygons( 
						data = county_spdf, 
						layerId = ~NAME, 
						fillColor = ~mypalette(POPCH_PCT), 
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
			if (tolower(input$geo_level_resp) == "facility") {
				respLeafletProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = county_spdf, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(POP2000), 
							stroke=TRUE, 
							fillOpacity = 0, 
							color="black", 
							weight=1.0, 
							group="county"
						) %>% 
						addCircleMarkers(data = resources$location %>% filter(tolower(location_category) == "wwtp"),
											 layerId = ~location_common_name, 
											 lat = ~location_lat, 
											 lng = ~location_lng, 
	#										 radius = 3500, 
											 radius = 5, 
											 stroke = FALSE,
											 weight = 4, 
											 opacity = 0.9,
	#										 color = ~alertPal(current_fold_change_smoothed), 
											 fill = TRUE,
	#										 fillColor = ~alertPal(current_fold_change_smoothed), 
											 fillColor = "red",
											 group = "facility", 
											 label = ~as.character(paste0(location_common_name, " (" , location_population_served, ")")), 
											 fillOpacity = 0.6)
			}
		}
		
		# Update the plots only if necessary
		if (controlRV$mapClick[1] != "WV") {
			controlRV$mapClick[1] <- "WV"

			# Update the plots to reflect the mapClick
			upodateAllPlots()
		}
		
	}, ignoreInit = TRUE)

	#
	# Change the date range to view
	#
  observeEvent(input$view_range_resp, {

		# Update some reactive elements
		controlRV$viewMonths[1] <- as.numeric(input$view_range_resp)
		
		# Update the plots
    updatePlotsResp()

	}, ignoreInit = TRUE)


	# 
	# React to plot click
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_resp_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)

})

