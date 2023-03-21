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
	watchLeafletProxy <- leafletProxy(mapId="watch_map", session)

	#
	# Initialize reactive values with some defaults
	#
	controlRV <- reactiveValues(
#								Dates = c(first_day, last_day),
								activeGeoLevel = GEOLEVELS_DEFAULT, 
								visibleTarget = TARGETS_DEFAULT, 
								mapClick = "WV", 
								widePlotMonths = 12, 
								narrowPlotMonths = 2, 
								clickLat=0, clickLng=0
	)

#	getTrendData <- function(target) {
#		
#	}

	
	#
	# Generate a basic plot on the given data frame
	#
	plotOne <- function(df_plot, date_win) {
		#print("plotOne called!")

		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(limits = c(today %m-% months(date_win), today)) + 
											#scale_color_manual(name = "Target", values = TARGET_COLORS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											#scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme() + 
											labs(x = NULL, y = NULL, color = NULL) + 
											geom_point(aes(x = date_to_plot, y = val, color = target_genetic_locus), shape = 1, size = 2, alpha=0.9) + 
											geom_line(aes(x = date_to_plot, y = val, color = target_genetic_locus), alpha=0.4, na.rm = TRUE)
		ggplotly(gplot)# %>% layout(legend = list(orientation = "h"))
	}
	

	# Generate a dataframe for a basic plot
	getPlotData <- function(date_win) {
		
		targ <- controlRV$visibleTarget
		df_targ <- df_result %>% filter(target == targ)
		dates <- c(today %m-% months(date_win), today)
		
		loc_name <- controlRV$mapClick
		if (loc_name == "WV") {
			df_plot <- df_targ %>% 
								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot, target_genetic_locus) %>% 
								 arrange(date_to_plot, target_genetic_locus) %>%
								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
		} else {
			loc_id <- unique((df_location %>% filter(location_common_name == loc_name))$location_id)
		
			df_plot <- left_join(
					df_targ %>% filter(location_id == loc_id & date_to_plot >= dates[1] & date_to_plot <= dates[2]), 
					df_location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
					by="location_id")
			df_plot <- df_plot %>% 
								 group_by(date_to_plot, target_genetic_locus) %>% 
								 arrange(date_to_plot, target_genetic_locus) %>%
								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
		}
		return(df_plot)
	}
	
	#
	# Render the WW table
	#
	output$tableWW <- renderDataTable(
		#df_result %>% filter(target == "SARS-CoV-2") %>% 
		#select(date_to_plot, location_id, target_genetic_locus, target_copies_flownorm_per_person) %>% 
		#rename(Date = date_to_plot, Location = location_id, Locus = target_genetic_locus, Copies = target_copies_flownorm_per_person) %>% 
		df_county %>% filter(county_lab_code != "") %>% select(county_name, county_population, county_lab_code) %>% rename(county = county_name, population = county_population, status = county_lab_code),
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
	# Render some text titles and legends
	#
	output$plotww_title = renderText(paste0(controlRV$visibleTarget, " in ", controlRV$mapClick, " wastewater", sep=""))

	output$plotww_locus_legend = renderUI(
					HTML(
					as.character(
						span(
							style="height: 20px; width: 20px; margin: 5px; color: #c5000b; background-color: #c5000b; border: 1 px solid black; border-radius: 3px;", 
							paste0("XX"))),
					as.character(span(style="font-size: 11px;", paste0("N1"))),
					as.character(
						span(
							style="font-size: 11px; height: 20px; width: 30px; color: #ffffff; background-color: #ffffff;", 
							paste0("XXX"))),
					as.character(
						span(
							style="height: 20px; width: 20px; margin: 5px; color: #c5000b; background-color: #c5000b; border: 1 px solid black; border-radius: 3px;", 
							paste0("XX"))),
					as.character(span(style="font-size: 11px;", paste0("N2")))
					)
	)

	output$plotww_wide_title = renderText(paste0("Past ", controlRV$widePlotMonths, " months", sep=""))

	output$plotww_narrow_title = renderText(paste0("Past ", controlRV$narrowPlotMonths, " months", sep=""))


	#
	# Render the pathogen plots
	#
	output$plotww_wide <- renderPlotly({
		df_plot <- getPlotData(controlRV$widePlotMonths)
		plotOne(df_plot, controlRV$widePlotMonths) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
	})	

	output$plotww_narrow <- renderPlotly({
		df_plot <- getPlotData(controlRV$narrowPlotMonths)
		plotOne(df_plot, controlRV$narrowPlotMonths) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
	})	


	#
	# Render the map
	#
	output$watch_map <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircleMarkers(data = df_location %>% filter(tolower(location_category) == "wwtp"),
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
					#dashArray = "",
					fillOpacity = 0.9,
					bringToFront = TRUE)
			)		
	})



	###########################
	#
	# OBSERVER FUNCTIONS
	#
	###########################



	#
	# Re-center the map (layer dependent)
	#
	observeEvent(input$center_map, {
		
		watchLeafletProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
#		watchLeafletProxy %>% flyTo(map_center$lng, map_center$lat, zoom = map_center$zoom, options = {animate = TRUE})
	}, ignoreInit = TRUE)

	
	#
	# Change the active map layer
	#
  observeEvent(input$geo_level, {

		# Update some reactive elements
		controlRV$activeGeoLevel <- input$geo_level
		
		# Reconfigure the map markers and shapes
		if (tolower(input$geo_level) == "county") {
			watchLeafletProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircleMarkers(data = df_location %>% filter(tolower(location_category) == "wwtp"),
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
			if (tolower(input$geo_level) == "facility") {
				watchLeafletProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = county_spdf, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(POPCH_PCT), 
							stroke=TRUE, 
							fillOpacity = 0, 
							color="black", 
							weight=1.0, 
							group="county"
						) %>% 
						addCircleMarkers(data = df_location %>% filter(tolower(location_category) == "wwtp"),
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
		if (controlRV$mapClick != "WV") {
			controlRV$mapClick <- "WV"

			# Update the plots to reflect the mapClick
			df_plotw <- getPlotData(controlRV$widePlotMonths)
			output$plotww_wide <- renderPlotly({plotOne(df_plotw, controlRV$widePlotMonths) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})

			df_plotn <- getPlotData(controlRV$narrowPlotMonths)
			output$plotww_narrow <- renderPlotly({plotOne(df_plotn, controlRV$narrowPlotMonths) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})

			# Update the plot title
			output$plotww_title = renderText(paste0(controlRV$visibleTarget, " in ", controlRV$mapClick, " wastewater", sep=""))
		}
		
	}, ignoreInit = TRUE)


	# 
	# React to map marker click
	#
  observeEvent(input$watch_map_marker_click, { 
		#print("Map MARKER click top")
		
    clickedLocation <- input$watch_map_marker_click$id
		print(paste0("Facility:", clickedLocation))
		
		loc_id <- unique((df_location %>% filter(location_common_name == clickedLocation))$location_id)

		if (loc_id %in% df_result$location_id) {
		
			# Update the reactive element
			controlRV$mapClick <- clickedLocation
		
			# Update the map marker colors to indicate clicked marker
		
			# Update the base plots
			df_plotw <- getPlotData(controlRV$widePlotMonths)
			output$plotww_wide <- renderPlotly({plotOne(df_plotw, controlRV$widePlotMonths) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})

			df_plotn <- getPlotData(controlRV$narrowPlotMonths)
			output$plotww_narrow <- renderPlotly({plotOne(df_plotn, controlRV$narrowPlotMonths) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})

			# Update the plot title
			output$plotww_title = renderText(paste0(controlRV$visibleTarget, " in ", controlRV$mapClick, " wastewater", sep=""))
		} else {
			print(paste0("No results for ", clickedLocation))
		}
				
  }, ignoreInit = TRUE)


	# 
	# React to map shape click
	#
  observeEvent(input$watch_map_shape_click, { 
    clickedLocation <- input$watch_map_shape_click$id
		print(paste0("County:", clickedLocation))
		
  }, ignoreInit = TRUE)


})

