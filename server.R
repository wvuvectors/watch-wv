#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


shinyServer(function(input, output, session) {
	
	# Leaflet proxies
	myWWTPLeafletProxy <- leafletProxy(mapId="map_wwtp", session)
	myUpstreamLeafletProxy <- leafletProxy(mapId="map_upstream", session)
	
	# Init reactive values with soem defaults
	wwtpRV <- reactiveValues(mapClick="All Facilities", clickLat=0, clickLng=0, Dates=c(FIRST_DATE_WWTP, LAST_DATE_WWTP), rollWin=3, ci="off")
	upstreamRV <- reactiveValues(mapClick="All Locations", clickLat=0, clickLng=0, Dates=c(FIRST_DATE_UPSTREAM, LAST_DATE_UPSTREAM), rollWin=3, ci="off")



	###########################
	#
	#			PLOT RENDERING
	#
	###########################
	
	plotWWTP = function(dates, facility, rollwin, ci) {
		#print(ci)
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
		if (facility == "All Facilities") {
#			plot_df <- df_wwtp_s %>% filter(week_ending >= fromDate & week_ending <= toDate)
			loc_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>% 
								summarize(mean.rnamass = mean(rnamass, na.rm = TRUE))
			loc_zoo <- zoo(loc_df$mean.rnamass, loc_df$day)
		} else {
#			plot_df <- df_wwtp_wk %>% filter(week_ending >= fromDate & week_ending <= toDate & location_common_name == facility)
			loc_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate & location_common_name == facility)
			loc_zoo <- zoo(loc_df$rnamass, loc_df$day)
		}
		
		loc_mean_zoo <- rollmean(loc_zoo, rollwin, fill=NA, align="right")
		loc_mean_df <- fortify(loc_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.rnamass"))

#		loc_se_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = se)
#		loc_se_df <- fortify(loc_se_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.se"))
		
		loc_ci90_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc_ci90_df <- fortify(loc_ci90_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci90"))

		loc_ci95_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci95)
		loc_ci95_df <- fortify(loc_ci95_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci95"))

		loc_ci99_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci99)
		loc_ci99_df <- fortify(loc_ci99_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci99"))

		loc_join1_df <- left_join(loc_mean_df, loc_ci90_df, by = c("day" = "day"))
		loc_join2_df <- left_join(loc_join1_df, loc_ci95_df, by = c("day" = "day"))
		loc_join3_df <- left_join(loc_join2_df, loc_ci99_df, by = c("day" = "day"))
		loc_rolled_df <- left_join(loc_df, loc_join3_df, by = c("day" = "day"))
		loc_rolled_df <- select(loc_rolled_df, -c(Series.x, Series.y, Series.x.x, Series.y.y))
		
		lims_x_date <- as.Date(strptime(c(fromDate, toDate), format = "%Y-%m-%d"))

		gplot <- ggplot(loc_rolled_df, aes(x = day, y = rollmean.rnamass)) + 
											#labs(y = "Mean Weekly Mass Load", x = "") + 
											labs(y = "Rolling Mean of RNA Mass Load", x = "") + 
											geom_point(color="red", shape = 16, size = 1, alpha=0.5) + 
											geom_line(show.legend = FALSE, color = "#000000") + 
		#									geom_errorbar(aes(ymin=rollmean.rnamass-rollmean.se, ymax=rollmean.rnamass+rollmean.se), width=.2, position=position_dodge(0.05)) + 
		#									scale_y_continuous(limits=c(0,NA)) + 
											scale_y_continuous(labels = comma) + 
											#geom_ribbon(aes(ymin=rollmean.rnamass-rollmean.ci, ymax=rollmean.rnamass+rollmean.ci), alpha=0.2) + 
#											scale_x_date(breaks = "1 month", date_labels = '%d-%b-%Y', limits = lims_x_date) + 
											scale_x_date(breaks = "2 weeks", labels = format_dates, limits = lims_x_date) + 
											my_theme() 
											#scale_fill_manual(values = reds7) +
											#ggtitle("Weekly Mean COVID, All WV Treatment Facilities") + 
											#theme(axis.text.x = element_text(angle = 60, hjust=0.9))

		if (ci) {
			gplot <- gplot + 
							 #geom_ribbon(aes(ymin=rollmean.rnamass-rollmean.ci99, ymax=rollmean.rnamass+rollmean.ci99), fill="red", alpha=0.2) + 
							 #geom_ribbon(aes(ymin=rollmean.rnamass-rollmean.ci95, ymax=rollmean.rnamass+rollmean.ci95), fill="orange", alpha=0.2) + 
							 geom_ribbon(aes(ymin=rollmean.rnamass-rollmean.ci90, ymax=rollmean.rnamass+rollmean.ci90), alpha=0.2)
		}
		
		ggplotly(gplot)
	}

	plotUpstream = function(dates, facility, rollwin, ci) {
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
		if (facility == "All Locations") {
#			plot_df <- df_wwtp_s %>% filter(week_ending >= fromDate & week_ending <= toDate)
			loc_df <- df_upstream %>% filter(week_starting >= fromDate & week_starting <= toDate) %>% 
								group_by(day) %>% 
								summarize(mean.n1n2 = mean(n1n2, na.rm = TRUE))
			loc_zoo <- zoo(loc_df$mean.n1n2, loc_df$day)
		} else {
#			plot_df <- df_wwtp_wk %>% filter(week_ending >= fromDate & week_ending <= toDate & location_common_name == facility)
			loc_df <- df_upstream %>% filter(week_starting >= fromDate & week_starting <= toDate & location_common_name == facility)
			loc_zoo <- zoo(loc_df$n1n2, loc_df$day)
		}
		
		loc_mean_zoo <- rollmean(loc_zoo, rollwin, fill=NA, align="right")
		loc_mean_df <- fortify(loc_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.n1n2"))

#		loc_se_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = se)
#		loc_se_df <- fortify(loc_se_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.se"))

		loc_ci90_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc_ci90_df <- fortify(loc_ci90_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci90"))

		loc_ci95_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci95)
		loc_ci95_df <- fortify(loc_ci95_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci95"))

		loc_ci99_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci99)
		loc_ci99_df <- fortify(loc_ci99_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci99"))

		loc_join1_df <- left_join(loc_mean_df, loc_ci90_df, by = c("day" = "day"))
		loc_join2_df <- left_join(loc_join1_df, loc_ci95_df, by = c("day" = "day"))
		loc_join3_df <- left_join(loc_join2_df, loc_ci99_df, by = c("day" = "day"))
		loc_rolled_df <- left_join(loc_df, loc_join3_df, by = c("day" = "day"))
		loc_rolled_df <- select(loc_rolled_df, -c(Series.x, Series.y, Series.x.x, Series.y.y))

		lims_x_date <- as.Date(strptime(c(fromDate, toDate), format = "%Y-%m-%d"))

		gplot <- ggplot(loc_rolled_df, aes(x = day, y = rollmean.n1n2)) + 
											labs(y = "Rolling Mean in Copies/L", x = "") + 
											geom_point(color="red", shape = 16, size = 1, alpha=0.5) + 
											geom_line(show.legend = FALSE, color = "#000000") + 
		#									geom_errorbar(aes(ymin=rollmean.rnamass-rollmean.se, ymax=rollmean.rnamass+rollmean.se), width=.2, position=position_dodge(0.05)) + 
		#									scale_y_continuous(limits=c(0,NA)) + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = "2 weeks", labels = format_dates, limits = lims_x_date) + 
											#scale_x_date(breaks = "2 weeks", date_labels = '%d-%b-%Y', limits = lims_x_date) + 
											my_theme()
											#theme(axis.text.x = element_text(angle = 60, hjust=0.9))
		if (ci) {
			gplot <- gplot + 
							 geom_ribbon(aes(ymin=rollmean.n1n2-rollmean.ci90, ymax=rollmean.n1n2+rollmean.ci90), alpha=0.2)
		}
		
		ggplotly(gplot)

	}
	



	###########################
	#
	#			MAPS
	#
	###########################
	
	# Render a base map
	generateMap <- function(data_in, center_lat, center_lng, zoom_level) {
		mymap = leaflet(data_in) %>% 
		addTiles() %>% 
		setView(lng = center_lng, lat = center_lat, zoom = zoom_level) %>% 
		addCircleMarkers(layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 weight = 0.2, 
										 radius = 10, 
										 color = "#945200", 
										 label = ~as.character(paste0(location_common_name, " (" , type, ")")), 
										 fillOpacity = 0.6)
		
		return(mymap)
	}
	
	output$map_wwtp <- renderLeaflet({
		generateMap(data_in = df_wwtp, center_lat = 40.0632001, center_lng = -82.2772167, zoom_level = 7) %>% 
#		addProviderTiles(providers$Thunderforest.TransportDark)
#							Jawg.Streets
#							Esri.NatGeoWorldMap
#							Esri.WorldTopoMap
#		addPolylines(data=county_sf, fill=FALSE, weight=3, color="#999999", layerId="countiesLayer") %>%
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer")
#		addLayersControl(
#			position = "bottomright",
#			overlayGroups = c("COVID-19", "Influenza", "RSV"),
#											options = layersControlOptions(collapsed = FALSE)
#		) %>%
#		hideGroup(c("Influenza", "RSV"))
#		fitBounds(~-100,-60,~60,70) %>%
#		addLegend("bottomright", pal = cv_pal, values = ~cv_large_countries$deaths_per_million,
#			title = "<small>Deaths per million</small>")
	})


	# Respond to off-marker map click
  observeEvent(input$map_wwtp_click, { 
		#print("Map click event top")
		if (input$map_wwtp_click$lat != wwtpRV$clickLat | input$map_wwtp_click$lng != wwtpRV$clickLng) {
			wwtpRV$clickLat <- 0
			wwtpRV$clickLng <- 0

			wwtpRV$mapClick <- "All Facilities"
		
			# Update map marker colors
			myWWTPLeafletProxy %>% 
				clearMarkers() %>% 
				addCircleMarkers(data = df_wwtp,
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 weight = 0.2, 
											 radius = 10, 
											 color = "#945200", 
											 label = ~as.character(paste0(location_common_name, " (" , type, ")")), 
											 fillOpacity = 0.6)
		
			# Update control panel graph
			output$plot_wwtp <- renderPlotly({
				plotWWTP(c(FIRST_DATE_WWTP, LAST_DATE_WWTP), wwtpRV$mapClick, 3, FALSE) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
			})
		
			# Update reactive text elements
			output$facility_name <- renderText(paste0("All Facilities", sep = " "))
			output$facility_samples <- renderText("")
			output$facility_capacity <- renderText(paste0("Total statewide capacity involved: ", formatC(CAP_TOTAL_WWTP, big.mark=","), " million gallons per day (MGD).", sep = " "))
			output$facility_popserved <- renderText(paste0("Total estimated population included: ", formatC(POP_TOTAL_WWTP, big.mark=","), ".", sep = " "))
			output$facility_counties <- renderText(paste0("Serving ", CTY_TOTAL_WWTP, " total WV counties.", sep = " "))
		}
		#print("Map click event bottom")
		
	})
	

	# Respond to click on WWTP map marker
  observeEvent(input$map_wwtp_marker_click, { 
		#print("Marker click event top")
		#print(input$map_wwtp_marker_click$id)
    clickedLocation <- input$map_wwtp_marker_click$id
		wwtpRV$mapClick <- clickedLocation
		
		wwtpRV$clickLat <- input$map_wwtp_marker_click$lat
		wwtpRV$clickLng <- input$map_wwtp_marker_click$lng
		
    # Update map marker colors
    myWWTPLeafletProxy %>% 
    	clearMarkers() %>% 
			addCircleMarkers(data = df_wwtp,
										 layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 weight = 0.2, 
										 radius = 10, 
										 color=~ifelse(location_common_name==clickedLocation, yes = "red", no = "#945200"), 
										 label = ~as.character(paste0(location_common_name, " (" , type, ")")), 
										 fillOpacity = 0.6)
		
		# Update control panel graph
		output$plot_wwtp <- renderPlotly({
			plotWWTP(c(FIRST_DATE_WWTP, LAST_DATE_WWTP), wwtpRV$mapClick, 3, FALSE) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
		
		# Update reactive text elements
		watch_fac <- df_wwtp %>% filter(location_common_name == clickedLocation)
		SAMPLE_TOTAL_FAC = formatC(n_distinct(watch_fac$"Sample ID"), big.mark=",")
		FIRST_DATE_FAC = format(min(watch_fac$day), format = "%d %b %Y")
		LAST_DATE_FAC = format(max(watch_fac$day), format = "%d %b %Y")
		COLL_TYPE = unique(watch_fac$collection_scheme)

		output$facility_name <- renderText(paste0(clickedLocation, sep = " "))
		output$facility_samples <- renderText(paste0(SAMPLE_TOTAL_FAC, " samples collected as ", COLL_TYPE, "s from ", clickedLocation, " processed between ", FIRST_DATE_FAC," and ", LAST_DATE_FAC, ".", sep = " "))
		output$facility_capacity <- renderText(paste0("Max capacity of this facility: ", formatC(watch_fac$capacity_mgd[1], big.mark=","), " million gallons per day (MGD).", sep = " "))
		output$facility_popserved <- renderText(paste0("Estimated population included: ", formatC(watch_fac$population_served[1], big.mark=","), ".", sep = " "))
		output$facility_counties <- renderText(paste0("This facility is located in ", watch_fac$counties_served[1], " county WV.", sep = " "))
		#print("Marker click event bottom")
		
  })



	# Render the upstream map
	output$map_upstream <- renderLeaflet({
		generateMap(data_in = df_upstream, center_lat = 39.6352701, center_lng = -80.0125177, zoom_level = 13) %>% 
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer")

#		addLayersControl(
#			position = "bottomright",
#			overlayGroups = c("COVID-19", "Influenza", "RSV"),
#											options = layersControlOptions(collapsed = FALSE)
#		) %>%
#		hideGroup(c("Influenza", "RSV"))
		
	})

	# Respond to off-marker map click
  observeEvent(input$map_upstream_click, { 
		#print("Map click event top")
		if (input$map_upstream_click$lat != upstreamRV$clickLat | input$map_upstream_click$lng != upstreamRV$clickLng) {
			upstreamRV$clickLat <- 0
			upstreamRV$clickLng <- 0

			upstreamRV$mapClick <- "All Locations"
		
			# Update map marker colors
			myUpstreamLeafletProxy %>% 
				clearMarkers() %>% 
				addCircleMarkers(data = df_upstream,
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 weight = 0.2, 
											 radius = 10, 
											 color = "#945200", 
											 label = ~as.character(paste0(location_common_name, " (" , type, ")")), 
											 fillOpacity = 0.6)
		
			# Update control panel graph
			output$plot_upstream <- renderPlotly({
				plotUpstream(c(FIRST_DATE_UPSTREAM, LAST_DATE_UPSTREAM), upstreamRV$mapClick, 3, FALSE) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
			})
		
			# Update reactive text elements
			output$facility_name_upstream <- renderText(paste0("All Locations", sep = " "))
			output$facility_samples_upstream <- renderText("")
		}
		#print("Map click event bottom")
		
	})

	# Respond to click on upstream map marker
  observeEvent(input$map_upstream_marker_click, { 
    clickedLocation <- input$map_upstream_marker_click$id
		upstreamRV$mapClick <- clickedLocation

		upstreamRV$clickLat <- input$map_upstream_marker_click$lat
		upstreamRV$clickLng <- input$map_upstream_marker_click$lng
		
    # Update map marker colors
    myUpstreamLeafletProxy %>% 
    	clearMarkers() %>% 
			addCircleMarkers(data = df_upstream,
										 layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 weight = 0.2, 
										 radius = 10, 
										 color=~ifelse(location_common_name==clickedLocation, yes = "red", no = "#945200"), 
										 label = ~as.character(paste0(location_common_name, " (" , type, ")")), 
										 fillOpacity = 0.6)
		
		# Update control panel graph
		output$plot_upstream <- renderPlotly({
			plotUpstream(c(FIRST_DATE_UPSTREAM, LAST_DATE_UPSTREAM), upstreamRV$mapClick, 3, FALSE) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
		
		# Update reactive text elements
		watch_upstream <- df_upstream %>% filter(location_common_name == clickedLocation)
		SAMPLE_TOTAL_EDULOC = formatC(n_distinct(watch_upstream$"Sample ID"), big.mark=",")
		FIRST_DATE_EDULOC = format(min(watch_upstream$day), format = "%d %b %Y")
		LAST_DATE_EDULOC = format(max(watch_upstream$day), format = "%d %b %Y")
		COLL_TYPE_UP = unique(watch_upstream$collection_scheme)

#
		output$facility_name_upstream <- renderText(paste0(clickedLocation, sep = " "))
		output$facility_samples_upstream <- renderText(paste0(SAMPLE_TOTAL_EDULOC, " samples collected as ", COLL_TYPE_UP, "s from ", clickedLocation, " processed between ", FIRST_DATE_EDULOC," and ", LAST_DATE_EDULOC, ".", sep = " "))
#		output$facility_capacity <- renderText(paste0("Max capacity of this facility: ", formatC(watch_fac$capacity_mgd[1], big.mark=","), " million gallons per day (MGD).", sep = " "))
#		output$facility_popserved <- renderText(paste0("Estimated population included: ", formatC(watch_fac$population_served[1], big.mark=","), ".", sep = " "))
#		output$facility_counties <- renderText(paste0("This facility is located in ", watch_fac$counties_served[1], " county WV.", sep = " "))
		
  })





	###########################
	#
	#			WWTP CONTROL PANELS
	#
	###########################
	
	#
	# Render default control elements
	#
	output$plot_wwtp <- renderPlotly({
		plotWWTP(c(FIRST_DATE_WWTP, LAST_DATE_WWTP), wwtpRV$mapClick, 3, FALSE) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
	})

	output$last_update <- renderText(paste0("Wastewater Nowcast for West Virginia (", format(LAST_DATE_WWTP, format = "%d %b %Y"), ")", sep = " "))
	output$all_facilities <- renderText(paste0(FACILITY_TOTAL_WWTP, " active treatment facilities.", sep = " "))
	output$all_samples <- renderText(paste0(formatC(SAMPLE_TOTAL_WWTP, big.mark=","), " total samples processed since ", format(FIRST_DATE_WWTP, format = "%d %b %Y"), ".", sep = " "))
	output$facility_name <- renderText(paste0("All Facilities", sep = " "))

	output$facility_capacity <- renderText(paste0("Total statewide capacity involved: ", formatC(CAP_TOTAL_WWTP, big.mark=","), " million gallons per day (MGD).", sep = " "))
	output$facility_popserved <- renderText(paste0("Total estimated population included: ", formatC(POP_TOTAL_WWTP, big.mark=","), ".", sep = " "))
	output$facility_counties <- renderText(paste0("Serving ", CTY_TOTAL_WWTP, " total WV counties.", sep = " "))
	
	# Open the big plot
	observeEvent(input$embiggen_open_wwtp,{
#		print(paste0("Embiggen! ", input$embiggen_open_wwtp, sep=""))
    shinyjs::show(id = "conditionalPanelWWTP")

		output$plot_wwtp_big <- renderPlotly({
			plotWWTP(wwtpRV$Dates, wwtpRV$mapClick, wwtpRV$rollWin, wwtpRV$ci) %>% config(displayModeBar = FALSE) 
		})
	
		output$facility_name_big <- renderText(paste0(wwtpRV$mapClick, sep = " "))

	})
	

	#
	# Render large plot control elements
	#
	
	# Respond to change in plot dates on big plot
  observeEvent(input$plot_dates_wwtp, {
		#print(input$plot_dates)
		dates <- input$plot_dates_wwtp
		wwtpRV$Dates <- dates
		#print(dates)
		
		output$plot_wwtp_big <- renderPlotly({
			plotWWTP(wwtpRV$Dates, wwtpRV$mapClick, wwtpRV$rollWin, wwtpRV$ci) %>% config(displayModeBar = FALSE) 
		})
  })

	# Respond to change in rolling window size on big plot
  observeEvent(input$plot_roll_wwtp, {
		wwtpRV$rollWin <- input$plot_roll_wwtp
		
		output$plot_wwtp_big <- renderPlotly({
			plotWWTP(wwtpRV$Dates, wwtpRV$mapClick, wwtpRV$rollWin, wwtpRV$ci) %>% config(displayModeBar = FALSE) 
		})
  })


	# Respond to change in CI on big plot
  observeEvent(input$plot_ci_wwtp, {
	  #print(input$plot_ci_wwtp)
		wwtpRV$ci <- input$plot_ci_wwtp
		
		output$plot_wwtp_big <- renderPlotly({
			plotWWTP(wwtpRV$Dates, wwtpRV$mapClick, wwtpRV$rollWin, wwtpRV$ci) %>% config(displayModeBar = FALSE) 
		})
  })


	observeEvent(input$embiggen_close_wwtp,{
#		print(paste0("Embiggen! ", input$embiggen_open_wwtp, sep=""))
    wwtpRV$Dates <- c(FIRST_DATE_WWTP, LAST_DATE_WWTP)
    wwtpRV$rollWin <- 3
    #wwtpRV$ci <- FALSE
    shinyjs::hide(id = "conditionalPanelWWTP")
	})



	###########################
	#
	#			UPSTREAM CONTROL PANELS
	#
	###########################
	
	#
	# Render default control panel elements
	#
	output$plot_upstream <- renderPlotly({
		plotUpstream(c(FIRST_DATE_UPSTREAM, LAST_DATE_UPSTREAM), upstreamRV$mapClick, 3, FALSE)  %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
	})
	
	output$last_update_upstream <- renderText(paste0("Sewer Network Wastewater Nowcast (", format(LAST_DATE_UPSTREAM, format = "%d %b %Y"), ")", sep = " "))
	output$all_facilities_upstream <- renderText(paste0(FACILITY_TOTAL_UPSTREAM, " active sewer network sites.", sep = " "))
	output$all_samples_upstream <- renderText(paste0(formatC(SAMPLE_TOTAL_UPSTREAM, big.mark=","), " total samples processed since ", format(FIRST_DATE_UPSTREAM, format = "%d %b %Y"), ".", sep = " "))
	output$facility_name_upstream <- renderText(paste0("All Locations", sep = " "))

	# Open the big plot
	observeEvent(input$embiggen_open_upstream,{
#		print(paste0("Embiggen! ", input$embiggen_open_wwtp, sep=""))
    shinyjs::show(id = "conditionalPanelUpstream")

		output$plot_upstream_big <- renderPlotly({
			plotUpstream(upstreamRV$Dates, upstreamRV$mapClick, upstreamRV$rollWin, upstreamRV$ci) %>% config(displayModeBar = FALSE)
		})
	
		output$facility_name_big_upstream <- renderText(paste0(upstreamRV$mapClick, sep = " "))

	})
	

	#
	# Render large upstream plot control elements
	#
	
	# Respond to change in plot dates on big plot
  observeEvent(input$plot_dates_upstream, {
		#print(input$plot_dates)
		dates <- input$plot_dates_upstream
		upstreamRV$Dates <- dates
		#print(dates)
		
		output$plot_upstream_big <- renderPlotly({
			plotUpstream(upstreamRV$Dates, upstreamRV$mapClick, upstreamRV$rollWin, upstreamRV$ci)  %>% config(displayModeBar = FALSE)
		})
  })

	# Respond to change in rolling window size on big plot
  observeEvent(input$plot_roll_upstream, {
		upstreamRV$rollWin <- input$plot_roll_upstream
		
		output$plot_upstream_big <- renderPlotly({
			plotUpstream(upstreamRV$Dates, upstreamRV$mapClick, upstreamRV$rollWin, upstreamRV$ci) %>% config(displayModeBar = FALSE) 
		})
  })


	# Respond to change in CI on big plot
  observeEvent(input$plot_ci_upstream, {
	  #print(input$plot_ci_wwtp)
		upstreamRV$ci <- input$plot_ci_upstream
		
		output$plot_upstream_big <- renderPlotly({
			plotUpstream(upstreamRV$Dates, upstreamRV$mapClick, upstreamRV$rollWin, upstreamRV$ci) %>% config(displayModeBar = FALSE) 
		})
  })


	observeEvent(input$embiggen_close_upstream,{
#		print(paste0("Embiggen! ", input$embiggen_open_wwtp, sep=""))
    upstreamRV$Dates <- c(FIRST_DATE_UPSTREAM, LAST_DATE_UPSTREAM)
    upstreamRV$rollWin <- 3
    #upstreamRV$ci <- FALSE
    shinyjs::hide(id = "conditionalPanelUpstream")
	})



})

