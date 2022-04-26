#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


shinyServer(function(input, output, session) {
	
	###########################
	#
	#			GENERAL USE
	#
	###########################

	# Leaflet proxy
	watchLeafletProxy <- leafletProxy(mapId="watch_map", session)
	
	# Init reactive values with some defaults
	controlRV <- reactiveValues(
								Dates = c(min(df_watch$week_starting), max(df_watch$week_starting)),
								rollWin = SMOOTHER_DEFAULT,
								ci = "off",
								visibleTargets = TARGETS_DEFAULT,
								activeLayer = "WWTP"
	)
	
	wwtpRV <- reactiveValues(
								mapClick="State of West Virginia", 
								clickLat=0, clickLng=0
	)
	
	upstreamRV <- reactiveValues(
								mapClick="State of West Virginia", 
								clickLat=0, clickLng=0
	)

	# Render a base map
	generateMap <- function(data_in, center_lat, center_lng, zoom_level) {
		mymap = leaflet(data_in) %>% 
		addTiles() %>% 
		setView(lng = center_lng, lat = center_lat, zoom = zoom_level) %>% 
		addCircleMarkers(data = data_in %>% filter(group == "WWTP"),
										 layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 radius = 10, 
										 stroke = TRUE,
										 weight = 2, 
										 color = "black", 
										 fill = TRUE,
										 fillColor = "black",
										 group = "WWTP", 
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.6) %>%
		addCircleMarkers(data = data_in %>% filter(group == "Sewer Network"),
										 lat = ~latitude, 
										 lng = ~longitude, 
										 radius = 10, 
										 stroke = TRUE,
										 weight = 2, 
										 color = "blue", 
										 fill = TRUE,
										 fillColor = "blue",
										 group = "Sewer Network", 
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.6) %>% 
#		addProviderTiles(providers$Thunderforest.TransportDark)
#							Jawg.Streets
#							Esri.NatGeoWorldMap
#							Esri.WorldTopoMap
#		addPolylines(data=county_sf, fill=FALSE, weight=3, color="#999999", layerId="countiesLayer") %>%
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer") %>%
		addLayersControl(
			position = "bottomright",
			baseGroups = sort(unique(data_in$group), decreasing=TRUE),
									 options = layersControlOptions(collapsed = FALSE)
		) %>%
		hideGroup(unique(data_in$group)) %>% 
		showGroup("WWTP")
#		fitBounds(~-100,-60,~60,70) %>%
#		addLegend("bottomright", pal = cv_pal, values = ~cv_large_countries$deaths_per_million,
#			title = "<small>Deaths per million</small>")
		
		return(mymap)
	}

	controller <- function(targets, counter_targets) {
		if (!missing(targets)) {
			for (targ in targets) {
				shinyjs::show(targ)
			}
		}
		if (!missing(counter_targets)) {
			for (targ in counter_targets) {
				shinyjs::hide(targ)
			}
		}
	}
	



	###########################
	#
	#			WWTP PANEL
	#
	###########################
	
	#
	# Base map, starting state
	#
	output$watch_map <- renderLeaflet({
		generateMap(data_in = df_watch, center_lat = 38.938532, center_lng = -81.4222577, zoom_level = 7)
	})

	#
	# WWTP plot
	#
	plotLoad = function(layer, facility, dates, targets, rollWin, ci) {

		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = "State of West Virginia" }
		if (missing(dates)) { dates = controlRV$Dates }
		if (missing(targets)) { targets = controlRV$visibleTargets }
		if (missing(rollWin)) { rollWin = controlRV$rollWin }
		if (missing(ci)) { ci = controlRV$ci }
		
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
				
		df_plot <- df_watch %>%
							 filter(group == layer & week_starting >= fromDate & week_starting <= toDate) %>%
							 group_by(day)
		
		rollcols <- c()
		for (target in targets) {
			src_colname <- paste0(target, ".load", sep="")
			dest_colname <- paste0(src_colname, ".mean", sep="")

			roll_colname <- paste0(src_colname, ".roll", rollWin, sep="")
			rollcols <- append(rollcols, roll_colname)

			ci_colname <- paste0(roll_colname, ".ci90", sep="")

			df_loc <- df_watch %>% 
								filter(group == layer & week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>%
								summarize("{dest_colname}" := mean(.data[[src_colname]], na.rm = TRUE))

			zoo_loc <- zoo(df_loc[[dest_colname]], df_loc$day)
			zoo_mean <- rollmean(zoo_loc, rollWin, fill=NA, align="right")
			df_mean <- fortify(zoo_mean, melt=TRUE, names=c(Index="day", Value=roll_colname))
			df_mean <- select(df_mean, -c("Series"))

			df_plot <- left_join(df_plot, df_mean, by = c("day" = "day"), copy=TRUE)

			zoo_ci <- rollapply(zoo_loc, width=rollWin, fill=NA, align="right", FUN = ci90)
			df_ci <- fortify(zoo_ci, melt=TRUE, names=c(Index="day", Value=ci_colname))
			df_ci <- select(df_ci, -c("Series"))

			df_plot <- left_join(df_plot, df_ci, by = c("day" = "day"), copy=TRUE)
		}
		
		lims_x_date <- as.Date(strptime(c(fromDate, toDate), format = "%Y-%m-%d"))

#		target_pal <- colorFactor(
#			palette = c("blue", "green", "dark orange"),
#			domain = as.factor(result_cols),
#			na.color = "#aaaaaa",
#			alpha = TRUE
#		)
		
		gplot <- ggplot(df_plot) + labs(y = "Rolling Mean of RNA Mass Load", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = "2 weeks", labels = format_dates, limits = lims_x_date) + 
											my_theme() 
											#ggtitle("Weekly Mean COVID, All WV Treatment Facilities") + 

		for (rcol in rollcols) {
			gplot <- gplot + geom_point(aes(x = day, y = .data[[rcol]]), color = "blue", shape = 1, size = 1, alpha=0.5) + 
							 				 geom_line(aes(x = day, y = .data[[rcol]]), color = "blue")
			cicol <- paste0(rcol, ".ci90", sep="")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=.data[[rcol]], ymin=.data[[rcol]]-.data[[cicol]], ymax=.data[[rcol]]+.data[[cicol]]), alpha=0.2)
			}
		}
		ggplotly(gplot)
	}
	
	
	# Metadata block reactions to map click or site selection
	md_blockset <- function(layer, facility) {
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = "State of West Virginia" }
		
		if (facility == "State of West Virginia") {
			df_facility <- df_watch %>% filter(group == layer)
		} else {
			df_facility <- df_watch %>% filter(location_common_name == facility)
		}
		
		num_facilities = n_distinct(df_facility$location_common_name)
		total_cap = sum(distinct(df_facility, location_common_name, capacity_mgd)$capacity_mgd)+1
		total_popserved = sum(distinct(df_facility, location_common_name, population_served)$population_served)
		num_counties = n_distinct(df_facility$counties_served)
		total_county_pop = sum(distinct(df_facility, counties_served, county_population)$county_population)

		total_samples = n_distinct(df_facility$"Sample ID")
		date_first_sampled = min(df_facility$day)
		date_last_sampled = max(df_facility$day)
		
		mean_daily_flow = mean(df_facility$daily_flow)
		
		df_last28days <- df_facility %>% filter(ymd(day) >= ymd(today) - 28)
		mean_collfreq = (n_distinct(df_last28days$"Sample ID"))/4
		
		output$alert_level <- renderText(ALERT_TXT)
		output$site_signal <- renderText(TREND_TXT)

		output$scope <- renderText(facility)
		output$sample_count <- renderText(paste0(formatC(total_samples, big.mark=","), " samples since ", date_first_sampled))

		output$county_population <- renderText(formatC(total_county_pop, big.mark=","))
		output$population_served <- renderText(formatC(total_popserved, big.mark=","))
		output$population_served_pct <- renderText(paste0(formatC(100*total_popserved/total_county_pop, big.mark=","), "% of counties", sep=""))

		output$mean_flow <- renderText(paste0(formatC(mean_daily_flow, big.mark=","), " MGD", sep=""))
		output$collection_frequency <- renderText(paste0(mean_collfreq, " samples/week", sep=""))
		output$last_update <- renderText(paste0(ymd(today)-ymd(date_last_sampled)-1, " days ago", sep=""))

		if (facility == "State of West Virginia") {
			output$scope_count <- renderText(paste0(num_facilities, " facilities", sep=""))
			output$counties_served <- renderText(paste0(num_counties, " counties", sep=""))
		} else {
			output$scope_count <- renderText("1 facility")
			output$counties_served <- renderText(paste0(unique(df_facility$counties_served), " county", sep=""))
		}
	}
		
	
	#
	# Init the starting plot and metadata content
	#
	output$watch_plot <- renderPlotly({
		#print("Hi from output watch_plot!")

		md_blockset(
			layer = "WWTP", 
			facility = "State of West Virginia"
		)
		#print("md_blockset() just ran!")

		plotLoad() %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
	})	

	
	#
	# map interactivity
	#
	
	# Respond to layer change
	observe({
		selected_group <- req(input$watch_map_groups)
		controlRV.activeLayer <- selected_group
		if (selected_group == "WWTP") {
			output$data_format <- renderText("Data for WWTPs are calculated as mean RNA mass loads (copies of target adjusted for average daily flow).")
		} else {
			output$data_format <- renderText("Data for Sewer Network sites are calculated as mean target copies/L (daily flow is not available for these sites).")
		}
	})


	# Respond to off-marker click
  observeEvent(input$watch_map_click, { 
		#print("Map click event top")
		
		# only respond if this click is in a new position on the map
		if (input$watch_map_click$lat != wwtpRV$clickLat | input$watch_map_click$lng != wwtpRV$clickLng) {
			wwtpRV$clickLat <- 0
			wwtpRV$clickLng <- 0

			wwtpRV$mapClick <- "State of West Virginia"
			
			data_in <- df_watch
			# Update map marker colors
			watchLeafletProxy %>% 
				clearMarkers() %>% 
				addCircleMarkers(data = data_in %>% filter(group == "WWTP"),
												 layerId = ~location_common_name, 
												 lat = ~latitude, 
												 lng = ~longitude, 
												 radius = 10, 
												 stroke = TRUE,
												 weight = 2, 
												 color = "black", 
												 fill = TRUE,
												 fillColor = "black",
												 group = "WWTP", 
												 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
												 fillOpacity = 0.6) %>%
				addCircleMarkers(data = data_in %>% filter(group == "Sewer Network"),
												 lat = ~latitude, 
												 lng = ~longitude, 
												 radius = 10, 
												 stroke = TRUE,
												 weight = 2, 
												 color = "blue", 
												 fill = TRUE,
												 fillColor = "blue",
												 group = "Sewer Network", 
												 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
												 fillOpacity = 0.6)
		
			# Reset the plot
			output$watch_plot <- renderPlotly({
				plotLoad(facility = wwtpRV$mapClick) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
			})
		
			# Update reactive text elements
			md_blockset()
			
		}
	}, ignoreInit = TRUE)
	

	# Respond to click on WWTP map marker
  observeEvent(input$watch_map_marker_click, { 
    clickedLocation <- input$watch_map_marker_click$id
		wwtpRV$mapClick <- clickedLocation
		
		wwtpRV$clickLat <- input$watch_map_marker_click$lat
		wwtpRV$clickLng <- input$watch_map_marker_click$lng
		
		data_in <- df_watch

    # Update map marker colors
    watchLeafletProxy %>% 
    	clearMarkers() %>% 
			addCircleMarkers(data = data_in %>% filter(group == "WWTP"),
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 radius = 10, 
											 stroke = TRUE,
											 weight = 2, 
											 color=~ifelse(location_common_name==clickedLocation, yes = "#73FDFF", no = "black"), 
											 fill = TRUE,
											 fillColor = "black",
											 group = "WWTP", 
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.6) %>%
			addCircleMarkers(data = data_in %>% filter(group == "Sewer Network"),
											 lat = ~latitude, 
											 lng = ~longitude, 
											 radius = 10, 
											 stroke = TRUE,
											 weight = 2, 
											 color=~ifelse(location_common_name==clickedLocation, yes = "#73FDFF", no = "blue"), 
											 fill = TRUE,
											 fillColor = "blue",
											 group = "Sewer Network", 
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.6)

		# Update plots
		output$watch_plot <- renderPlotly({
				plotLoad(facility = wwtpRV$mapClick) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
		
		# Update reactive text elements
		md_blockset(facility = clickedLocation)
		
  }, ignoreInit = TRUE)







	###########################
	#
	# CONTROL PANEL EVENTS
	#
	###########################


	onevent("click", "alert_panel", controller(targets=c("alert_level_info")))
	onevent("click", "site_status_panel", controller(targets=c("site_status_info")))
	onevent("click", "scope_panel", controller(targets=c("scope_info")))
	onevent("click", "population_panel", controller(targets=c("population_info")))
	onevent("click", "networkpop_panel", controller(targets=c("networkpop_info")))
	onevent("click", "daily_flow_panel", controller(targets=c("daily_flow_info")))
	onevent("click", "collection_panel", controller(targets=c("collection_info")))
	onevent("click", "last_date_panel", controller(targets=c("last_date_info")))

	observeEvent(input$alert_level_info_close,{
    controller(counter_targets=c("alert_level_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$site_status_info_close,{
    controller(counter_targets=c("site_status_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$scope_info_close,{
    controller(counter_targets=c("scope_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$population_info_close,{
    controller(counter_targets=c("population_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$networkpop_info_close,{
    controller(counter_targets=c("networkpop_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$daily_flow_info_close,{
    controller(counter_targets=c("daily_flow_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$collection_info_close,{
    controller(counter_targets=c("collection_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$last_date_info_close,{
    controller(counter_targets=c("last_date_info"))
	}, ignoreInit = TRUE)
	

	#onevent("mouseenter", "site_status_panel", shinyjs::inlineCSS(list(.site_status_panel = "background-color: #FFFB00")))
	#onevent("mouseleave", "site_status_panel", shinyjs::inlineCSS(list(.site_status_panel = "background-color: #f0fff0")))

	#
	# Respond to change in targets to plot
	#
  observeEvent(input$targets_control, {
		controlRV$visibleTargets <- input$targets_control
		#updatePrettyCheckboxGroup(session = session, inputId = "targets_upstream", selected = controlRV$visibleTargets)
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = wwtpRV$mapClick) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

	#
	# Respond to change in dates to view
	#
	observeEvent(input$dates_control, {
		#print(paste0("Observed: ", input$dates_wwtp, sep = ""))
		#print(paste0("Min: ", min(df_watch$week_starting), ". Max: ", max(df_watch$week_starting), sep=""))
		
		controlRV$Dates <- input$dates_control
		#updateSliderTextInput(session = session, inputId = "dates_upstream", selected = controlRV$Dates)
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = wwtpRV$mapClick) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)


	#
	# Respond to change in rolling window size
	#
  observeEvent(input$roll_control, {
		controlRV$rollWin <- input$roll_control
		#updateSliderInput(session = session, inputId = "roll_upstream", value = controlRV$rollWin)
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = wwtpRV$mapClick) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)


})

