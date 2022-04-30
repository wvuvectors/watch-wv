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
#								Dates = c(min(df_watch$week_starting), max(df_watch$week_starting)),
								Dates = c(first_day, last_day),
								rollWin = SMOOTHER_DEFAULT,
								ci = FALSE,
								visibleTargets = TARGETS_DEFAULT,
								activeLayer = "WWTP",
								mapClick = "State of West Virginia",
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
										 layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 radius = 10, 
										 stroke = TRUE,
										 weight = 2, 
										 color = "dark orange", 
										 fill = TRUE,
										 fillColor = "dark orange",
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
	


	#
	# Base map, starting state
	#
	output$watch_map <- renderLeaflet({
		center <- MAP_CENTERS %>% filter(layer == "WWTP")
		generateMap(data_in = df_watch, center_lat = center$lat, center_lng = center$lng, zoom_level = center$zoom)
	})

	#
	# Main plot
	#
	plotLoad = function(layer, facility, dates, targets, rollWin, ci) {
		#print("plotLoad called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		if (missing(dates)) { dates = controlRV$Dates }
		if (missing(targets)) { targets = controlRV$visibleTargets }
		if (missing(rollWin)) { rollWin = controlRV$rollWin }
		if (missing(ci)) { ci = controlRV$ci }
		
#		print(layer)
#		print(facility)
#		print(dates)
#		print(targets)
#		print(rollWin)
#		print(ci)
		
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
		
		if (facility == "State of West Virginia") {
			df_plot <- df_watch %>%
								 filter(group == layer & day >= fromDate & day <= toDate) %>%
								 group_by(day)
		} else {
		
			df_plot <- df_watch %>%
								 filter(location_common_name == facility & day >= fromDate & day <= toDate) %>%
								 group_by(day)
		}
		
		df_facility <- df_plot
		anchor_date <- ymd(max(df_plot$day) - 28)
		
		for (target in targets) {
			if (layer == "Sewer Network") {
				src_colname <- target
			} else {
				src_colname <- paste0(target, ".load", sep="")
			}
			dest_colname <- paste0(src_colname, ".mean", sep="")

			roll_colname <- paste0(src_colname, ".roll", rollWin, sep="")
			ci_colname <- paste0(roll_colname, ".ci", sep="")
						
			if (facility == "State of West Virginia") {
				df_loc <- df_watch %>% 
									filter(group == layer & day >= fromDate & day <= toDate) %>%
									group_by(day) %>%
									summarize("{dest_colname}" := mean(.data[[src_colname]], na.rm = TRUE))
			} else {
				df_loc <- df_watch %>% 
									filter(location_common_name == facility & day >= fromDate & day <= toDate) %>%
									group_by(day) %>%
									summarize("{dest_colname}" := mean(.data[[src_colname]], na.rm = TRUE))
			}
			
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
		
		date_step <- "2 weeks"
		if (toDate - fromDate < 60) {
			date_step <- "2 days"
		}
		if (toDate - fromDate < 15) {
			date_step <- "1 day"
		}

		target_pal <- colorFactor(
			palette = TARGETS_DF$target_color,
			domain = TARGETS_DF$target_value,
			ordered = TRUE,
			na.color = "#aaaaaa",
			alpha = TRUE
		)
		
		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = date_step, labels = format_dates, limits = lims_x_date) + 
											my_theme()

		for (target in targets) {
			if (layer == "Sewer Network") {
				rcol <- paste0(target, ".roll", rollWin, sep="")
				src_col <- target
			} else {
				rcol <- paste0(target, ".load.roll", rollWin, sep="")
				src_col <- paste0(target, ".load", sep="")
			}
			df_tmp <- df_facility %>% group_by(day) %>% filter(day >= anchor_date) %>% select(c(day, src_col))
			#print(df_tmp)
			mean4weeks <- mean(df_tmp[[src_col]])
			ci4weeks <- ci99(df_tmp[[src_col]])
			print(mean4weeks)

			gplot <- gplot + geom_point(aes(x = day, y = .data[[rcol]]), color = target_pal(target), shape = 1, size = 1, alpha=0.5) + 
							 				 geom_line(aes(x = day, y = .data[[rcol]]), color = target_pal(target)) + 
											 geom_hline(yintercept=mean4weeks, linetype="dashed", color="#dddddd", size=0.5) + 
											 geom_ribbon(aes(x=day, y=mean4weeks, ymin=mean4weeks-ci4weeks, ymax=mean4weeks+ci4weeks), alpha=0.2)
			if (ci) {
				cicol <- paste0(rcol, ".ci", sep="")
				gplot <- gplot + geom_ribbon(aes(x=day, y=.data[[rcol]], ymin=.data[[rcol]]-.data[[cicol]], ymax=.data[[rcol]]+.data[[cicol]]), alpha=0.2)
			}
		}
		ggplotly(gplot)
	}
	

	#
	# Plot of daily flow
	#
	plotFlow = function(layer, facility) {
		#print("plotFlow called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		
#		print(layer)
#		print(facility)

		if (facility == "State of West Virginia") {
			#capacity <- sum(unique((df_watch %>% filter(group == layer))$capacity_mgd))
			df_plot <- df_watch %>% filter(group == layer) %>% group_by(week_starting) %>% 
				summarize(mean_flow = mean(daily_flow), 
									se_flow = sd(daily_flow) / sqrt(n())
				)
		} else {
			capacity <- unique((df_watch %>% filter(location_common_name == facility))$capacity_mgd)
			df_plot <- df_watch %>% filter(location_common_name == facility) %>% group_by(week_starting) %>% 
				summarize(mean_flow = mean(daily_flow), 
									se_flow = sd(daily_flow) / sqrt(n())
				)
		}

		lims_x_date <- as.Date(strptime(c(first_day, last_day), format = "%Y-%m-%d"))

		gplot <- ggplot(df_plot, aes(x = week_starting, y = mean_flow)) + labs(y = "Weekly mean flow (MGD)", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = "1 month", labels = format_dates, limits = lims_x_date) + 
											geom_errorbar(aes(ymin=mean_flow-se_flow, ymax=mean_flow+se_flow), color="#cccccc", width=.1, alpha=0.5, position=position_dodge(0.05)) + 
											geom_line(color="#aaaaaa", alpha=0.7) + 
											geom_point(aes(color=mean_flow), size=2) + 
											scale_color_gradient(low = "#E7C6B9", high = "#E71417") + 
											ggtitle("Daily flow (MGD), averaged per week") + 
											my_theme()
		if (facility != "State of West Virginia") {
			gplot <- gplot + geom_hline(yintercept=capacity, linetype="dashed", color="#dddddd", size=0.5)
		}
		ggplotly(gplot)
	}

	plotCollection = function(layer, facility) {
		#print("plotCollection called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		
#		print(layer)
#		print(facility)

		if (facility == "State of West Virginia") {
			df_plot <- df_watch %>% filter(group == layer) %>% group_by(week_starting) %>% summarize(count = n())
		} else {
			df_plot <- df_watch %>% filter(location_common_name == facility) %>% group_by(week_starting) %>% summarize(count = n())
		}

		lims_x_date <- as.Date(strptime(c(first_day, last_day), format = "%Y-%m-%d"))

		gplot <- ggplot(df_plot) + labs(y = "Sample count", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = "1 month", labels = format_dates, limits = lims_x_date) + 
											geom_col(aes(x = week_starting, y = count), alpha=0.5) + 
											ggtitle("Number of samples collected per week") + 
											my_theme()
		ggplotly(gplot)
	}
		
	
	# Metadata block reactions to map click or site selection
	md_blockset <- function(layer, facility, rollWin) {
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		if (missing(rollWin)) { rollWin = controlRV$rollWin }
		
		#print(paste0("facility from md_blockset is ", facility, sep=""))
		
		if (facility == "State of West Virginia") {
			df_facility <- df_watch %>% filter(group == layer)
		} else {
			df_facility <- df_watch %>% filter(location_common_name == facility)
		}
				
		num_facilities = n_distinct(df_facility$location_common_name)
		if (num_facilities == 1) {
			facility_text = "facility"
		} else {
			facility_text = "facilities"
		}

		total_cap = sum(distinct(df_facility, location_common_name, capacity_mgd)$capacity_mgd)+1

		total_popserved = sum(distinct(df_facility, location_common_name, population_served)$population_served)

		num_counties = n_distinct(df_facility$counties_served)
		if (num_counties == 1) {
			county_text = "county"
		} else {
			county_text = "counties"
		}
		total_county_pop = sum(distinct(df_facility, counties_served, county_population)$county_population)

		if (total_popserved == -1) {
			total_popserved = "Unknown"
			pct_served <- "(pct unknown)"
		} else {
			pct_served <- paste0(formatC(100*total_popserved/total_county_pop, big.mark=","), "% of ", county_text, sep="")
		}

		total_samples = n_distinct(df_facility$"Sample ID")
		date_first_sampled = min(df_facility$day)
		date_last_sampled = max(df_facility$day)
		
		mean_daily_flow = mean(df_facility$daily_flow)
		
		df_last28days <- df_facility %>% filter(ymd(day) >= ymd(today) - 28)
		mean_collfreq = (n_distinct(df_last28days$"Sample ID"))/4
		
		if (layer == "Sewer Network") {
			layer = "sewer network"
		} else {
			if (facility == "State of West Virginia") {
				layer = "WWTPs"
			}
		}
		
		output$plot_title = renderText(paste0("Showing the ", rollWin, "-day rolling mean for ", facility, " ", layer, sep=""))

		if (layer == "WWTPs") {
			output$data_format <- renderText("Calculated as mean target mass load (copies of target adjusted for average daily flow)")
		} else {
			output$data_format <- renderText("Calculated as mean target copies/L (daily flow is not available at sewers)")
		}

		output$alert_level <- renderText(ALERT_TXT)
		output$site_signal <- renderText(TREND_TXT)

		output$scope <- renderText(facility)
		output$sample_count <- renderText(paste0(formatC(total_samples, big.mark=","), " samples since ", format(as.Date(date_first_sampled), format="%d %b %Y")))

		output$county_population <- renderText(formatC(total_county_pop, big.mark=","))
		output$population_served <- renderText(formatC(total_popserved, big.mark=","))
		output$population_served_pct <- renderText(pct_served)

		output$mean_flow <- renderText(paste0(formatC(mean_daily_flow, big.mark=","), " MGD", sep=""))
		output$collection_frequency <- renderText(paste0(mean_collfreq, " samples/week", sep=""))
		
		update_mod <- "days"
		if (ymd(today)-ymd(date_last_sampled) < 2) {
			update_mod <- "day"
		}
		output$last_update <- renderText(paste0(ymd(today)-ymd(date_last_sampled), " ", update_mod, " ago", sep=""))

		if (facility == "State of West Virginia") {
			output$scope_count <- renderText(paste0(num_facilities, " ", facility_text, sep=""))
			output$counties_served <- renderText(paste0(num_counties, " ", county_text, sep=""))
		} else {
			output$scope_count <- renderText(paste0(num_facilities, " ", facility_text, sep=""))
			output$counties_served <- renderText(paste0(unique(df_facility$counties_served), " county", sep=""))
		}
	}
		
	
	#
	# Init the starting plots and metadata content
	#
	output$watch_plot <- renderPlotly({
		#print("watch_plot top")

		md_blockset()
		#print("md_blockset() just ran!")

		plotLoad() %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		#print("watch_plot bottom")
	})	

	output$focus_plot <- renderPlotly({
		#print("focus_plot top")

		output$focus_plot_title = renderText(paste0("Last 4 weeks of signal from ", controlRV$mapClick, " (", controlRV$activeLayer, ")", sep=""))

		if (controlRV$mapClick == "State of West Virginia") {
			df_facility <- df_watch %>% filter(group == controlRV$activeLayer)
		} else {
			df_facility <- df_watch %>% filter(location_common_name == controlRV$mapClick)
		}
				
		dates <- c(max(df_facility$day) - 30, max(df_facility$day))
		plotLoad(dates = dates) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")

		#print("focus_plot bottom")
	})	

	output$collection_plot <- renderPlotly({
		#print("collection_plot top")

		plotCollections() %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("collection_plot bottom")
	})	

	output$flow_plot <- renderPlotly({
		#print("flow_plot top")

		output$flow_plot_title = renderText(paste0(controlRV$mapClick, " (", controlRV$activeLayer, ")", sep=""))

		plotFlow() %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("flow_plot bottom")
	})	
	
	
	#
	# map interactivity
	#
	
	# Respond to layer change
	observe({
		#print("layer change top")
		selected_group <- req(input$watch_map_groups)
		# only respond if the layer actually changed
		if (selected_group != controlRV$activeLayer) {
			controlRV$activeLayer <- selected_group
			controlRV$mapClick <- "State of West Virginia"
			controlRV$clickLat <- 0
			controlRV$clickLng <- 0
		
			md_blockset(layer = selected_group)
		}
		#print("layer change bottom")
	})


	# Respond to off-marker click
  observeEvent(input$watch_map_click, { 
		#print("Map click event top")
		
		# only respond if this click is in a new position on the map
		if (input$watch_map_click$lat != controlRV$clickLat | input$watch_map_click$lng != controlRV$clickLng) {
			controlRV$clickLat <- 0
			controlRV$clickLng <- 0

			controlRV$mapClick <- "State of West Virginia"
			
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
												 layerId = ~location_common_name, 
												 lat = ~latitude, 
												 lng = ~longitude, 
												 radius = 10, 
												 stroke = TRUE,
												 weight = 2, 
												 color = "dark orange", 
												 fill = TRUE,
												 fillColor = "dark orange",
												 group = "Sewer Network", 
												 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
												 fillOpacity = 0.6)
		
			# Reset the plot
			output$watch_plot <- renderPlotly({
				plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
			})
		
			# Update reactive text elements
			md_blockset()
		#print("Map click event bottom")
			
		}
	}, ignoreInit = TRUE)
	

	# Respond to click on WWTP map marker
  observeEvent(input$watch_map_marker_click, { 
		#print("Map MARKER click top")
    clickedLocation <- input$watch_map_marker_click$id
		controlRV$mapClick <- clickedLocation
		
		controlRV$clickLat <- input$watch_map_marker_click$lat
		controlRV$clickLng <- input$watch_map_marker_click$lng
		
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
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 radius = 10, 
											 stroke = TRUE,
											 weight = 2, 
											 color=~ifelse(location_common_name==clickedLocation, yes = "#73FDFF", no = "dark orange"), 
											 fill = TRUE,
											 fillColor = "dark orange",
											 group = "Sewer Network", 
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.6)

		# Update plots
		output$watch_plot <- renderPlotly({
				plotLoad(facility = clickedLocation) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
		
		# Update reactive text elements
		md_blockset(facility = clickedLocation)
		#print("Map MARKER click bottom")
		
  }, ignoreInit = TRUE)

	
	# Re-center map
	observeEvent(input$center_map, {
		map_center <- MAP_CENTERS %>% filter(layer == controlRV$activeLayer)
		
		watchLeafletProxy %>% setView(map_center$lng, map_center$lat, zoom = map_center$zoom)
#		watchLeafletProxy %>% flyTo(map_center$lng, map_center$lat, zoom = map_center$zoom, options = {animate = TRUE})
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
	
#	observeEvent(input$scope_info_close,{
#    controller(counter_targets=c("scope_info"))
#	}, ignoreInit = TRUE)
	
#	observeEvent(input$population_info_close,{
#    controller(counter_targets=c("population_info"))
#	}, ignoreInit = TRUE)
	
#	observeEvent(input$networkpop_info_close,{
#    controller(counter_targets=c("networkpop_info"))
#	}, ignoreInit = TRUE)
	
	observeEvent(input$daily_flow_info_close,{
    controller(counter_targets=c("daily_flow_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$collection_info_close,{
    controller(counter_targets=c("collection_info"))
	}, ignoreInit = TRUE)
	
	observeEvent(input$last_date_info_close,{
    controller(counter_targets=c("last_date_info"))
	}, ignoreInit = TRUE)
	

	#
	# Respond to change in targets to plot
	#
  observeEvent(input$targets_control, {
		controlRV$visibleTargets <- input$targets_control
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

	#
	# Respond to change in dates to view
	#
	observeEvent(input$dates_control, {
		#print(paste0("Observed: ", input$dates_wwtp, sep = ""))
		#print(paste0("Min: ", min(df_watch$week_starting), ". Max: ", max(df_watch$week_starting), sep=""))
		
#		controlRV$Dates <- format(mdy(input$dates_control), format="%Y-%m-%d")
		controlRV$Dates <- input$dates_control
		#print(controlRV$Dates)
				
		output$watch_plot <- renderPlotly({
			plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)


	#
	# Respond to change in rolling window size
	#
  observeEvent(input$roll_control, {
		controlRV$rollWin <- as.numeric(input$roll_control)
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
		
		layer <- controlRV$activeLayer
		if (layer == "Sewer Network") {
			layer = "sewer network"
		} else {
			if (controlRV$mapClick == "State of West Virginia") {
				layer = "WWTPs"
			}
		}
		output$plot_title = renderText(paste0("Showing the ", controlRV$rollWin, "-day rolling mean for ", controlRV$mapClick, " ", layer, sep=""))

  }, ignoreInit = TRUE)


	#
	# Respond to change in CI visibility
	#
  observeEvent(input$ci_control, {
  	#print(input$ci_control)
		controlRV$ci <- input$ci_control
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

})

