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

	# Leaflet proxies
	myWWTPLeafletProxy <- leafletProxy(mapId="map_wwtp", session)
	myUpstreamLeafletProxy <- leafletProxy(mapId="map_upstream", session)
	
	# Init reactive values with some defaults
	controlRV <- reactiveValues(
								Dates = c(min(df_watch$week_starting), max(df_watch$week_starting)),
								rollWin = SMOOTHER_DEFAULT,
								ci = "off",
								visibleTargets = TARGETS_DEFAULT
	)
	
	wwtpRV <- reactiveValues(
								mapClick="State of West Virginia", 
								clickLat=0, clickLng=0
	)
	
	upstreamRV <- reactiveValues(
								mapClick="All Locations", 
								clickLat=0, clickLng=0
	)

	# Render a base map
	generateMap <- function(data_in, center_lat, center_lng, zoom_level) {
		mymap = leaflet(data_in) %>% 
		addTiles() %>% 
		setView(lng = center_lng, lat = center_lat, zoom = zoom_level) %>% 
		addCircleMarkers(layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 radius = 10, 
										 stroke = TRUE,
										 weight = 2, 
										 color = "black", 
										 fill = TRUE,
										 fillColor = "black",
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.6)
		
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
	# WWTP map
	#
	output$map_wwtp <- renderLeaflet({
		generateMap(data_in = df_wwtp, center_lat = 38.938532, center_lng = -81.4222577, zoom_level = 8) %>% 
#		addProviderTiles(providers$Thunderforest.TransportDark)
#							Jawg.Streets
#							Esri.NatGeoWorldMap
#							Esri.WorldTopoMap
#		addPolylines(data=county_sf, fill=FALSE, weight=3, color="#999999", layerId="countiesLayer") %>%
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer")
#		addLayersControl(
#			position = "bottomright",
#			overlayGroups = sort(unique(df_watch$group), decreasing=TRUE),
#											options = layersControlOptions(collapsed = FALSE)
#		)
#		hideGroup(TARGETS) %>% 
#		showGroup(TARGETS_DEFAULT)
#		fitBounds(~-100,-60,~60,70) %>%
#		addLegend("bottomright", pal = cv_pal, values = ~cv_large_countries$deaths_per_million,
#			title = "<small>Deaths per million</small>")
	})

	#
	# WWTP plot
	#
	plotWWTP = function(dates, facility, rollwin, ci, targets) {
#		print(dates)
#		print(facility)
#		print(rollwin)
#		print(ci)
#		print(targets)
		
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
		if (facility == "State of West Virginia") {

			loc_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>% 
								summarize(mean.rnamass = mean(rnamass, na.rm = TRUE))
			loc_zoo <- zoo(loc_df$mean.rnamass, loc_df$day)

			loc1_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>% 
								summarize(mean.rnamass.n1 = mean(rnamass.n1, na.rm = TRUE))
			loc1_zoo <- zoo(loc1_df$mean.rnamass.n1, loc1_df$day)
			
			loc2_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>% 
								summarize(mean.rnamass.n2 = mean(rnamass.n2, na.rm = TRUE))
			loc2_zoo <- zoo(loc2_df$mean.rnamass.n2, loc2_df$day)

		} else {
			loc_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate & location_common_name == facility)
			loc_zoo <- zoo(loc_df$rnamass, loc_df$day)
			loc1_zoo <- zoo(loc_df$rnamass.n1, loc_df$day)
			loc2_zoo <- zoo(loc_df$rnamass.n2, loc_df$day)
		}
		
		# calculate the rolling means
		loc_mean_zoo <- rollmean(loc_zoo, rollwin, fill=NA, align="right")
		loc_mean_df <- fortify(loc_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.rnamass"))

		loc1_mean_zoo <- rollmean(loc1_zoo, rollwin, fill=NA, align="right")
		loc1_mean_df <- fortify(loc1_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.rnamass.n1"))

		loc2_mean_zoo <- rollmean(loc2_zoo, rollwin, fill=NA, align="right")
		loc2_mean_df <- fortify(loc2_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.rnamass.n2"))
		
		# join the mean dataframes
		loc_mean_j1_df <- left_join(loc_mean_df, loc1_mean_df, by = c("day" = "day"))
		loc_mean_jfinal_df <- left_join(loc_mean_j1_df, loc2_mean_df, by = c("day" = "day"))

		# calculate the CIs
		loc_ci_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc_ci_df <- fortify(loc_ci_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci"))

		loc1_ci_zoo <- rollapply(loc1_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc1_ci_df <- fortify(loc1_ci_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci.n1"))

		loc2_ci_zoo <- rollapply(loc2_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc2_ci_df <- fortify(loc2_ci_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci.n2"))
		
		# join the CI dataframes
		loc_ci_j1_df <- left_join(loc_ci_df, loc1_ci_df, by = c("day" = "day"))
		loc_ci_jfinal_df <- left_join(loc_ci_j1_df, loc2_ci_df, by = c("day" = "day"))

		# join all the dataframes to the original (loc_df)
		loc_join1_df <- left_join(loc_mean_jfinal_df, loc_ci_jfinal_df, by = c("day" = "day"))
		loc_rolled_df <- left_join(loc_df, loc_join1_df, by = c("day" = "day"))
#		loc_rolled_df <- select(loc_rolled_df, -c(Series.x, Series.y, Series.x.x, Series.y.y))
		
		lims_x_date <- as.Date(strptime(c(fromDate, toDate), format = "%Y-%m-%d"))

		gplot <- ggplot(loc_rolled_df) + labs(y = "Rolling Mean of RNA Mass Load", x = "") + 
		#									geom_errorbar(aes(ymin=rollmean.rnamass-rollmean.se, ymax=rollmean.rnamass+rollmean.se), width=.2, position=position_dodge(0.05)) + 
		#									scale_y_continuous(limits=c(0,NA)) + 
											scale_y_continuous(labels = comma) + 
#											scale_x_date(breaks = "1 month", date_labels = '%d-%b-%Y', limits = lims_x_date) + 
											scale_x_date(breaks = "2 weeks", labels = format_dates, limits = lims_x_date) + 
											my_theme() 
											#scale_fill_manual(values = reds7) +
											#ggtitle("Weekly Mean COVID, All WV Treatment Facilities") + 
											#theme(axis.text.x = element_text(angle = 60, hjust=0.9))

		if ("n1n2" %in% targets) {
			gplot <- gplot + geom_point(aes(x = day, y = rollmean.rnamass), color="#000000", shape = 1, size = 1, alpha=0.5) + 
							 				 geom_line(aes(x = day, y = rollmean.rnamass), color = "#000000")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=rollmean.rnamass, ymin=rollmean.rnamass-rollmean.ci, ymax=rollmean.rnamass+rollmean.ci), alpha=0.2)
			}
		}
		if ("n1" %in% targets) {
			gplot <- gplot + geom_point(aes(x = day, y = rollmean.rnamass.n1), color="#03A049", shape = 0, size = 1, alpha=0.5) + 
											 geom_line(aes(x = day, y = rollmean.rnamass.n1), color = "#03A049")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=rollmean.rnamass.n1, ymin=rollmean.rnamass.n1-rollmean.ci.n1, ymax=rollmean.rnamass.n1+rollmean.ci.n1), color="#03A049", alpha=0.2)
			}
		}
		if ("n2" %in% targets) {
			gplot <- gplot + geom_point(aes(x = day, y = rollmean.rnamass.n2), color="#9437FF", shape = 2, size = 1, alpha=0.5) + 
											 geom_line(aes(x = day, y = rollmean.rnamass.n2), color = "#9437FF")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=rollmean.rnamass.n2, ymin=rollmean.rnamass.n2-rollmean.ci.n2, ymax=rollmean.rnamass.n2+rollmean.ci.n2), color="#9437FF", alpha=0.2) 
			}
		}
		
		ggplotly(gplot)
	}
	
	md_blockset <- function(facility) {
	
		if (missing(facility)) {
			facility = "State of West Virginia"
		}
	
	
		# Generalizable
		#
		output$alert_level <- renderText(ALERT_TXT)
		output$site_signal <- renderText(TREND_TXT)
		output$scope <- renderText(facility)
	
		if (facility == "State of West Virginia") {
			output$scope_count <- renderText(paste0(FACILITY_TOTAL_WWTP, " facilities"))
			output$sample_count <- renderText(paste0(formatC(SAMPLE_TOTAL_WWTP, big.mark=","), " samples since ", FIRST_DATE_WWTP))

			output$counties_served <- renderText(paste0(CTY_TOTAL_WWTP, " counties"))
			output$county_population <- renderText(formatC(POP_TOTAL_CTY, big.mark=","))

			output$population_served <- renderText(formatC(POP_TOTAL_WWTP, big.mark=","))
			output$population_served_pct <- renderText(paste0(formatC(100*POP_TOTAL_WWTP/POP_TOTAL_CTY, big.mark=","), "% of counties"))

			output$mean_flow <- renderText(paste0(formatC(mean(df_watch$daily_flow, big.mark=",")), " MGD"))

			df_lastmon <- df_wwtp %>% filter(ymd(day) >= ymd(today) - 28)

			output$collection_frequency <- renderText(paste0((n_distinct(df_lastmon$"Sample ID")/4), " samples/week"))
			#output$collection_scheme <- renderText("composite samples")

			output$last_update <- renderText(paste0(ymd(today)-ymd(LAST_DATE_WWTP)-1, " days ago"))

		} else {
			df_facility <- df_wwtp %>% filter(location_common_name == facility)
			SAMPLE_TOTAL = formatC(n_distinct(df_facility$"Sample ID"), big.mark=",")
			FIRST_DATE = format(min(df_facility$day), format = "%b %d, %Y")
			LAST_DATE = max(df_facility$day)
			COUNTY = unique(df_facility$counties_served)
			COUNTY_POP = unique(df_facility$county_population)
			POP_SERVED = unique(df_facility$population_served)
		
			#COLL_TYPE = unique(watch_fac$collection_scheme)

			output$scope_count <- renderText("1 facility")
			output$sample_count <- renderText(paste0(formatC(SAMPLE_TOTAL, big.mark=","), " samples since ", FIRST_DATE))

			output$counties_served <- renderText(paste0(COUNTY, " county"))
			output$county_population <- renderText(formatC(COUNTY_POP, big.mark=","))

			output$population_served <- renderText(formatC(POP_SERVED, big.mark=","))
			output$population_served_pct <- renderText(paste0(formatC(100*POP_SERVED/COUNTY_POP, big.mark=","), "% of county"))
	
			output$mean_flow <- renderText(paste0(formatC(mean(df_facility$daily_flow, big.mark=",")), " MGD"))

			df_lastmon <- df_facility %>% filter(ymd(day) >= ymd(today) - 28)
			output$collection_frequency <- renderText(paste0((n_distinct(df_lastmon$"Sample ID")/4), " samples/week"))
			#output$collection_scheme <- renderText("composite samples")

			output$last_update <- renderText(paste0(ymd(today)-ymd(LAST_DATE)-1, " days ago"))
		}
	}
		
	
	#
	# Init basic elements
	#

	output$plot_wwtp <- renderPlotly({
		plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
	})
	
	md_blockset()
	

	
	#
	# map interactivity
	#
	
	# Respond to off-marker click
  observeEvent(input$map_wwtp_click, { 
		#print("Map click event top")
		
		# only respond if this click is in a new position on the map
		if (input$map_wwtp_click$lat != wwtpRV$clickLat | input$map_wwtp_click$lng != wwtpRV$clickLng) {
			wwtpRV$clickLat <- 0
			wwtpRV$clickLng <- 0

			wwtpRV$mapClick <- "State of West Virginia"
		
			# Update map marker colors
			myWWTPLeafletProxy %>% 
				clearMarkers() %>% 
				addCircleMarkers(data = df_wwtp,
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 radius = 10, 
											 stroke = TRUE,
											 weight = 5, 
											 color = "black", 
											 fill = TRUE,
											 fillColor = "black",
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.6)
		
			# Reset WWTP plots
			output$plot_wwtp <- renderPlotly({
				plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
			})
		
#			output$plot_wwtp_big <- renderPlotly({
#				plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE)
#			})

			# Update reactive text elements
			md_blockset()
			
		}
	}, ignoreInit = TRUE)
	

	# Respond to click on WWTP map marker
  observeEvent(input$map_wwtp_marker_click, { 
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
										 radius = 10, 
										 stroke = TRUE,
										 weight = 5, 
										 color=~ifelse(location_common_name==clickedLocation, yes = "#73FDFF", no = "black"), 
										 fill = TRUE,
										 fillColor = "black",
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.6)
		
		# Update plots
		output$plot_wwtp <- renderPlotly({
				plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
		
		# Update reactive text elements
		md_blockset(clickedLocation)
		
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
  observeEvent(input$targets_wwtp, {
		controlRV$visibleTargets <- input$targets_wwtp
		updatePrettyCheckboxGroup(session = session, inputId = "targets_upstream", selected = controlRV$visibleTargets)
		
		output$plot_wwtp <- renderPlotly({
			plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

  observeEvent(input$targets_upstream, {
		controlRV$visibleTargets <- input$targets_upstream
		updatePrettyCheckboxGroup(session = session, inputId = "targets_wwtp", selected = controlRV$visibleTargets)
		
		output$plot_upstream <- renderPlotly({
			plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

	
	#
	# Respond to change in dates to view
	#
	observeEvent(input$dates_wwtp, {
		#print(paste0("Observed: ", input$dates_wwtp, sep = ""))
		#print(paste0("Min: ", min(df_watch$week_starting), ". Max: ", max(df_watch$week_starting), sep=""))
		
		controlRV$Dates <- input$dates_wwtp
		updateSliderTextInput(session = session, inputId = "dates_upstream", selected = controlRV$Dates)
		
		output$plot_wwtp <- renderPlotly({
			plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

	observeEvent(input$dates_upstream, {
		controlRV$Dates <- input$dates_upstream
		updateSliderTextInput(session = session, inputId = "dates_wwtp", selected = controlRV$Dates)
		
		output$plot_upstream <- renderPlotly({
			plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)


	#
	# Respond to change in rolling window size
	#
  observeEvent(input$roll_wwtp, {
		controlRV$rollWin <- input$roll_wwtp
		updateSliderInput(session = session, inputId = "roll_upstream", value = controlRV$rollWin)
		
		output$plot_wwtp <- renderPlotly({
			plotWWTP(controlRV$Dates, wwtpRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)

  observeEvent(input$roll_upstream, {
		controlRV$rollWin <- input$roll_upstream
		updateSliderInput(session = session, inputId = "roll_wwtp", value = controlRV$rollWin)
		
		output$plot_upstream <- renderPlotly({
			plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)












	#
	# Sewer network default info panel elements
	#
	output$plot_upstream <- renderPlotly({
			plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
	})

	output$last_update_upstream <- renderText(paste0("Sewer Network Wastewater Nowcast (", format(LAST_DATE_UPSTREAM, format = "%d %b %Y"), ")", sep = " "))
	output$all_facilities_upstream <- renderText(paste0(FACILITY_TOTAL_UPSTREAM, " active sewer network sites.", sep = " "))
	output$all_samples_upstream <- renderText(paste0(formatC(SAMPLE_TOTAL_UPSTREAM, big.mark=","), " total samples processed since ", format(FIRST_DATE_UPSTREAM, format = "%d %b %Y"), ".", sep = " "))
	output$facility_name_upstream <- renderText(paste0("All Locations", sep = " "))


	plotUpstream = function(dates, facility, rollwin, ci, targets) {
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
		if (facility == "All Locations") {

			loc_df <- df_upstream %>% filter(week_starting >= fromDate & week_starting <= toDate) %>% 
								group_by(day) %>% 
								summarize(mean.n1n2 = mean(n1n2, na.rm = TRUE))
			loc_zoo <- zoo(loc_df$mean.n1n2, loc_df$day)
			loc1_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>% 
								summarize(mean.n1 = mean(n1, na.rm = TRUE))
			loc1_zoo <- zoo(loc1_df$mean.n1, loc1_df$day)
			
			loc2_df <- df_wwtp %>% filter(week_starting >= fromDate & week_starting <= toDate) %>%
								group_by(day) %>% 
								summarize(mean.n2 = mean(n2, na.rm = TRUE))
			loc2_zoo <- zoo(loc2_df$mean.n2, loc2_df$day)
		} else {
			loc_df <- df_upstream %>% filter(week_starting >= fromDate & week_starting <= toDate & location_common_name == facility)
			loc_zoo <- zoo(loc_df$n1n2, loc_df$day)
			loc1_zoo <- zoo(loc_df$n1, loc_df$day)
			loc2_zoo <- zoo(loc_df$n2, loc_df$day)
		}
		
		# calculate the rolling means
		loc_mean_zoo <- rollmean(loc_zoo, rollwin, fill=NA, align="right")
		loc_mean_df <- fortify(loc_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.n1n2"))

		loc1_mean_zoo <- rollmean(loc1_zoo, rollwin, fill=NA, align="right")
		loc1_mean_df <- fortify(loc1_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.n1"))

		loc2_mean_zoo <- rollmean(loc2_zoo, rollwin, fill=NA, align="right")
		loc2_mean_df <- fortify(loc2_mean_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.n2"))
		
		# join the mean dataframes
		loc_mean_j1_df <- left_join(loc_mean_df, loc1_mean_df, by = c("day" = "day"))
		loc_mean_jfinal_df <- left_join(loc_mean_j1_df, loc2_mean_df, by = c("day" = "day"))

		# calculate the CIs
		loc_ci_zoo <- rollapply(loc_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc_ci_df <- fortify(loc_ci_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci.n1n2"))

		loc1_ci_zoo <- rollapply(loc1_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc1_ci_df <- fortify(loc1_ci_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci.n1"))

		loc2_ci_zoo <- rollapply(loc2_zoo, width=rollwin, fill=NA, align="right", FUN = ci90)
		loc2_ci_df <- fortify(loc2_ci_zoo, melt=TRUE, names=c(Index="day", Value="rollmean.ci.n2"))
		
		# join the CI dataframes
		loc_ci_j1_df <- left_join(loc_ci_df, loc1_ci_df, by = c("day" = "day"))
		loc_ci_jfinal_df <- left_join(loc_ci_j1_df, loc2_ci_df, by = c("day" = "day"))

		# join all the dataframes to the original (loc_df)
		loc_join1_df <- left_join(loc_mean_jfinal_df, loc_ci_jfinal_df, by = c("day" = "day"))
		loc_rolled_df <- left_join(loc_df, loc_join1_df, by = c("day" = "day"))
#		loc_rolled_df <- select(loc_rolled_df, -c(Series.x, Series.y, Series.x.x, Series.y.y))
		
		lims_x_date <- as.Date(strptime(c(fromDate, toDate), format = "%Y-%m-%d"))

		gplot <- ggplot(loc_rolled_df) + labs(y = "Rolling Mean of Copies/L", x = "") + 
		#									geom_errorbar(aes(ymin=rollmean.rnamass-rollmean.se, ymax=rollmean.rnamass+rollmean.se), width=.2, position=position_dodge(0.05)) + 
		#									scale_y_continuous(limits=c(0,NA)) + 
											scale_y_continuous(labels = comma) + 
#											scale_x_date(breaks = "1 month", date_labels = '%d-%b-%Y', limits = lims_x_date) + 
											scale_x_date(breaks = "2 weeks", labels = format_dates, limits = lims_x_date) + 
											my_theme() 
											#scale_fill_manual(values = reds7) +
											#ggtitle("Weekly Mean COVID, All WV Treatment Facilities") + 
											#theme(axis.text.x = element_text(angle = 60, hjust=0.9))

		if ("n1n2" %in% targets) {
			gplot <- gplot + geom_point(aes(x = day, y = rollmean.n1n2), color="#000000", shape = 1, size = 1, alpha=0.5) + 
							 				 geom_line(aes(x = day, y = rollmean.n1n2), color = "#000000")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=rollmean.n1n2, ymin=rollmean.n1n2-rollmean.ci.n1n2, ymax=rollmean.n1n2+rollmean.ci.n1n2), alpha=0.2)
			}
		}
		if ("n1" %in% targets) {
			gplot <- gplot + geom_point(aes(x = day, y = rollmean.n1), color="#03A049", shape = 0, size = 1, alpha=0.5) + 
											 geom_line(aes(x = day, y = rollmean.n1), color = "#03A049")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=rollmean.n1, ymin=rollmean.n1-rollmean.ci.n1, ymax=rollmean.n1+rollmean.ci.n1), color="#03A049", alpha=0.2)
			}
		}
		if ("n2" %in% targets) {
			gplot <- gplot + geom_point(aes(x = day, y = rollmean.n2), color="#9437FF", shape = 2, size = 1, alpha=0.5) + 
											 geom_line(aes(x = day, y = rollmean.n2), color = "#9437FF")
			if (ci != "off") {
				gplot <- gplot + geom_ribbon(aes(x=day, y=rollmean.n2, ymin=rollmean.n2-rollmean.ci.n2, ymax=rollmean.n2+rollmean.ci.n2), color="#9437FF", alpha=0.2) 
			}
		}

		ggplotly(gplot)

	}

	#
	# Sewer network map
	#
	output$map_upstream <- renderLeaflet({
		generateMap(data_in = df_upstream, center_lat = 39.6352701, center_lng = -80.0125177, zoom_level = 13) %>% 
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer")
#		addLayersControl(
#			position = "bottomright",
#			overlayGroups = TARGETS,
#											options = layersControlOptions(collapsed = FALSE)
#		) %>%
#		hideGroup(TARGETS) %>% 
#		showGroup(TARGETS_DEFAULT)
		
	})

	#
	# map interactivity
	#
	
	# Respond to off-marker click
  observeEvent(input$map_upstream_click, { 
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
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.6)
		
			# Update plots
			output$plot_upstream <- renderPlotly({
				plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
			})
			output$plot_upstream_big <- renderPlotly({
				plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE)
			})
		
			# Update reactive text elements
			output$facility_name_upstream <- renderText(paste0("All Locations", sep = " "))
			output$facility_samples_upstream <- renderText("")
		}
	}, ignoreInit = TRUE)


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
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.6)
		
		# Update plots
		output$plot_upstream <- renderPlotly({
			plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
		output$plot_upstream_big <- renderPlotly({
			plotUpstream(controlRV$Dates, upstreamRV$mapClick, controlRV$rollWin, controlRV$ci, controlRV$visibleTargets) %>% config(displayModeBar = FALSE)
		})
		
		# Update reactive text elements
		watch_upstream <- df_upstream %>% filter(location_common_name == clickedLocation)
		SAMPLE_TOTAL_EDULOC = formatC(n_distinct(watch_upstream$"Sample ID"), big.mark=",")
		FIRST_DATE_EDULOC = format(min(watch_upstream$day), format = "%d %b %Y")
		LAST_DATE_EDULOC = format(max(watch_upstream$day), format = "%d %b %Y")
		COLL_TYPE_UP = unique(watch_upstream$collection_scheme)

		output$facility_name_upstream <- renderText(paste0(clickedLocation, sep = " "))
		output$facility_samples_upstream <- renderText(paste0(SAMPLE_TOTAL_EDULOC, " samples collected as ", COLL_TYPE_UP, "s from ", clickedLocation, " processed between ", FIRST_DATE_EDULOC," and ", LAST_DATE_EDULOC, ".", sep = " "))
#		output$facility_capacity <- renderText(paste0("Max capacity of this facility: ", formatC(watch_fac$capacity_mgd[1], big.mark=","), " million gallons per day (MGD).", sep = " "))
#		output$facility_popserved <- renderText(paste0("Estimated population included: ", formatC(watch_fac$population_served[1], big.mark=","), ".", sep = " "))
#		output$facility_counties <- renderText(paste0("This facility is located in ", watch_fac$counties_served[1], " county WV.", sep = " "))
  }, ignoreInit = TRUE)


})

