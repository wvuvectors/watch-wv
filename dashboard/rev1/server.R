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
#								Dates = c(min(df_watch$week_starting), max(df_watch$week_starting)),
								Dates = c(first_day, last_day),
								rollWin = SMOOTHER_DEFAULT,
								ci = FALSE,
								activeInfection = INFECTIONS_DEFAULT,
								visibleTargets = TARGETS_DEFAULT,
								activeLayer = "WWTP",
								mapClick = "All facilities",
								clickLat=0, clickLng=0
	)
	
	
	#
	# Determine alert status
	#
	getAlertStatus <- function(facility, layer) {
		
		if (missing(facility)) {
			facility <- controlRV$mapClick
		}
		if (missing(layer)) {
			layer <- controlRV$activeLayer
		}
		
		if (facility == "All facilities") {
			df_alert <- df_signal %>% 
									filter(group == layer) %>% 
									summarize(change_val = mean(current_fold_change_smoothed, na.rm = TRUE),
														trend_val = mean(current_signal_trend, na.rm = TRUE))
		} else {
			df_alert <- df_signal %>% 
									filter(location_common_name == facility) %>% 
									summarize(change_val = current_fold_change_smoothed,
														trend_val = current_signal_trend)
		}
		
		change_bin <- cut(df_alert$change_val, SIGNAL_BINS$fold_change_smoothed, include.lowest=TRUE, labels=FALSE)
		trend_bin <- cut(df_alert$trend_val, SIGNAL_BINS$signal_trend, include.lowest=TRUE, labels=FALSE)
		
#		print(paste0("Change Val     : ", df_alert$change_val, sep=""))
#		print(paste0("Change Bin   : ", change_bin, sep=""))
#		print(paste0("Trend Val     : ", df_alert$trend_val, sep=""))
#		print(paste0("Trend Bin   : ", trend_bin, sep=""))

		color <- SIGNAL_CODES$color[change_bin]
		delta <- SIGNAL_CODES$change[change_bin]
		trend <- SIGNAL_CODES$trend[trend_bin]
		alert <- c(color, delta, trend)
		
		return(alert)

	}
	
	
	#
	# Create the map
	#
	generateMap <- function(center_lat, center_lng, zoom_level) {

		df_map <- df_watch %>% 
							group_by(location_common_name) %>% 
							summarize(longitude = longitude, latitude = latitude, level = level)
		df_map <- left_join(df_map, df_signal, by="location_common_name")
		
		mymap = leaflet(df_map) %>% 
		addTiles() %>% 
		setView(lng = center_lng, lat = center_lat, zoom = zoom_level) %>% 
		addCircleMarkers(data = df_map %>% filter(group == "WWTP"),
										 layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 radius = 10, 
										 stroke = TRUE,
										 weight = 4, 
										 opacity = 0.9,
										 color = ~alertPal(current_fold_change_smoothed), 
										 fill = TRUE,
										 fillColor = ~alertPal(current_fold_change_smoothed), 
#										 fillColor = "black",
										 group = "WWTP", 
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.9) %>%
		addCircleMarkers(data = df_map %>% filter(group == "Sewer Network"),
										 layerId = ~location_common_name, 
										 lat = ~latitude, 
										 lng = ~longitude, 
										 radius = 10, 
										 stroke = TRUE,
										 weight = 4, 
										 opacity = 0.9,
										 color = ~alertPal(current_fold_change_smoothed), 
										 fill = TRUE,
										 fillColor = ~alertPal(current_fold_change_smoothed), 
										 group = "Sewer Network", 
										 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
										 fillOpacity = 0.9) %>% 
#		addProviderTiles(providers$Stadia.AlidadeSmooth) %>%
#		addMiniMap(
#		tiles = providers$Esri.NatGeoWorldMap,
#		toggleDisplay = TRUE)
#		addProviderTiles(providers$Thunderforest.TransportDark)
#							Jawg.Streets
#							Esri.NatGeoWorldMap
#							Esri.WorldTopoMap
#		addPolylines(data=county_sf, fill=FALSE, weight=3, color="#999999", layerId="countiesLayer") %>%
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer") %>%
		addLayersControl(
			position = "bottomleft",
			baseGroups = sort(unique(df_map$group), decreasing=TRUE),
									 options = layersControlOptions(collapsed = FALSE)
		) %>%
		hideGroup(unique(df_map$group)) %>% 
		showGroup("WWTP")
#		fitBounds(~-100,-60,~60,70) %>%
#		addLegend("bottomright", pal = cv_pal, values = ~cv_large_countries$deaths_per_million,
#			title = "<small>Deaths per million</small>")
		
		return(mymap)
	}


	#
	# Plot of pathogen load
	#
	plotLoad <- function(layer, facility, dates, targets) {
		#print("plotLoad called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		if (missing(dates)) { dates = controlRV$Dates }
		if (missing(targets)) { targets = controlRV$visibleTargets }
		
#		print(layer)
#		print(facility)
#		print(dates)
#		print(targets)
		
		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))
		
		if (facility == "All facilities") {
			df_facility <- df_watch %>%
								 		 filter(group == layer & day >= fromDate & day <= toDate) %>%
								 		 group_by(day)
		} else {
			df_facility <- df_watch %>%
								 		 filter(location_common_name == facility & day >= fromDate & day <= toDate) %>%
								 		 group_by(day)
		}
		
		anchor_date <- ymd(max(df_facility$day) - 28)
		
		for (target in targets) {
			if (layer == "Sewer Network") {
				plot_roll <- paste0(target, ".day5.mean", sep="")
				plot_val <- target
				plot_delta <- paste0("fold_change_smoothed", sep="")
			} else {
				plot_roll <- paste0(target, ".loadcap.day5.mean", sep="")
				plot_val <- paste0(target, ".loadcap", sep="")
				plot_delta <- paste0("fold_change_smoothed", sep="")
			}

			if (facility == "All facilities") {
				if (layer == "Sewer Network") {
					df_plot <- df_watch %>% 
										filter(group == layer & day >= fromDate & day <= toDate) %>%
										group_by(day) %>%
										arrange(day) %>%
										summarize(roll := mean(.data[[plot_roll]], na.rm = TRUE),
															val := mean(.data[[plot_val]], na.rm = TRUE),
															delta := mean(.data[[plot_delta]], na.rm = TRUE))
				} else {
					df_plot <- df_watch %>% 
										filter(group == layer & day >= fromDate & day <= toDate) %>%
										group_by(day) %>%
										arrange(day) %>%
										summarize(roll := mean(.data[[plot_roll]], na.rm = TRUE),
															val := mean(.data[[plot_val]], na.rm = TRUE),
															delta := mean(.data[[plot_delta]], na.rm = TRUE))
				}
			} else {
				df_plot <- df_watch %>% 
									filter(location_common_name == facility & day >= fromDate & day <= toDate) %>% 
									group_by(day) %>% 
									arrange(day) %>%
									summarize(roll := mean(.data[[plot_roll]], na.rm = TRUE),
												 		val := mean(.data[[plot_val]], na.rm = TRUE),
												 		delta := mean(.data[[plot_delta]], na.rm = TRUE))
			}
		}
		
		
		lims_x_date <- as.Date(strptime(c(fromDate-1, toDate+1), format = "%Y-%m-%d"))
		
		date_step <- "2 weeks"
		if (toDate - fromDate < 60) {
			date_step <- "2 days"
		}
		if (toDate - fromDate < 15) {
			date_step <- "1 day"
		}

		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = date_step, labels = format_dates, limits = lims_x_date) + 
											scale_color_manual(name = "Target", values = TARGET_COLORS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme()

		for (target in targets) {

			gplot <- gplot + geom_point(aes(x = day, y = val, color = target), shape = 1, size = 2, alpha=0.1) + 
							 				 #geom_point(aes(x = day, y = delta), color = "#991111", shape = 2, size = 2, alpha=0.4) + 
							 				 #geom_line(aes(x = day, y = val), color = target_pal(target), alpha=0.1) + 
							 				 geom_col(aes(x = day, y = val, fill = target), alpha=0.1, na.rm = TRUE) + 
							 				 #geom_point(aes(x = day, y = roll), color = target_pal(target), shape = 1, size = 1, alpha=0.5) + 
							 				 geom_line(aes(x = day, y = roll, color = target))
											 #geom_hline(yintercept=mean_rib, linetype="dashed", color="#dddddd", size=0.5) + 
											 #geom_ribbon(aes(x=day, y=.data[[mean_col]], ymin=.data[[mean_col]]-.data[[ci_col]], ymax=.data[[mean_col]]+.data[[ci_col]]), alpha=0.2)
#			if (ci) {
#				cicol <- paste0(rcol, ".ci", sep="")
#				gplot <- gplot + geom_ribbon(aes(x=day, y=.data[[rcol]], ymin=.data[[rcol]]-.data[[cicol]], ymax=.data[[rcol]]+.data[[cicol]]), alpha=0.2)
#			}
		}
		ggplotly(gplot)# %>% layout(legend = list(orientation = "h"))
	}
	

	#
	# Plot of fold change
	#
	plotFoldChange <- function(layer, facility, dates, targets) {
		#print("plotFoldChange called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		if (missing(dates)) { dates = controlRV$Dates }
		if (missing(targets)) { targets = controlRV$visibleTargets }
		
#		print(layer)
#		print(facility)

		fromDate <- as.Date(ymd(dates[1]))
		toDate <- as.Date(ymd(dates[2]))

		lims_x_date <- as.Date(strptime(c(fromDate-1, toDate+1), format = "%Y-%m-%d"))
		
		date_step <- "2 weeks"
		if (toDate - fromDate < 60) {
			date_step <- "2 days"
		}
		if (toDate - fromDate < 15) {
			date_step <- "1 day"
		}
		
		if (facility == "All facilities") {
			#capacity <- sum(unique((df_watch %>% filter(group == layer))$capacity_mgd))
			df_plot <- df_watch %>% 
								 filter(group == layer & day >= fromDate & day <= toDate) %>% 
								 group_by(day) %>%
								 summarize(fold_change_smoothed = mean(fold_change_smoothed, na.rm = TRUE),
								 					 fold_change = mean(fold_change, na.rm = TRUE))
		} else {
			df_plot <- df_watch %>% 
								 filter(location_common_name == facility & day >= fromDate & day <= toDate) %>% 
								 group_by(day) %>%
								 summarize(fold_change_smoothed = fold_change_smoothed,
								 					 fold_change = fold_change)
		}

		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date(breaks = date_step, labels = format_dates, limits = lims_x_date) + 
											scale_color_manual(name = "Target", values = TARGET_COLORS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme()

		for (target in targets) {

			gplot <- gplot + geom_ribbon(aes(x = day, ymin = SIGNAL_BINS$fold_change_smoothed[1], ymax = SIGNAL_BINS$fold_change_smoothed[2]-1), outline.type = "upper", fill = NA, color = SIGNAL_CODES$color[1], alpha = 0.4) + 
							 				 geom_ribbon(aes(x = day, ymin = SIGNAL_BINS$fold_change_smoothed[2], ymax = SIGNAL_BINS$fold_change_smoothed[3]-1), outline.type = "upper", fill = NA, color = SIGNAL_CODES$color[2], alpha = 0.4) + 
							 				 geom_ribbon(aes(x = day, ymin = SIGNAL_BINS$fold_change_smoothed[3], ymax = SIGNAL_BINS$fold_change_smoothed[4]-1), outline.type = "upper", fill = NA, color = SIGNAL_CODES$color[3], alpha = 0.4) + 
							 				 geom_ribbon(aes(x = day, ymin = SIGNAL_BINS$fold_change_smoothed[4], ymax = SIGNAL_BINS$fold_change_smoothed[5]-1), outline.type = "upper", fill = NA, color = SIGNAL_CODES$color[4], alpha = 0.4) + 
							 				 geom_ribbon(aes(x = day, ymin = SIGNAL_BINS$fold_change_smoothed[5], ymax = max(fold_change, na.rm = TRUE)), outline.type = "upper", fill = NA, color = SIGNAL_CODES$color[5], alpha = 0.4) + 
							 				 geom_col(aes(x = day, y = fold_change, fill = target), alpha=0.4, na.rm = TRUE) + 
							 				 #geom_point(aes(x = day, y = delta), color = "#991111", shape = 2, size = 2, alpha=0.4) + 
							 				 #geom_line(aes(x = day, y = val), color = target_pal(target), alpha=0.1) + 
							 				 #geom_point(aes(x = day, y = fold_change_smoothed, color = target), shape = 1, size = 2, alpha=0.4) + 
							 				 geom_line(aes(x = day, y = fold_change_smoothed, color = target, group = 1))
		}

		ggplotly(gplot)
	}


	#
	# Plot of daily flow
	#
	plotFlow <- function(layer, facility) {
		#print("plotFlow called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		
#		print(layer)
#		print(facility)

		if (facility == "All facilities") {
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
											plot_theme()
		if (facility != "All facilities") {
			gplot <- gplot + geom_hline(yintercept=capacity, linetype="dashed", color="#dddddd", size=0.5)
		}
		ggplotly(gplot)
	}


	#
	# Plot of collection frequency
	#
	plotCollections <- function(layer, facility) {
		#print("plotCollections called!")
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		
#		print(layer)
#		print(facility)

		if (facility == "All facilities") {
			df_plot <- df_watch %>% drop_na(day_received) %>% filter(group == layer) %>% group_by(day_received) %>% summarize(total = n())
#			df_plot <- df_watch %>% filter(group == layer) %>% group_by(day_received)
		} else {
			df_plot <- df_watch %>% drop_na(day_received) %>% filter(location_common_name == facility) %>% group_by(day_received) %>% summarize(total = n())
#			df_plot <- df_watch %>% filter(location_common_name == facility) %>% group_by(day_received)
		}

		df_plot <- df_plot %>% mutate(year = year(day_received),
											 						month = month(day_received, label = TRUE),
																	wkday = fct_relevel(wday(day_received, label=TRUE), c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")),
																	rday = day(day_received),
																	wk = format(day_received, "%W"),
																	total = total) %>%
													select(year, month, wkday, rday, wk, total)
		df_plot <- aggregate(total ~ year + month + wkday, data=df_plot,sum)

#  select(year, month, wkday, day, wk, 5, returns) %>% 
    gplot <- ggplot(df_plot, aes(month, wkday, fill=total)) +
      geom_tile(color="black") +
      #geom_text(aes(label=rday), size=3) + 
      #geom_bar() + 
      labs(x='', y='') + #, title="Test") +
      scale_fill_gradient(low = "#aecf00", high = "#314004", na.value = NA) +
			facet_grid(year~month, scales="free", space="free", drop=FALSE) + 
      plot_theme()
#      facet_grid(year~month, scales="free", space="free")
#		lims_x_date <- as.Date(strptime(c(first_day, last_day), format = "%Y-%m-%d"))
#		gplot <- ggplot(df_plot) + labs(y = "Sample count", x = "") + 
#											scale_y_continuous(labels = comma) + 
#											scale_x_date(breaks = "1 month", labels = format_dates, limits = lims_x_date) + 
#											geom_col(aes(x = day_received, y = count), alpha=0.5) + 
#											ggtitle("Samples received from this facility") + 
#											plot_theme()
		ggplotly(gplot)
	}
		

	plotLoadSpread <- function(layer, facility) {

		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		
		if (layer == "WWTP") {
			df_plot <- df_watch %>% 
								 filter(group == layer) %>% 
								 group_by(location_common_name) %>% 
								 mutate(plot_val = n1n2.loadcap)
			max_y <- 1000
		} else {
			df_plot <- df_watch %>% 
								 filter(group == layer) %>% 
								 group_by(location_common_name) %>% 
								 mutate(plot_val = n1n2)
			max_y <- 5000
		}

		gplot <- ggplot(df_plot) + 
						 geom_boxplot(outlier.shape=NA, aes(x=location_common_name, y=plot_val)) + 
						 scale_y_continuous(limits=c(0,max_y))

		ggplotly(gplot)
	}
	
	#
	# Plot of 2way insight
	#
	plotTwoway <- function(target1, target2) {
		#print("plotTwoway called!")
		
		df_plot <- df_watch %>% arrange(day) %>% 
														mutate(t1 := .data[[target1]],
																	 t2 := .data[[target2]]) %>% 
														select(t1, t2, group)
		df_plot <- df_plot[!is.na(df_plot$t1), ]												
		df_plot <- df_plot[!is.na(df_plot$t2), ]												

    #cc <- cor.test(df_plot$t1, df_plot$t2, method = "pearson") #cc$p.value, cc$estimate
    
    gplot <- ggplot(df_plot, aes(x = t1, y = t2, color = group)) +
      			 geom_point(shape = 1, size = 2, alpha = 0.6) + 
      			 scale_x_log10(labels = comma) + 
      			 scale_y_log10(labels = comma) + 
      			 labs(x="", y="") + 
      			 geom_smooth(formula = y ~ x, method = lm, na.rm = TRUE, se = FALSE, size = 0.25, linetype = 3, color = "#999999") + 
      			 plot_theme()
      			 
		ggplotly(gplot)
	}
		

	#
	# Q-Q plot of 2way insight
	#
	plotTwowayQQ <- function(target) {
		#print("plotTwowayQQ called!")
		
		df_plot <- df_watch %>% arrange(day) %>% 
														mutate(targ := .data[[target]]) %>% 
														select(targ, group)
														    
    gplot <- ggplot(df_plot, aes(sample = targ, color = factor(group))) +
      			 geom_qq(na.rm = TRUE) + 
      			 labs(x="", y="") + 
      			 geom_qq_line(size = 0.25, linetype = 3, color = "#999999") + 
      			 ggtitle(paste0("Q-Q plot of ", toupper(target), sep="")) + 
      			 plot_theme()
      			 
		ggplotly(gplot)
	}
		

	#
	# Correlation coefficient of 2way insight
	#
	calcTwowayCorrelation <- function(target1, target2) {
		df_corr <- df_watch %>% arrange(day) %>% 
														mutate(t1 := .data[[target1]],
																	 t2 := .data[[target2]]) %>% 
														select(t1, t2, group)
		df_corr <- df_corr[!is.na(df_corr$t1), ]												
		df_corr <- df_corr[!is.na(df_corr$t2), ]												
														
    cc <- cor.test(df_corr$t1, df_corr$t2, method = "kendall") #cc$p.value, cc$estimate
    return(cc)
	}
	
	
	#
	# Write the metadata block (reaction to map click or site selection)
	#
	writeMetadata <- function(layer, facility, rollWin, dates) {
		
		if (missing(layer)) { layer = controlRV$activeLayer }
		if (missing(facility)) { facility = controlRV$mapClick }
		if (missing(rollWin)) { rollWin = controlRV$rollWin }
		if (missing(dates)) { dates = controlRV$Dates }
				
		if (facility == "All facilities") {

			df_facility <- df_watch %>% filter(group == layer)
			
			df_alerts <- df_signal %>% 
									 filter(group == layer) %>% 
									 summarize(signal = mean(current_signal_strength_smoothed, na.rm = TRUE),
														 fold_above = mean(current_fold_change_smoothed, na.rm = TRUE),
														 fold_below = mean(current_signal_strength_smoothed/max_signal_strength_smoothed, na.rm = TRUE),
														 trend = mean(current_signal_trend, na.rm = TRUE),
														 scaled_trend = mean(current_scaled_signal_trend, na.rm = TRUE))
										 					 
			num_facilities = n_distinct(df_facility$location_common_name)
			facility_text <- paste0(num_facilities, " ", layer, "s", sep="")

		} else {

			df_facility <- df_watch %>% filter(location_common_name == facility)

			df_alerts <- df_signal %>% 
									 filter(location_common_name == facility) %>% 
									 summarize(signal = current_signal_strength_smoothed,
														 fold_above = current_fold_change_smoothed,
														 fold_below = current_signal_strength_smoothed/max_signal_strength_smoothed,
														 trend = current_signal_trend,
														 scaled_trend = current_scaled_signal_trend)

			num_facilities = 1
			facility_text <- facility

		}
		
		print(paste0("Facility Name   : ", facility, sep=""))
		print(paste0("Signal Strength : ", df_alerts$signal, sep=""))
		print(paste0("Fold Above     : ", df_alerts$fold_above, sep=""))
		print(paste0("Fold Below     : ", df_alerts$fold_below, sep=""))
		print(paste0("Trend           : ", df_alerts$trend, sep=""))
		print(paste0("Scaled Trend    : ", df_alerts$scaled_trend, sep=""))
		
		# get the text code for the fold-change value
		#fold_change_bin <- cut(df_alerts$fold_change, SIGNAL_BINS$fold_change_smoothed, include.lowest=TRUE, labels=FALSE)
		#fold_change_txt <- SIGNAL_CODES$change[fold_change_bin]

		# get the text code for the trend value
		#trend_bin <- cut(df_alerts$trend, SIGNAL_BINS$signal_trend, include.lowest=TRUE, labels=FALSE)
		#trend_txt <- SIGNAL_CODES$trend[trend_bin]
		
		total_cap = sum(distinct(df_facility, location_common_name, capacity_mgd)$capacity_mgd)+1

		num_counties = n_distinct(df_facility$counties_served)
		if (num_counties == 1) {
			county_text = "county"
		} else {
			county_text = "counties"
		}

		total_county_pop = sum(distinct(df_facility, counties_served, county_population)$county_population)

		total_popserved = sum(distinct(df_facility, location_common_name, population_served)$population_served)
		if (total_popserved == -1) {
			total_popserved = "Unknown"
		}

		total_samples = n_distinct(df_facility$"Sample ID")
		date_first_sampled = min(df_facility$day)
		date_last_sampled = max(df_facility$day)
		
		mean_daily_flow = mean(df_facility$daily_flow)
		
		if (layer == "Sewer Network") {
			layer_txt <- "sewer network"
		} else {
			if (facility == "All facilities") {
				layer_txt <- "WWTPs"
			}
		}
		
		# Plot title and legend
		output$plot_title = renderText(paste0("Levels of wastewater ", controlRV$activeInfection, " since July 2021 (", controlRV$mapClick, ")", sep=""))

		if (controlRV$activeLayer != "Sewer Network") {
			output$data_format <- renderText("Calculated as copies of target per person and adjusted for average daily flow. Columns and circles show the daily signal; the solid line shows the 5-day rolling mean.")
		} else {
			output$data_format <- renderText("Calculated as copies of target per liter of wastewater. Columns and circles show the daily signal; the solid line shows the 5-day rolling mean.")
		}

		# Get alert status
		alert_codes <- getAlertStatus()	# 1=color, 2=delta, 3=trend

		output$site_change_hdr <- renderText(paste0("The level of ", controlRV$activeInfection, " is", sep=""))
		#output$site_change_txt <- renderText(paste0(alert_codes[2], sep=""))
		output$site_change_above <- renderText(paste0(prettyNum(df_alerts$fold_above, scientific=FALSE, big.mark=",", digits=2), "X above the lowest value", sep=""))
		output$site_change_below <- renderText(paste0("and ", prettyNum(df_alerts$fold_below*100, scientific=FALSE, big.mark=",", digits=2), "% of the highest value", sep=""))

		output$site_trend_hdr <- renderText(paste0("Rate of ", controlRV$activeInfection, " change:", sep=""))
		output$site_trend_txt <- renderText(paste0(prettyNum(df_alerts$trend, scientific=FALSE, big.mark=",", digits=3), sep=""))
		output$site_trend_level <- renderText(paste0("over the last ", SIGNAL_TREND_WINDOW, " sampling days", sep=""))

		output$alert_hdr <- renderText(paste0("Community spread of ", controlRV$activeInfection, " is", sep=""))
		output$alert_change <- renderText(paste0(alert_codes[2], sep=""))
		output$alert_trend <- renderText(paste0(" and trending ", alert_codes[3], sep=""))
		bgchange <- paste0("document.getElementById('alert_panel').style.backgroundColor = '", alert_codes[1], "';", sep="")
		runjs(bgchange)

		# Report scope, population served, and total number of samples from this facility
		output$scope <- renderText(facility_text)
		output$population_served <- renderText(paste0("Serving ", formatC(total_popserved, big.mark=","), " residents", sep=""))
		output$sample_count <- renderText(paste0(formatC(total_samples, big.mark=","), " samples since ", format(as.Date(date_first_sampled), format="%d %b %Y")))

		# Report county of facility, or total counties served, and county population estimates
		if (num_counties > 1) {
			output$counties_served <- renderText(paste0(num_counties, " ", county_text, sep=""))
		} else {
			output$counties_served <- renderText(paste0(unique(df_facility$counties_served), sep=""))
		}
		output$county_population <- renderText(formatC(total_county_pop, big.mark=","))

		# Report mean daily flow
		output$mean_flow <- renderText(paste0(formatC(mean_daily_flow, big.mark=","), " MGD", sep=""))

		# Report weekly collection frequency
		df_last28days <- df_facility %>% filter(ymd(day) >= ymd(today) - 28)
		mean_collfreq = (n_distinct(df_last28days$"Sample ID"))/4
		output$collection_frequency <- renderText(paste0("Receiving ", mean_collfreq, " samples/week over the last 4 weeks.", sep=""))
		
		# Report the most recent day that data is included for the current location(s)
		update_mod <- "days"
		if (ymd(today)-ymd(date_last_sampled) < 2) {
			update_mod <- "day"
		}
		output$last_update <- renderText(paste0(ymd(today)-ymd(date_last_sampled), " ", update_mod, " ago", sep=""))
		output$last_update_stamp <- renderText(paste0("(Collected on ", format(as.Date(date_last_sampled), format="%d %b %Y"), ")", sep=""))
		
	}

		
	#
	# Toggle visibility of panels
	#
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


	#
	# Render the map
	#
	output$watch_map <- renderLeaflet({
		center <- MAP_CENTERS %>% filter(layer == "WWTP")
		generateMap(center_lat = center$lat, center_lng = center$lng, zoom_level = center$zoom)
	})


	#
	# Render the main plot and update associated metadata
	#
	output$watch_plot <- renderPlotly({
		#print("watch_plot top")

		writeMetadata()
		#print("writeMetadata() just ran!")

		plotLoad() %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		#print("watch_plot bottom")
	})	


	#
	# Render the focus plot and update associated metadata
	#
	output$focus_plot <- renderPlotly({
		#print("focus_plot top")

		output$focus_plot_title = renderText(paste0("Levels of wastewater ", controlRV$activeInfection, " over the last month (", controlRV$mapClick, ")", sep=""))

		if (controlRV$activeLayer != "Sewer Network") {
			output$focus_data_format <- renderText("Calculated as copies of target per person, and adjusted for average daily flow. Columns and circles show the daily signal; the solid line shows the 5-day rolling mean. Rate of change is the slope of the best-fit line over the most recent 5 sampling days.")
		} else {
			output$focus_data_format <- renderText("Calculated as copies of target per liter of wastewater. Columns and circles show the daily signal; the solid line shows the 5-day rolling mean. Rate of change is the slope of the best-fit line over the most recent 5 sampling days.")
		}

		if (controlRV$mapClick == "All facilities") {
			df_facility <- df_watch %>% filter(group == controlRV$activeLayer)
		} else {
			df_facility <- df_watch %>% filter(location_common_name == controlRV$mapClick)
		}
				
		dates <- c(max(df_facility$day) - 30, max(df_watch$day))
		plotLoad(dates = dates) %>% config(displayModeBar = FALSE) #%>% style(hoverinfo = "skip")

		#print("focus_plot bottom")
	})	


	#
	# Render the fold change plot
	#
	output$change_plot <- renderPlotly({
		#print("trend_plot top")

		output$change_plot_title = renderText(paste0("Fold change in ", controlRV$activeInfection, " over the last month (", controlRV$mapClick, ")", sep=""))

		if (controlRV$mapClick == "All facilities") {
			df_facility <- df_watch %>% filter(group == controlRV$activeLayer)
		} else {
			df_facility <- df_watch %>% filter(location_common_name == controlRV$mapClick)
		}
		dates <- c(max(df_facility$day) - 30, max(df_facility$day))

		plotFoldChange(dates = dates) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("trend_plot bottom")
	})	


	#
	# Render the collections plot
	#
	output$collection_plot <- renderPlotly({
		#print("collection_plot top")

		plotCollections() %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("collection_plot bottom")
	})	


	#
	# Render the daily flow plot and associated metadata
	#
	output$flow_plot <- renderPlotly({
		#print("flow_plot top")

		output$flow_plot_title = renderText(paste0(controlRV$mapClick, " (", controlRV$activeLayer, ")", sep=""))

		plotFlow() %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("flow_plot bottom")
	})	
	
	
	#
	# Render the insights-2way plot and associated metadata
	#
	output$insights_plot_2way <- renderPlotly({
		#print("insights_plot_2way top")
		
		target1 <- "n1"
		target2 <- "n2"

		output$plot_2way_title = renderText(paste0(controlRV$activeInfection, " ", toupper(target1), " vs. ", toupper(target2), " (log-transformed)", sep=""))

		cc <- calcTwowayCorrelation(target1, target2)		
#		output$insights_2way_cor = renderText(paste0("The Kendall correlation coefficient between ", controlRV$activeInfection, " ", toupper(target1), " and ", toupper(target2), " is ", prettyNum(cc$estimate, digits = 4), ", with a p-value of ", prettyNum(cc$p.value, digits = 3, format = "e"), "). These data are not normally distributed.", sep=""))
		output$insights_2way_cor = renderText(paste0("The Kendall correlation coefficient between ", controlRV$activeInfection, " ", toupper(target1), " and ", toupper(target2), " is ", prettyNum(cc$estimate, digits = 4), ". These data are not normally distributed.", sep=""))

		plotTwoway(target1="n1", target2="n2") %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("insights_plot_2way bottom")
	})	
	
	#
	# Render the insights-2way plot for target 1
	#
	output$insights_plot_2wayQQ_1 <- renderPlotly({
		#print("insights_plot_2way top")
		target <- "n1"
		
		plotTwowayQQ(target=target) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("insights_plot_2way bottom")
	})	
	

	#
	# Render the insights-2way plot for target 2
	#
	output$insights_plot_2wayQQ_2 <- renderPlotly({
		#print("insights_plot_2way top")
		target <- "n2"
		
#		output$plot_2wayQQ_title = renderText(paste0("Quartile (Q-Q) plot of ", target1, " vs. ", target2, sep=""))

		plotTwowayQQ(target=target) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
		#print("insights_plot_2way bottom")
		
	})	



	###########################
	#
	# MAP EVENTS
	#
	###########################


	#
	# React to layer change on map
	#
	observe({
		#print("layer change top")
		selected_group <- req(input$watch_map_groups)
		# only respond if the layer actually changed
		if (selected_group != controlRV$activeLayer) {
			controlRV$activeLayer <- selected_group
			controlRV$mapClick <- "All facilities"
			controlRV$clickLat <- 0
			controlRV$clickLng <- 0
		
			writeMetadata(layer = selected_group)
		}
		#print("layer change bottom")
	})


	#
	# React to off-marker click on map
	#
  observeEvent(input$watch_map_click, { 
		#print("Map click event top")
		
		# only respond if this click is in a new position on the map
		if (input$watch_map_click$lat != controlRV$clickLat | input$watch_map_click$lng != controlRV$clickLng) {
			controlRV$clickLat <- 0
			controlRV$clickLng <- 0

			controlRV$mapClick <- "All facilities"
			
			df_map <- df_watch %>% 
								group_by(location_common_name) %>% 
								summarize(longitude = longitude, latitude = latitude, level = level)
			df_map <- left_join(df_map, df_signal, by="location_common_name")
			
			# Update map marker colors
			watchLeafletProxy %>% 
				clearMarkers() %>% 
				addCircleMarkers(data = df_map %>% filter(group == "WWTP"),
												 layerId = ~location_common_name, 
												 lat = ~latitude, 
												 lng = ~longitude, 
												 radius = 10, 
												 stroke = TRUE,
												 weight = 4, 
												 opacity = 0.9,
												 color = ~alertPal(current_fold_change_smoothed), 
												 fill = TRUE,
												 fillColor = ~alertPal(current_fold_change_smoothed), 
												 group = "WWTP", 
												 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
												 fillOpacity = 0.9) %>%
				addCircleMarkers(data = df_map %>% filter(group == "Sewer Network"),
												 layerId = ~location_common_name, 
												 lat = ~latitude, 
												 lng = ~longitude, 
												 radius = 10, 
												 stroke = TRUE,
												 weight = 4, 
												 opacity = 0.9,
												 color = ~alertPal(current_fold_change_smoothed), 
												 fill = TRUE,
												 fillColor = ~alertPal(current_fold_change_smoothed), 
												 group = "Sewer Network", 
												 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
												 fillOpacity = 0.9)
		
			# Reset the plot
			output$watch_plot <- renderPlotly({
				plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
			})
		
			# Update reactive text elements
			writeMetadata()
		#print("Map click event bottom")
			
		}
	}, ignoreInit = TRUE)
	
	
	#
	# React to on-marker click on map
	#
  observeEvent(input$watch_map_marker_click, { 
		#print("Map MARKER click top")
    clickedLocation <- input$watch_map_marker_click$id
		controlRV$mapClick <- clickedLocation
		
		controlRV$clickLat <- input$watch_map_marker_click$lat
		controlRV$clickLng <- input$watch_map_marker_click$lng
		
		df_map <- df_watch %>% 
							group_by(location_common_name) %>% 
							summarize(longitude = longitude, latitude = latitude, level = level)
		df_map <- left_join(df_map, df_signal, by="location_common_name")
		
    # Update map marker colors
    watchLeafletProxy %>% 
    	clearMarkers() %>% 
			addCircleMarkers(data = df_map %>% filter(group == "WWTP"),
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 radius = 10, 
											 stroke = TRUE,
											 weight = 4, 
											 opacity = 0.9,
											 color = ~alertPal(current_fold_change_smoothed),
											 fill = TRUE,
											 fillColor = ~ifelse(location_common_name==clickedLocation, yes = "#73FDFF", no = alertPal(current_fold_change_smoothed)), 
											 group = "WWTP", 
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.9) %>%
			addCircleMarkers(data = df_map %>% filter(group == "Sewer Network"),
											 layerId = ~location_common_name, 
											 lat = ~latitude, 
											 lng = ~longitude, 
											 radius = 10, 
											 stroke = TRUE,
											 weight = 4, 
											 opacity = 0.9,
											 color = ~alertPal(current_fold_change_smoothed),
											 fill = TRUE,
											 fillColor = ~ifelse(location_common_name==clickedLocation, yes = "#73FDFF", no = alertPal(current_fold_change_smoothed)), 
											 group = "Sewer Network", 
											 label = ~as.character(paste0(location_common_name, " (" , level, ")")), 
											 fillOpacity = 0.9)

		# Update plots
		output$watch_plot <- renderPlotly({
				plotLoad(facility = clickedLocation) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
		
		# Update reactive text elements
		writeMetadata(facility = clickedLocation)
		#print("Map MARKER click bottom")
		
  }, ignoreInit = TRUE)

	
	#
	# Re-center the map (layer dependent)
	#
	observeEvent(input$center_map, {
		map_center <- MAP_CENTERS %>% filter(layer == controlRV$activeLayer)
		
		watchLeafletProxy %>% setView(map_center$lng, map_center$lat, zoom = map_center$zoom)
#		watchLeafletProxy %>% flyTo(map_center$lng, map_center$lat, zoom = map_center$zoom, options = {animate = TRUE})
	}, ignoreInit = TRUE)





	###########################
	#
	# METADATA PANEL EVENTS
	#
	###########################


	#
	# React to clicks on the metadata panels
	#
	onevent("click", "site_status_panel", togglePanels(on=c("site_change_info")))
	onevent("click", "site_trend_panel", togglePanels(on=c("site_focus_info")))
	onevent("click", "alert_panel", togglePanels(on=c("alert_level_info")))
	#onevent("click", "scope_panel", togglePanels(on=c("scope_info")))
	#onevent("click", "population_panel", togglePanels(on=c("population_info")))
	onevent("click", "daily_flow_panel", togglePanels(on=c("daily_flow_info")))
	onevent("click", "collection_panel", togglePanels(on=c("collection_info")))
	onevent("click", "last_date_panel", togglePanels(on=c("last_date_info")))

	
	#
	# React to click on the site focus panel close button
	#
	observeEvent(input$site_focus_info_close,{
    togglePanels(off=c("site_focus_info"))
	}, ignoreInit = TRUE)
	
	#
	# React to click on the site change panel close button
	#
	observeEvent(input$site_change_info_close,{
    togglePanels(off=c("site_change_info"))
	}, ignoreInit = TRUE)
	
	#
	# React to click on the alert level panel close button
	#
	observeEvent(input$alert_level_info_close,{
    togglePanels(off=c("alert_level_info"))
	}, ignoreInit = TRUE)
	
	#
	# React to click on the daily flow panel close button
	#
	observeEvent(input$daily_flow_info_close,{
    togglePanels(off=c("daily_flow_info"))
	}, ignoreInit = TRUE)
	
	#
	# React to click on the collection panel close button
	#
	observeEvent(input$collection_info_close,{
    togglePanels(off=c("collection_info"))
	}, ignoreInit = TRUE)
	
	#
	# React to click on the last date panel close button
	#
	observeEvent(input$last_date_info_close,{
    togglePanels(off=c("last_date_info"))
	}, ignoreInit = TRUE)
	



	###########################
	#
	# CONTROL PANEL EVENTS
	#
	###########################


	#
	# React to change in targets to plot
	#
  observeEvent(input$targets_control, {
		controlRV$visibleTargets <- input$targets_control
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)


	#
	# React to change in dates to view
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
	# React to change in rolling window size
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
			if (controlRV$mapClick == "All facilities") {
				layer = "WWTPs"
			}
		}
		output$plot_title = renderText(paste0("Levels of wastewater ", controlRV$activeInfection, " since July 2021 (", controlRV$mapClick, ")", sep=""))

  }, ignoreInit = TRUE)


	#
	# React to change in CI visibility
	#
  observeEvent(input$ci_control, {
  	#print(input$ci_control)
		controlRV$ci <- input$ci_control
		
		output$watch_plot <- renderPlotly({
			plotLoad(facility = controlRV$mapClick) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		})
  }, ignoreInit = TRUE)
	
})

