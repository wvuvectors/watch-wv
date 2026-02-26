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
	# Assign Leaflet proxies.
	#
	leafletProxyCovid <- leafletProxy(mapId="map_covid", session)
	leafletProxyFluA <- leafletProxy(mapId="map_flua", session)
	leafletProxyFluB <- leafletProxy(mapId="map_flub", session)
	leafletProxyRsv <- leafletProxy(mapId="map_rsv", session)
	

	#
	# Get the mapProxy by tab index.
	#
	getMapProxy <- function(indx) {
		if (indx == 1) {
			return(leafletProxyCovid)
		} else if (indx == 2) {
			return(leafletProxyFluA)
		} else if (indx == 3) {
			return(leafletProxyFluB)
		} else if (indx == 4) {
			return(leafletProxyRsv)
		}
		return(leafletProxyCovid)
	}
	
	
	#
	# Initialize reactive values with some defaults.
	#
	controlRV <- reactiveValues(
		activeMapColor = c("COVID", "FluA", "FluB", "RSV"), 
		mapClick = c("WV", "WV", "WV", "WV"), 
		trendLines = c(TRUE, TRUE),
		viewMonths = c(VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY), 
		mapClickLat = c(0, 0, 0, 0),
		mapClickLng = c(0, 0, 0, 0),
		mapIndex = 1,
		clicked_shape_id = NULL
	)
	

	# 
	# Watch for change of tab and update the mapIndex and plot targets. This observer also 
	# fires on init so use INIT = FALSE in all other observer functions.
	#
  observe({ 
    if (tolower(input$nav) == "covid") {
      #print("tab = covid")
      controlRV$mapIndex <- 1
			
			updateAllPlots()
			updateSelectionDetails()
			updateAbundanceFreshness()
			updateAlertBlocks()

    } else if (tolower(input$nav) == "flua") {
      #print("tab = fluA")
      controlRV$mapIndex <- 2
			
			updateAllPlots()
			updateSelectionDetails()
			updateAbundanceFreshness()
			updateAlertBlocks()

    } else if (tolower(input$nav) == "flub") {
      #print("tab = fluB")
      controlRV$mapIndex <- 3

			updateAllPlots()
			updateSelectionDetails()
			updateAbundanceFreshness()
			updateAlertBlocks()

    } else if (tolower(input$nav) == "rsv") {
      #print("tab = rsv")
      controlRV$mapIndex <- 4

			updateAllPlots()
			updateSelectionDetails()
			updateAbundanceFreshness()
			updateAlertBlocks()
    }
  })
	

	# The following three observers keep the map and controls synched across tabs.
	
	# 
	# Update the map colorBy when the map content changes.
	#
	observe({
		updateSelectInput(
			session, 
			inputId = "map_color", 
			selected = controlRV$activeMapColor[controlRV$mapIndex]
		)
	}) 

	# 
	# Update the plot view range when the map content changes.
	#
	observe({
		updateSelectInput(
			session, 
			inputId = "view_range", 
			selected = controlRV$viewMonths[controlRV$mapIndex]
		)
	}) 


	# Generic function to toggle state of shiny panels using shinyjs.
	#
	# Accepts a vector of panels to turn on and a vector of panels to turn off.
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


	# Open data warning window(s).
	#
	fireDataWarnings <- function(target_indx) {
		if (missing(target_indx)) {
			target_indx <- controlRV$mapIndex
		}
		plotSuffix <- tolower(DISEASES[target_indx])
		togglePanels(on=c(paste0("missing_data_", plotSuffix, "_popup", sep="")))
	}
	
	# Close data warning window(s).
	#
	killDataWarnings <- function(target_indx) {
		if (missing(target_indx)) {
			target_indx <- controlRV$mapIndex
		}
		plotSuffix <- tolower(DISEASES[target_indx])
		togglePanels(off=c(paste0("missing_data_", plotSuffix, "_popup", sep="")))
	}
	
	
	getColorByIndex <- function(menuOpt) {
		if (missing(menuOpt)) {
			i <- controlRV$mapIndex
			menuOpt <- controlRV$activeMapColor[i]
		}
		
		if (tolower(menuOpt) == "covid") {
			return(1)
		} else if (tolower(menuOpt) == "flua") {
			return(2)
		} else if (tolower(menuOpt) == "flub") {
			return(3)
		} else if (tolower(menuOpt) == "rsv") {
			return(4)
		} else {
			return(1)
		}
		
	}
	
	
	getAbundanceLevel <- function(loc_id, target_index, yr, wk) {
	  #print("##### getAbundanceLevel called!")
		
		if (missing(loc_id)) {
			loc_id <- controlRV$mapClick[controlRV$mapIndex]
		}

		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}

		df_base <- df_rss %>% filter(location_id == loc_id & target == TARGETS[target_index])
		
		# Get the single row that represents the queried date, or the most recent if no date was 
		# specified. This is stored in df_this.
		if (missing(yr) | missing(wk)) {
			df_this <- df_base %>% filter(epi_date == max(df_base$epi_date))
		} else {
			df_this <- df_base %>% filter(epi_year == yr & epi_week == wk)
		}
		
		level <- 1	# Ensures that even if there is some data snafu, we return a valid set.
		
		# Only set the level if exactly one row was returned, and the level is not missing.
		if (nrow(df_this) == 1 & !is.na(df_this$abundance_level)) {
			level <- as.numeric(df_this$abundance_level)
			# If the data exists but is stale, set the level to 1 anyway.
			if (isStale(df_this$primary_date) == TRUE) {
				level <- 1
			}
		}
		
		df_level <- data.frame(
			level = level, 
			label = ALEVEL_STRINGS[level], 
			color = ALEVEL_COLORS[level], 
			detail = ALEVEL_DESCRIPTIONS[level]
		)
		return(df_level)
	}
	
	
	getTrendLevel <- function(loc_id, target_index, yr, wk) {
	  #print("##### getTrendLevel called!")
		
		if (missing(loc_id)) {
			loc_id <- controlRV$mapClick[controlRV$mapIndex]
		}

		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}

		df_base <- df_rss %>% filter(location_id == loc_id & target == TARGETS[target_index])
		
		# Get the single row that represents the queried date, or the most recent if no date was 
		# specified. This is stored in df_this.
		if (missing(yr) | missing(wk)) {
			df_this <- df_base %>% filter(epi_date == max(df_base$epi_date))
		} else {
			df_this <- df_base %>% filter(epi_year == yr & epi_week == wk)
		}
		
		level <- 1	# Ensures that even if there is some data snafu, we return a valid set.
		
		# Only set the level if exactly one row was returned, and the level is not missing.
		if (nrow(df_this) == 1 & !is.na(df_this$trend_level)) {
			level <- as.numeric(df_this$trend_level)
			# If the data exists but is stale, set the level to 1 anyway.
			if (isStale(df_this$primary_date) == TRUE) {
				level <- 1
			}
		}
		
		df_level <- data.frame(
			label = TLEVEL_STRINGS[level], 
			color = TLEVEL_COLORS[level], 
			detail = TLEVEL_DESCRIPTIONS[level]
		)
		return(df_level)
	}
	
	
	getDominantVariant <- function(loc_id, target_index, yr, wk) {
		# do stuff in here!
		return(c("KP.1", "#ffffff"))
	}


	# Accepts a location id and the index of a target.
	#
	# Return a dataframe of abundance data for the current target at the given location.
	#
	getAbundanceData <- function(loc_id, target_index) {
	  #print("##### getAbundanceData called!")
		
		if (missing(loc_id)) {
			loc_id <- controlRV$mapClick[controlRV$mapIndex]
		}

		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}
		
		# Pretty simple except when different epi weeks have the same primary_date. This occurs 
		# at the end of some years. We take the mean of those events.
		df_this <- df_rss %>% 
							 filter(location_id == loc_id & target == TARGETS[target_index]) %>% 
							 arrange(primary_date) %>% 
							 group_by(primary_date) %>% 
							 summarize(
							 	mean_abundance = mean(mean_abundance, na.rm=TRUE),
							 	rolling_mean = mean(rolling_mean, na.rm = TRUE),
							 	epi_week = last(epi_week)
							 )
		df_this$mean_abundance[is.nan(df_this$mean_abundance)] <- NA
		df_this$rolling_mean[is.nan(df_this$rolling_mean)] <- NA
		
		#df_this$primary_date <- lubridate::ymd(df_this$primary_date)
		
		#print("Content of df_this from getAbundanceData:")
		#View(df_this)
		return(df_this)
	}


	# Accepts a dataframe and an integer representing # of months to plot.
	#
	# Generate a ggplotly object of the abundance data within the date window.
	#
	plotAbundance <- function(df_abundance, months_to_plot, target_index) {
	  #print("##### plotAbundance called!")
    #View(df_abundance)

		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}
		
		if (missing(months_to_plot)) {
			months_to_plot <- controlRV$viewMonths
		}
		
		dlab <- case_when(
			months_to_plot == 1 ~ DATE_LABELS[1],
			months_to_plot == 3 ~ DATE_LABELS[2],
			months_to_plot == 6 ~ DATE_LABELS[3],
			months_to_plot == 12 ~ DATE_LABELS[4],
			months_to_plot == 24 ~ DATE_LABELS[5]
		)

		dbrk <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS[1],
			months_to_plot == 3 ~ DATE_BREAKS[2],
			months_to_plot == 6 ~ DATE_BREAKS[3],
			months_to_plot == 12 ~ DATE_BREAKS[4],
			months_to_plot == 24 ~ DATE_BREAKS[5]
		)

		dbrk_minor <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS_MINOR[1],
			months_to_plot == 3 ~ DATE_BREAKS_MINOR[2],
			months_to_plot == 6 ~ DATE_BREAKS_MINOR[3],
			months_to_plot == 12 ~ DATE_BREAKS_MINOR[4],
			months_to_plot == 24 ~ DATE_BREAKS_MINOR[5]
		)
		
		end_date <- this_week
		date_limits <- c(end_date %m-% months(months_to_plot), end_date)

		most_recent_sample_date <- max(df_abundance$primary_date, na.rm=TRUE)
		gos_dates <- c(most_recent_sample_date-12, most_recent_sample_date+3)
		
		# Calculate a 4-period rolling mean, aligned to the right, filling with NA
		df_abundance$roll4 <- rollmean(df_abundance$mean_abundance, k = 4, fill = NA, align = "right")
		# Calculate the mean signal for each target over the last 3 months and last year.
		trend03 <- calcTrend(df_abundance, 3)
		trend12 <- calcTrend(df_abundance, 12)
		
		df_plot <- df_abundance %>% filter(primary_date >= date_limits[1] & primary_date <= date_limits[2])

		if (nrow(df_plot) > 0) {
			
			top_of_tent_win <- max(df_plot$mean_abundance, na.rm=TRUE) * 1.1
			
			t3_y <- trend03 - (0.25 * trend03)
			t12_y <- trend12 + (0.25 * trend12)
			if (trend03 > trend12) {
				t3_y <- trend03 + (0.25 * trend03)
				t12_y <- trend12 - (0.25 * trend12)
			}
			
			gplot <- ggplot(df_plot) + 
				#labs(y = "", x = "") + 
				scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
				scale_x_date(date_breaks = dbrk, date_minor_breaks = dbrk_minor, date_labels = dlab) + 
				plot_theme() + 
				labs(x = NULL, y = NULL, color = NULL) + 
				annotate("rect", xmin = gos_dates[1], xmax = gos_dates[2], ymin = 0, ymax = top_of_tent_win, alpha = 0.7, fill = "#dfdfdf") + 
				annotate("text", x = gos_dates[1]-15, y = top_of_tent_win, color = "#474747", size = 2.5, family = "Arial", label = "Data subject to change") + 
				geom_hline(aes(yintercept=trend03), color=TRENDL_03_COLOR, linetype="dotted", alpha=0.8, linewidth=0.8) + 
				annotate("text", x = gos_dates[1]-100, y = t3_y, color=TRENDL_03_COLOR, size = 2.5, family = "Arial", label = "3 month average") + 
				geom_hline(aes(yintercept=trend12), color=TRENDL_12_COLOR, linetype="dotted", alpha=0.6, linewidth=0.8) + 
				annotate("text", x = gos_dates[1]-30, y = t12_y, color=TRENDL_12_COLOR, size = 2.5, family = "Arial", label = "12 month average") + 
				geom_point(aes(x = primary_date, y = mean_abundance, text=paste0("Epi week ", epi_week, " (", prettyNum(mean_abundance, big.mark=",", digits=1), ")")), na.rm = TRUE, shape = 21, color="#333333", fill="#444444", size = 1, alpha=0.4) + 
				geom_area(aes(x = primary_date, y = roll4), outline.type="upper", alpha=0.5, fill = "#B7B1D6", linewidth=0.5)
		} else {
			gplot <- ggplot()
			fireDataWarnings(target_index)
		}
		
		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	# Accepts a location id and the index of a target.
	#
	# Return a dataframe of fold change and trend slope data for the current target at the given location.
	#
	getChangeData <- function(loc_id, target_index) {
	  #print("##### getChangeData called!")
		
		if (missing(loc_id)) {
			loc_id <- controlRV$mapClick[controlRV$mapIndex]
		}

		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}
		
		# Pretty simple except when different epi weeks have the same primary_date. This occurs 
		# at the end of some years. We take the mean of those events.
		df_a <- df_rss %>% 
							 filter(location_id == loc_id & target == TARGETS[target_index]) %>% 
							 arrange(primary_date) %>% 
							 group_by(primary_date) %>% 
							 summarize(
							 	val = mean(100*abundance_fold_change, na.rm=TRUE),
							 	level = last(abundance_level),
							 	epi_week = last(epi_week),
							 	epi_year = last(epi_year)
							 )
		df_a$category <- "fold_change"
		df_a <- df_a %>% mutate(
			ccolor = ALEVEL_COLORS[level], 
			clabel = ALEVEL_STRINGS[level])

		df_b <- df_rss %>% 
							 filter(location_id == loc_id & target == TARGETS[target_index]) %>% 
							 arrange(primary_date) %>% 
							 group_by(primary_date) %>% 
							 summarize(
							 	val = mean(trend_slope, na.rm=TRUE),
							 	level = last(trend_level),
							 	epi_week = last(epi_week),
							 	epi_year = last(epi_year)
							 )
		df_b$category <- "trend"
		df_b <- df_b %>% mutate(
			ccolor = TLEVEL_COLORS[level], 
			clabel = TLEVEL_STRINGS[level])
		
		df_this <- rbind(df_a, df_b)
		df_this$val[is.nan(df_this$val)] <- NA
		
		return(df_this)
	}


	# Accepts a dataframe and an integer representing # of months to plot.
	#
	# Generate a ggplotly object of the fold change and trend slope within the date window.
	#
	plotChangeVals <- function(df_changes, months_to_plot, target_index) {
	  #print("##### plotChangeVals called!")
		#View(df_changes)
		
		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}
		
		if (missing(months_to_plot)) {
			months_to_plot <- controlRV$viewMonths
		}
		
		dlab <- case_when(
			months_to_plot == 1 ~ DATE_LABELS[1],
			months_to_plot == 3 ~ DATE_LABELS[2],
			months_to_plot == 6 ~ DATE_LABELS[3],
			months_to_plot == 12 ~ DATE_LABELS[4],
			months_to_plot == 24 ~ DATE_LABELS[5]
		)

		dbrk <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS[1],
			months_to_plot == 3 ~ DATE_BREAKS[2],
			months_to_plot == 6 ~ DATE_BREAKS[3],
			months_to_plot == 12 ~ DATE_BREAKS[4],
			months_to_plot == 24 ~ DATE_BREAKS[5]
		)

		dbrk_minor <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS_MINOR[1],
			months_to_plot == 3 ~ DATE_BREAKS_MINOR[2],
			months_to_plot == 6 ~ DATE_BREAKS_MINOR[3],
			months_to_plot == 12 ~ DATE_BREAKS_MINOR[4],
			months_to_plot == 24 ~ DATE_BREAKS_MINOR[5]
		)
		
		end_date <- this_week
		date_limits <- c(end_date %m-% months(months_to_plot), end_date)

		most_recent_sample_date <- max(df_changes$val, na.rm=TRUE)
		
		df_plot <- df_changes %>% filter(primary_date >= date_limits[1] & primary_date <= date_limits[2])

		if (nrow(df_plot) > 0) {
			
			gplot <- ggplot(df_plot) + 
				facet_wrap(~category, nrow=2, scales = "free_y") + 
				scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
				scale_x_date(date_breaks = dbrk, date_minor_breaks = dbrk_minor, date_labels = dlab) + 
				scale_color_identity() + 
				scale_fill_identity() + 
				plot_theme() + 
				labs(x = NULL, y = NULL, color = NULL, title="TOP: Percent fold change compared to 3 month average. BOTTOM: Slope of the 4-week trend line.") + 
#				annotate("text", x = date_limits[1]+5, y = top_of_tent_win, color = "#474747", size = 2.5, family = "Arial", label = "Data subject to change") + 
				geom_col(aes(x = primary_date, y = val, color = ccolor, fill = ccolor, text=paste0("Week ", epi_week, " (", epi_year ,"): ", clabel, sep="")))
#				geom_point(aes(x = primary_date, y = val, color = ccolor, fill = ccolor, text=paste0("Epi week ", epi_week, sep="")), na.rm = TRUE, shape = 21, size = 1, alpha=0.8)
#				geom_area(aes(x = primary_date, y = roll4), outline.type="upper", alpha=0.3, fill = "#B7B1D6", linewidth=0.5)
		} else {
			gplot <- ggplot()
			fireDataWarnings(target_index)
		}
		
		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	#
	# Update all the plots (usually on reaction to map click or site selection).
	#
  updateAllPlots <- function() {
	  #print("##### updateAllPlots called!")

		killDataWarnings(controlRV$mapIndex)
		
		loc_id <- controlRV$mapClick[controlRV$mapIndex]

		if (loc_id == "WV") {
			plot_location <- "West Virginia"
		} else {
			plot_location <- paste0(loc_id, " county, WV", sep="")
		}
		
		target_index <- controlRV$mapIndex
		
		plotSuffix <- tolower(DISEASES[target_index])

		df_abundance <- getAbundanceData(loc_id, target_index)
		aplot_title <- paste0(DISEASE_LABELS[target_index], " abundance in ", plot_location, " wastewater")
		aplot <- plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex], target_index)

		df_changes <- getChangeData(loc_id, target_index)
		cplot <- plotChangeVals(df_changes, controlRV$viewMonths[controlRV$mapIndex], target_index)
		
		if (target_index == 1) {
			output$aplot_title_covid <- renderText(aplot_title)
			output$aplot_covid <- renderPlotly({aplot %>% config(displayModeBar = FALSE)})
			output$cplot_covid <- renderPlotly({cplot %>% config(displayModeBar = FALSE)})
		} else if (target_index == 2) {
			output$aplot_title_flua <- renderText(aplot_title)
			output$aplot_flua <- renderPlotly({aplot %>% config(displayModeBar = FALSE)})
			output$cplot_flua <- renderPlotly({cplot %>% config(displayModeBar = FALSE)})
		} else if (target_index == 3) {
			output$aplot_title_flub <- renderText(aplot_title)
			output$aplot_flub <- renderPlotly({aplot %>% config(displayModeBar = FALSE)})
			output$cplot_flub <- renderPlotly({cplot %>% config(displayModeBar = FALSE)})
		} else if (target_index == 4) {
			output$aplot_title_rsv <- renderText(aplot_title)
			output$aplot_rsv <- renderPlotly({aplot %>% config(displayModeBar = FALSE)})
			output$cplot_rsv <- renderPlotly({cplot %>% config(displayModeBar = FALSE)})
		}
			
  }

  
	#
	# Update the alert status elements (usually on reaction to map click or site selection).
	#
	updateAlertBlocks <- function() {
	  #print("##### updateAlertBlocks called!")
		
		# Get the current map click
		loc_id <- controlRV$mapClick[controlRV$mapIndex]
		
		this_region <- "West Virginia"
		if (loc_id != "WV") {
			this_region <- paste0(loc_id, " county", sep="")
		}
		
		target_index <- controlRV$mapIndex
		target_name <- TARGETS[target_index]
		disease_name <- DISEASES[target_index]
					
		abundance_head <- paste0(disease_name, " Level", sep="")
		trend_head <- paste0(disease_name, " Trend", sep="")
		
		plotSuffix <- tolower(DISEASES[target_index])
		#print(paste0("target index is ", target_index, " and plot suffix is ", plotSuffix, sep=""))
		
		# Generate the value for the abundance block.
		df_alevel <- getAbundanceLevel(loc_id, target_index)

		# Generate the value for the trend block.
		df_tlevel <- getTrendLevel(loc_id, target_index)

		alert_details_text <- getAlertDetail(disease_name, this_region, df_alevel, df_tlevel)
		
		frunner <- paste0(
			"document.getElementById('abundance_text_", plotSuffix, "').style.backgroundColor = '", df_alevel$color, "';",
			"document.getElementById('trend_text_", plotSuffix, "').style.backgroundColor = '", df_tlevel$color, "';", sep="")
		runjs(frunner)
		
		if (target_index == 1) {
			# Write the text blocks. Be nice to construct these variable names from visPlot but can't figure that out.
			output$abundance_head_covid <- renderText(abundance_head)
			output$abundance_covid <- renderText(df_alevel$label)

			output$trend_head_covid <- renderText(trend_head)
			output$trend_covid <- renderText(df_tlevel$label)

			output$alert_details_covid <- renderText(alert_details_text)
		
		} else if (target_index == 2) {
			output$abundance_head_flua <- renderText(abundance_head)
			output$abundance_flua <- renderText(df_alevel$label)

			output$trend_head_flua <- renderText(trend_head)
			output$trend_flua <- renderText(df_tlevel$label)

			output$alert_details_flua <- renderText(alert_details_text)
		
		} else if (target_index == 3) {
			output$abundance_head_flub <- renderText(abundance_head)
			output$abundance_flub <- renderText(df_alevel$label)

			output$trend_head_flub <- renderText(trend_head)
			output$trend_flub <- renderText(df_tlevel$label)

			output$alert_details_flub <- renderText(alert_details_text)
		
		} else if (target_index == 4) {
			output$abundance_head_rsv <- renderText(abundance_head)
			output$abundance_rsv <- renderText(df_alevel$label)

			output$trend_head_rsv <- renderText(trend_head)
			output$trend_rsv <- renderText(df_tlevel$label)

			output$alert_details_rsv <- renderText(alert_details_text)
		
		}
	}
	
	
	#
	# Update the selection details (usually on reaction to map click or site selection).
	#
	updateSelectionDetails <- function(loc_id) {
	  #print("##### updateSelectionDetails called!")
	  #print(paste0("##### mapClick is ", controlRV$mapClick[controlRV$mapIndex], sep=""))
	  
	  if (missing(loc_id)) {
			loc_id <- controlRV$mapClick[controlRV$mapIndex]	# This is a county name, or WV
		}
			
	  #print(paste0("##### loc_id is ", loc_id, sep=""))

		pop_served <- -1
		pop_location <- -1
		pct_served <- -1
				
		if (loc_id == "WV") {
			pop_served <- sum(df_facilities$location_population_served)
			pop_location <- sum(df_counties$county_population)
			num_counties <- length(COUNTIES)
			num_facilities <- length(FACILITIES)
			county_text <- "counties"

			if (pop_served > 0 & pop_location > 0) {
				pct_served <- 100 * pop_served / pop_location
			}

			selection_details <- paste0(
				"West Virginia is represented by ", num_facilities, " active WaTCH facilities serving approximately ", 
				prettyNum(pop_served, big.mark=","), " residents (", 
				prettyNum(pct_served, digits=1) ,"% of the state) across ", 
				num_counties, " counties.", sep=""
			)

		} else {
			facilities_this <- df_facilities %>% filter(loc_id == location_counties_served)
			#View(facilities_this)
			pop_served <- sum(facilities_this$location_population_served)
			pop_location <- as.numeric((df_counties %>% filter(county_id == loc_id))$county_population)
			num_facilities <- length(unique(facilities_this$location_id))
			region_name <- paste0(loc_id, " county", sep="")

			if (pop_served > 0 & pop_location > 0) {
				pct_served <- 100 * pop_served / pop_location
			}
			if (num_facilities == 1) {
				facility_text <- "facility"
			} else {
				facility_text <- "facilities"
			}

			selection_details <- paste0(
				loc_id, " county is represented by ", num_facilities, " active WaTCH ", facility_text, " serving approximately ", 
				prettyNum(pop_served, big.mark=","), " residents (", 
				prettyNum(pct_served, digits=1) ,"% of the county).", sep=""
			)
		}
		
		target_index <- controlRV$mapIndex

		# Update the UI.
		if (target_index == 1) {
			output$selection_details_covid <- renderText(selection_details)
		} else if (target_index == 2) {
			output$selection_details_flua <- renderText(selection_details)
		} else if (target_index == 3) {
			output$selection_details_flub <- renderText(selection_details)
		} else if (target_index == 4) {
			output$selection_details_rsv <- renderText(selection_details)
		}

	}
	

	#
	# Update the abundance data freshness report to the user.
	#
	updateAbundanceFreshness <- function() {
	  #print("##### updateDataFreshness called!")
	  
	  loc_id <- controlRV$mapClick[controlRV$mapIndex]	# This is a county name, or WV
		this_region <- "West Virginia"
		if (loc_id != "WV") {
			this_region <- paste0(loc_id, " county", sep="")
		}
	  
		target_index <- controlRV$mapIndex
		target_name <- TARGETS[target_index]
					
		# Calculate the freshness of the abundance data.
		#	
		df_this <- df_rss %>% filter(location_id == loc_id & target == target_name)
		df_latest <- df_this %>% filter(epi_date == max(epi_date, na.rm = TRUE))
		
		this_freshness <- paste0(
			"The latest abundance data for ", this_region, " was collected during the week that began ", 
			printy_dates(df_latest$primary_date), " (epi week ", df_latest$epi_week, "). ", sep="")

		# Calculate the completeness of the abundance data.
		#	
		num_facilities <- -1
		all_facilities <- -1
		if (loc_id == "WV") {
			df_comp <- df_all %>% filter(location_id %in% FACILITIES & epi_date == df_latest$epi_date)
			num_facilities <- length(unique(df_comp$location_id))
			all_facilities <- length(FACILITIES)
		} else {
			df_countyfacs <- df_facilities %>% filter(location_counties_served == loc_id)
			df_comp <- df_all %>% filter(location_id %in% df_countyfacs$location_id & epi_date == df_latest$epi_date)
			num_facilities <- length(unique(df_comp$location_id))
			all_facilities <- length(unique(df_countyfacs$location_id))
		}
		
		if (num_facilities > 1) {
			site_suffix <- "sites"
		} else {
			site_suffix <- "site"
		}
		this_completeness <- paste0(
			"It includes samples from ", num_facilities, " reporting ", site_suffix, " (", 
			prettyNum(100*(num_facilities/all_facilities), digits=1), "% of ", this_region, " participants).", 
			sep="")
		
		text2print <- paste0(this_freshness, this_completeness, sep = "")

		# Update the UI.
		if (target_index == 1) {
			output$selection_freshness_covid <- renderText(text2print)
		} else if (target_index == 2) {
			output$selection_freshness_flua <- renderText(text2print)
		} else if (target_index == 3) {
			output$selection_freshness_flub <- renderText(text2print)
		} else if (target_index == 4) {
			output$selection_freshness_rsv <- renderText(text2print)
		}
			
	}

		
	getDownloadTableRSS <- function(targetIndex) {
		df_this <- df_rss %>% 
			filter(target == TARGETS[targetIndex]) %>% 
			select(collection_week = primary_date, location_id, target, mean_abundance)

		return(df_this)
	}
		

	output$download_data_covid <- downloadHandler(
		filename = "watch-wv_COVID.csv",
		content = function(file) {
			write.csv(getDownloadTableRSS(1), file, row.names = FALSE)
		}
	)

	output$download_data_flua <- downloadHandler(
		filename = "watch-wv_FluA.csv",
		content = function(file) {
			write.csv(getDownloadTableRSS(2), file, row.names = FALSE)
		}
	)

	output$download_data_flub <- downloadHandler(
		filename = "watch-wv_FluB.csv",
		content = function(file) {
			write.csv(getDownloadTableRSS(3), file, row.names = FALSE)
		}
	)

	output$download_data_rsv <- downloadHandler(
		filename = "watch-wv_RSV.csv",
		content = function(file) {
			write.csv(getDownloadTableRSS(4), file, row.names = FALSE)
		}
	)


	# Accepts a map (target) index.
	#
	#	Reset the zoom and position of the given map.
	#
	resetMap <- function() {
		#print("resetMap called!")
		mapProxy <- getMapProxy(controlRV$mapIndex)
		mapProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
	}
	

	# Accepts a map marker.
	#
	# Respond to a click on the given map marker in the given map.
	#
	clickMapMarker <- function(clicked) {
    #print("##### clickMapMarker called!")
    if (length(clicked) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- clicked$id
			#loc_id <- unique((df_rss %>% filter(location_id == clickedLocation))$location_id)
    }
		controlRV$mapClick[controlRV$mapIndex] <- clickedLocation    
		
# 		if (clickedLocation %in% df_active_loc$location_id | clickedLocation == "WV") {
# 			# Update the reactive element
# 			controlRV$mapClick[controlRV$mapIndex] <- clickedLocation
# 		} else {
# 			print(paste0("No data for ", clickedLocation))
# 		}

		updateAllPlots()
		updateSelectionDetails()
		updateAbundanceFreshness()
		updateAlertBlocks()

	}


	# Accepts a map shape.
	#
	# Respond to a click on the given map shape in the given map.
	#
	clickMapShape <- function(clicked) {
    print("##### clickMapShape called!")

		if (length(clicked) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- clicked$id
			controlRV$mapClickLat[controlRV$mapIndex] <- clicked$lat
			controlRV$mapClickLng[controlRV$mapIndex] <- clicked$lng
    }
		    
    # Reset style of previously clicked shape, if any.
    prev_id <- controlRV$mapClick[controlRV$mapIndex]
    if (prev_id != "WV") {
      getMapProxy(controlRV$mapIndex) %>%
      	removeShape(prev_id) %>% 
				addPolygons(
					data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex] & county_id == prev_id), 
					layerId = ~county_id,
					fill = TRUE,
					stroke = TRUE, 
					fillOpacity = 0.7, 
					opacity = 1.0,
					weight = 2, 
					fillColor = ~ALEVEL_COLORS[this_alevel], 
					color = ~TLEVEL_COLORS[this_tlevel], 
					group="county",
					label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")")), 
					highlightOptions = highlightOptions(
						weight = 4,
						fillOpacity = 0.9))
    }

    # Highlight the newly clicked shape.
    if (clickedLocation != "WV") {
      getMapProxy(controlRV$mapIndex) %>%
      	removeShape(clickedLocation) %>% 
				addPolygons(
					data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex] & county_id == clickedLocation), 
					layerId = ~county_id,
					fill = TRUE,
					stroke = TRUE, 
					fillOpacity = 1.0, 
					opacity = 1.0,
					weight = 3, 
					fillColor = ~ALEVEL_COLORS[this_alevel], 
					color = "#000000", 
					group="county",
					label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")"))
        )
 
    }
    		
		# Update the tracking variable
		controlRV$mapClick[controlRV$mapIndex] <- clickedLocation    

    updateAllPlots()
		#updateSelectionDetails()
		updateAbundanceFreshness()
		updateAlertBlocks()

	}
	

	# Accepts a map click.
	#
	# Respond to a generic off-marker click in the given map.
	#
	clickMapOffMarker <- function(clicked) {
    #print("##### clickMapOffMarker called!")

		# Only respond if this click is in a new position on the map, or if the mapShape observer 
		# hasn't set the map coords to the new shape.
		if (clicked$lat != controlRV$mapClickLat[controlRV$mapIndex] | clicked$lng != controlRV$mapClickLng[controlRV$mapIndex]) {
			print("New position detected!")
			controlRV$mapClickLat[controlRV$mapIndex] <- 0
			controlRV$mapClickLng[controlRV$mapIndex] <- 0

			# Reset style of previously clicked shape, if any.
			prev_id <- controlRV$mapClick[controlRV$mapIndex]
			if (prev_id != "WV") {
				getMapProxy(controlRV$mapIndex) %>%
					removeShape(prev_id) %>% 
					addPolygons(
						data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex] & county_id == prev_id), 
						layerId = ~county_id,
						fill = TRUE,
						stroke = TRUE, 
						fillOpacity = 0.7, 
						opacity = 1.0,
						weight = 2, 
						fillColor = ~ALEVEL_COLORS[this_alevel], 
						color = ~TLEVEL_COLORS[this_tlevel], 
						group="county",
						label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")")), 
						highlightOptions = highlightOptions(
							weight = 3,
							fillOpacity = 0.9))

    	}
		
			# Update the reactive element.
			controlRV$mapClick[controlRV$mapIndex] <- "WV"

			updateAllPlots()
			updateSelectionDetails()
			updateAbundanceFreshness()
			updateAlertBlocks()
		}
	}
	
	
	# Accepts a map zoom level.
	#
	# Respond to a change in the active map zoom level.
	#
	zoomMap <- function(zoom) {
		#print("mapZoom")
		
		mapProxy <- getMapProxy(controlRV$mapIndex)
	}

	
	#
	# Render the maps.
	#
	output$map_covid <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addPolygons(
					data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex]), 
					layerId = ~county_id,
					fill = TRUE,
					stroke = TRUE, 
					fillOpacity = 0.7, 
					opacity = 1.0,
					weight = 2, 
					fillColor = ~ALEVEL_COLORS[this_alevel], 
					color = ~TLEVEL_COLORS[this_tlevel], 
					group="county",
					label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")")), 
					highlightOptions = highlightOptions(
						weight = 3,
						fillOpacity = 0.9))
	})


	output$map_flua <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addPolygons( 
					data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex]), 
					layerId = ~county_id,
					fill = TRUE,
					stroke = TRUE, 
					fillOpacity = 0.7, 
					opacity = 1.0,
					weight = 2, 
					fillColor = ~ALEVEL_COLORS[this_alevel], 
					color = ~TLEVEL_COLORS[this_tlevel], 
					group="county",
					label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")")), 
					highlightOptions = highlightOptions(
						weight = 3,
						fillOpacity = 0.9))
	})


	output$map_flub <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addPolygons( 
					data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex]), 
					layerId = ~county_id,
					fill = TRUE,
					stroke = TRUE, 
					fillOpacity = 0.7, 
					opacity = 1.0,
					weight = 2, 
					fillColor = ~ALEVEL_COLORS[this_alevel], 
					color = ~TLEVEL_COLORS[this_tlevel], 
					group="county",
					label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")")), 
					highlightOptions = highlightOptions(
						weight = 3,
						fillOpacity = 0.9))
	})


	output$map_rsv <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addPolygons( 
					data = map_county_spdf %>% filter(target == TARGETS[controlRV$mapIndex]), 
					layerId = ~county_id,
					fill = TRUE,
					stroke = TRUE, 
					fillOpacity = 0.7, 
					opacity = 1.0,
					weight = 2, 
					fillColor = ~ALEVEL_COLORS[this_alevel], 
					color = ~TLEVEL_COLORS[this_tlevel], 
					group="county",
					label = ~as.character(paste0(county_id, " county (", ALEVEL_STRINGS[this_alevel], " & ", TLEVEL_STRINGS[this_tlevel],")")), 
					highlightOptions = highlightOptions(
						weight = 3,
						fillOpacity = 0.9))
	})


	
	###########################
	#
	# OBSERVER FUNCTIONS - MAP CLICKS
	# These must be specific to each map.
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_covid_marker_click, { 
    #print("##### Map MARKER click top")
    #clickMapMarker(input$map_covid_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flua_marker_click, { 
    #print("##### Map MARKER click top")
    #clickMapMarker(input$map_flua_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flub_marker_click, { 
    #print("##### Map MARKER click top")
    #clickMapMarker(input$map_flub_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_rsv_marker_click, { 
    #print("##### Map MARKER click top")
    #clickMapMarker(input$map_rsv_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)


	# 
	# React to map shape click
	#
  observeEvent(input$map_covid_shape_click, {
    #print("##### Map Shape Click!")
    clickMapShape(input$map_covid_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flua_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_flua_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flub_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_flub_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_rsv_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_rsv_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)



	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	#
  observeEvent(input$map_covid_click, { 
		#print("##### Map off-target click top")
    clickMapOffMarker(input$map_covid_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flua_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_flua_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flub_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_flub_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_rsv_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_rsv_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)


	#
	# React to a change in the map zoom level.
	#
#	observeEvent(input$map_covid_zoom, {
#		print("##### Map zoom top")
#		zoomMap(input$map_covid_zoom)
#	}, ignoreNULL = FALSE, ignoreInit = TRUE)



	###########################
	#
	# OBSERVER FUNCTIONS - MAP AND PLOT CONTROLS
	# These are generic across tabs.
	#
	###########################

	#
	# Reset the map.
	#
	observeEvent(input$map_reset, {
		#print("Reset map!")
		resetMap()
	}, ignoreInit = TRUE)



	#
	# Change the plot date range to view.
	#
  observeEvent(input$view_range, {

		# Update some reactive elements
		controlRV$viewMonths[controlRV$mapIndex] <- as.numeric(input$view_range)
		
		# Update the plots
    updateAllPlots()

	}, ignoreInit = TRUE)




})







#