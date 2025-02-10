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
		activeMapColor = c("Trend", "Trend", "Trend"), 
		mapClick = c("WV", "WV", "WV"), 
		trendLines = c(TRUE, TRUE),
		viewMonths = c(VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY), 
		mapClickLat = c(0, 0, 0),
		mapClickLng = c(0, 0, 0),
		mapIndex = 1,
		targetVec = c(1)
	)
	
	# 
	# Watch for change of tab and update the mapIndex and plot targets
	#
  observe({ 
    if(input$nav == "COVID") {
      #print("tab = covid")
      controlRV$mapIndex <- 1
      controlRV$targetVec <- c(1)

			updateAllPlots()
			updateSelectionInfo()
			updateStatus()

    } else if(input$nav == "Influenza"){
      #print("tab = flu")
      controlRV$mapIndex <- 2
      controlRV$targetVec <- c(2,3)

			updateAllPlots()
			updateSelectionInfo()
			updateStatus()
    } else if(input$nav == "RSV"){
      #print("tab = rsv")
      controlRV$mapIndex <- 3
      controlRV$targetVec <- c(4)

			updateAllPlots()
			updateSelectionInfo()
			updateStatus()
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
	# Update the map geoLevel when the map content changes.
	#
	observe({
		updateSelectInput(
			session, 
			inputId = "geo_level", 
			selected = controlRV$activeGeoLevel[controlRV$mapIndex]
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


	getAbundance <- function(loc_id, target_index) {
	  print("##### getAbundance called!")

		val <- (dflist_alerts[[target_index]] %>% filter(region_name == loc_id))$abundance_pct_change
		
		if (is.na(val)) {
			return(c(ALERT_LEVEL_STRINGS[5], ALERT_LEVEL_COLORS[5]))
		}
		
		if (val < ALERT_LEVEL_THRESHOLDS[1]) {
			# LOW
			txt <- ALERT_LEVEL_STRINGS[1]
			color <- ALERT_LEVEL_COLORS[1]
		} else if (val >= ALERT_LEVEL_THRESHOLDS[1] & val < ALERT_LEVEL_THRESHOLDS[2]) {
			# MODERATE
			txt <- ALERT_LEVEL_STRINGS[2]
			color <- ALERT_LEVEL_COLORS[2]
		} else if (val >= ALERT_LEVEL_THRESHOLDS[2] & val < ALERT_LEVEL_THRESHOLDS[3]) {
			# HIGH
			txt <- ALERT_LEVEL_STRINGS[3]
			color <- ALERT_LEVEL_COLORS[3]
		} else if (val >= ALERT_LEVEL_THRESHOLDS[3]) {
			# VERY HIGH
			txt <- ALERT_LEVEL_STRINGS[4]
			color <- ALERT_LEVEL_COLORS[4]
		} else {
			# UNKNOWN
			txt <- ALERT_LEVEL_STRINGS[5]
			color <- ALERT_LEVEL_COLORS[5]
		}
		
		return(c(txt, color))
	}
	
	
	getTrend <- function(loc_id, target_index) {
	# 	i <- 3
	# 	txt <- ALERT_LEVEL_STRINGS[i]
	# 	color <- ALERT_LEVEL_COLORS[i]
		print((dflist_alerts[[target_index]] %>% filter(region_name == loc_id))$trend)
		
		txt <- (dflist_alerts[[target_index]] %>% filter(region_name == loc_id))$trend
		color <- ALERT_LEVEL_COLORS[5]

		for (i in 1:length(TREND_STRINGS)) {
			tstr <- TREND_STRINGS[i]
			if (txt == tstr) {
				color <- ALERT_LEVEL_COLORS[i]
				break
			}
		}
		
		return(c(txt, color))
	}
	
	
	getDominantVariant <- function(loc_id, target_index) {
		i <- 4
		txt <- ALERT_LEVEL_STRINGS[i]
		color <- ALERT_LEVEL_COLORS[i]
		# do stuff in here!
	
		return(c("KP.1", "#ffffff"))
	}

	# Accepts a location id and the index of the target
	#
	# Return a dataframe of abundance data for the current bio target at the given location.
	#
	getAbundanceData <- function(loc_id, target_index) {
	  print("##### getAbundanceData called!")
	  #print(controlRV$mapIndex)
		
		if (missing(loc_id)) {
			loc_id <- controlRV$mapClick[controlRV$mapIndex]
		}

		if (missing(target_index)) {
			target_index <- controlRV$mapIndex
		}

		df_base <- dflist_rs[[target_index]]
		#View(df_base)
		
		if (loc_id == "WV") {
			# get all active facilities
			loc_ids <- (df_active_loc %>% filter(location_category == "wwtp"))$location_id
		} else {
			if ((df_regions %>% filter(region_name == loc_id))$region_geolevel == "county") {
				# roll up county active facilities
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			} else {
				# just the single facility, but need it as a vector
				loc_ids <- (df_active_loc %>% filter(location_id == loc_id))$location_id
			}
		}

		df_this <- df_base %>% 
							 filter(location_id %in% loc_ids) %>%
							 group_by(date_to_plot) %>% 
							 arrange(date_to_plot) %>%
							 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE),
												date_primary := date_primary)
		#View(df_this)
		return(df_this)
	}


	# Accepts a location id.
	#
	# Return a dataframe of variant data for the current bio target at the given location.
	#
	getVariantData <- function(loc_id) {
	  print("##### getVariantData called!")
	
		if (loc_id == "WV") {
			# get all active facilities
			loc_ids <- (df_active_loc %>% filter(location_category == "wwtp"))$location_id
		} else {
			if ((df_regions %>% filter(region_name == loc_id))$region_geolevel == "county") {
				# roll up county active facilities
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			} else {
				# just the single facility, but need it as a vector
				loc_ids <- (df_active_loc %>% filter(location_id == loc_id))$location_id
			}
		}
		
		df_trans <- df_seqr %>% filter(location_id %in% loc_ids)

		t1 <- df_trans %>% group_by(date_to_plot, color_group) %>% tally(variant_proportion)	# n = sum of variant prop across all samples in a date
		t2 <- df_trans %>% group_by(date_to_plot) %>% mutate(location_count = n_distinct(location_id, na.rm = TRUE)) %>% select(date_to_plot, location_count) %>% distinct()
		t3 <- df_trans %>% group_by(date_to_plot) %>% mutate(collection_count = n_distinct(sample_collection_end_datetime, na.rm = TRUE)) %>% select(date_to_plot, collection_count) %>% distinct()
		t23 <- merge(t2, t3, by = "date_to_plot")
		df_this <- merge(t1, t23, by = "date_to_plot")
		df_this$total_prop <- (df_this$n / df_this$location_count) / df_this$collection_count
		df_this$total_pct <- as.numeric(formatC(100*df_this$total_prop, format="f", digits=2))

		return(df_this)

	}


	# Accepts a dataframe and an integer representing # of months to plot.
	#
	# Generate a ggplotly object of the abundance data within the date window.
	#
	plotAbundance <- function(df_abundance, months_to_plot) {
	  print("##### plotAbundance called!")
    #View(df_abundance)

		dlab <- case_when(
			months_to_plot == 1 ~ DATE_LABELS[1],
			months_to_plot == 3 ~ DATE_LABELS[2],
			months_to_plot == 6 ~ DATE_LABELS[3],
			months_to_plot == 12 ~ DATE_LABELS[4],
			months_to_plot == 24 ~ DATE_LABELS[5],
		)

		dbrk <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS[1],
			months_to_plot == 3 ~ DATE_BREAKS[2],
			months_to_plot == 6 ~ DATE_BREAKS[3],
			months_to_plot == 12 ~ DATE_BREAKS[4],
			months_to_plot == 24 ~ DATE_BREAKS[5],
		)
		
		end_date <- this_week
		date_limits <- c(end_date %m-% months(months_to_plot), end_date)

		most_recent_sample_date <- max(df_abundance$date_to_plot, na.rm=TRUE)
		gos_dates <- c(most_recent_sample_date-12, most_recent_sample_date+3)
		
		df_plot <- df_abundance %>% filter(date_to_plot >= date_limits[1] & date_to_plot <= date_limits[2])
		largest_val <- max(df_plot$val+10, na.rm=TRUE)
		
		# Calculate the mean signal for each target over the last 3 months and last year.
		trend03 <- calcTrend(df_abundance, 3)
		trend12 <- calcTrend(df_abundance, 12)

#		geompt_msg <- paste0(date_to_plot, sep="")
		
		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
											scale_x_date(date_breaks = dbrk, date_labels = dlab, limits = date_limits) + 
											#scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme() + 
											labs(x = NULL, y = NULL, color = NULL) + 
											annotate("rect", xmin = gos_dates[1], xmax = gos_dates[2], ymin = 0, ymax = largest_val, alpha = 0.7, fill = "#bfbfbf") + 
											annotate("text", x = gos_dates[1]+10, y = largest_val+10, color = "#474747", size = 2.5, family = "Arial", label = "Data subject to change.") + 
											geom_hline(aes(yintercept=trend03, text=paste0("Most recent 3 months: ", prettyNum(trend03, big.mark=",", digits=1), " particles/person", sep="")), color=TRENDL_03_COLOR, linetype="dotted", alpha=0.8, linewidth=0.5) + 
											geom_hline(aes(yintercept=trend12, text=paste0("Most recent year:", prettyNum(trend12, big.mark=",", digits=1), " particles/person", sep="")), color=TRENDL_12_COLOR, linetype="dotted", alpha=0.6, linewidth=0.5) + 
											geom_point(aes(x = date_to_plot, y = val, text=paste0("Week of ", printy_dates(date_to_plot-7), " - ", printy_dates(date_to_plot), "\n", prettyNum(val, big.mark=",", digits=1), " particles/person (on average).", sep="")), shape = 1, size = 2, alpha=0.9) + 
											geom_line(aes(x = date_to_plot, y = val), alpha=0.4, na.rm = TRUE)
		
		#gplot$x$data[[1]]$hoverinfo <- "none"	# Supposed to get rid of the popup on the rect annotation but doesn't work

		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	# Accepts a dataframe and an integer representing # of months to plot.
	#
	# Generate a ggplotly object of the variant data within the date window.
	#
	plotVariants <- function(df_variants, months_to_plot) {
		#print("##### plotVariants called!")

		dlab <- case_when(
			months_to_plot == 1 ~ DATE_LABELS[1],
			months_to_plot == 3 ~ DATE_LABELS[2],
			months_to_plot == 6 ~ DATE_LABELS[3],
			months_to_plot == 12 ~ DATE_LABELS[4],
			months_to_plot == 24 ~ DATE_LABELS[5],
		)

		dbrk <- case_when(
			months_to_plot == 1 ~ DATE_BREAKS[1],
			months_to_plot == 3 ~ DATE_BREAKS[2],
			months_to_plot == 6 ~ DATE_BREAKS[3],
			months_to_plot == 12 ~ DATE_BREAKS[4],
			months_to_plot == 24 ~ DATE_BREAKS[5],
		)
				
		#end_date <- max(df_variants$date_to_plot, na.rm=TRUE)
    end_date <- this_week
		date_limits <- c(end_date %m-% months(months_to_plot), end_date)
		df_plot <- df_variants %>% filter(date_to_plot >= date_limits[1] & date_to_plot <= date_limits[2])
    
		gplot <- ggplot(df_plot, aes(fill=color_group, y=total_pct, x=date_to_plot)) + labs(y = "", x = "") + 
#							geom_bar(position="stack", stat="identity", aes(fill=factor(color_group), text=paste0("Week that starts ", printy_dates(date_to_plot), "\nVariant family: ", color_group, "\nProportion: ", prettyNum(total_pct, digits=2), "%", sep=""))) + 
							geom_bar(position="stack", stat="identity", aes(fill=factor(color_group), text=paste0("Week of ", printy_dates(date_to_plot-7), " - ", printy_dates(date_to_plot), "\nVariant family: ", color_group, "\nProportion: ", prettyNum(total_pct, digits=2), "%", sep=""))) + 
							scale_fill_brewer(type="div", palette = "RdYlBu", direction = -1, na.value = "#a8a8a8") + 
							labs(x="", y="", fill=NULL) + 
		          scale_x_date(date_breaks = dbrk, date_labels = dlab, limits = date_limits) + 
#							scale_y_continuous(name = NULL, limits = c(0, 110), breaks = c(0, 25, 50, 75, 100)) +  DOESN'T WORK FOR SOME REASON!?
							plot_theme() + 
							theme(legend.position = "right", legend.title=element_blank())

		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}


  updateAllPlots <- function() {
	  print("##### updateAllPlots called!")

# Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in dplyr 1.1.0.
# Please use `reframe()` instead.
# When switching from `summarise()` to `reframe()`, remember that `reframe()` always returns an ungrouped data
#  frame and adjust accordingly.

		loc_id <- controlRV$mapClick[controlRV$mapIndex]

		if (loc_id == "WV") {
			plot_location <- "West Virginia"
		} else {
			if ((df_regions %>% filter(region_name == loc_id))$region_geolevel == "county") {
				#this_county <- (df_active_loc %>% filter(location_id == loc_id))$location_counties_served
				this_state <- "WV"
				plot_location <- paste0(loc_id, " county, ", this_state, sep="")
			} else {
				plot_location <- unique((df_active_loc %>% filter(location_id == loc_id))$location_common_name)
			}
		}
		
		lapply(1:length(controlRV$targetVec), function(i) {
		#for (target_index in controlRV$targetVec) {
			target_index <- controlRV$targetVec[i]
			
			plotSuffix <- tolower(DISEASES[target_index])
			
	    df_abundance <- getAbundanceData(loc_id, target_index)

			if (target_index == 1) {
		    df_variants <- getVariantData(loc_id)
				output$plot_title_covid <- renderText(paste0(toupper(plotSuffix), " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_covid <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", toupper(plotSuffix), "!"))
				}

				output$plotsq_title_covid <- renderText(paste0(toupper(plotSuffix), " variants in ", plot_location, " wastewater"))
				if (nrow(df_variants) > 0) {
					output$plotsq_covid <- renderPlotly({
						plotVariants(df_variants, controlRV$viewMonths[target_index]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plotsq_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No variant data for ", toupper(plotSuffix), "!"))
				}

			} else if (target_index == 2) {
				output$plot_title_flua <- renderText(paste0(toupper(plotSuffix), " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_flua <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_flua <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", toupper(plotSuffix), "!"))
				}
			} else if (target_index == 3) {
				output$plot_title_flub <- renderText(paste0(toupper(plotSuffix), " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_flub <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_flub <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", toupper(plotSuffix), "!"))
				}
			} else if (target_index == 4) {
				output$plot_title_rsv <- renderText(paste0(toupper(plotSuffix), " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_rsv <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_rsv <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", TARGETS[target_index], "!"))
				}
			}
			
		}) # end of lapply
  }

  
	#
	# Update the alert status elements (usually on reaction to map click or site 
	# selection).
	#
	updateStatus <- function() {
	  print("##### updateStatus called!")
	  # Update abundance and trend text blocks, and the dominant variant text block (if applicable).

#	  print("0 updateStatus called!")
		
		# Get the current map click
		loc_id <- controlRV$mapClick[controlRV$mapIndex]
		
#	  print("1 updateStatus called!")
		
		lapply(1:length(controlRV$targetVec), function(i) {
			target_index <- controlRV$targetVec[i]
			
			plotSuffix <- tolower(DISEASES[target_index])
#			print(paste0("target index is ", target_index, " and plot suffix is ", plotSuffix, sep=""))
			
			# Generate the value for the abundance block.
			vec_abundance <- getAbundance(loc_id, target_index)
	
#	  	print("2 updateStatus called!")
	
			# Generate the value for the trend block.
			vec_trend <- getTrend(loc_id, target_index)
	
#		  print("3 updateStatus called!")
	
			frunner <- paste0(
				"document.getElementById('abundance_text_", plotSuffix, "').style.backgroundColor = '", vec_abundance[2], "';",
				"document.getElementById('trend_text_", plotSuffix, "').style.backgroundColor = '", vec_trend[2], "';", sep="")
			runjs(frunner)
			
			if (target_index == 1) {
				# Generate the value for the variant block.
				vec_variant <- getDominantVariant(loc_id, target_index)
				frunner2 <- paste0(
					"document.getElementById('variant_text_", plotSuffix, "').style.backgroundColor = '", vec_variant[2], "';", sep="")
				runjs(frunner2)
				
				# Write the text blocks. Be nice to construct these variable names from visPlot but can't figure that out.
				output$abundance_covid <- renderText(vec_abundance[1])
				output$trend_covid <- renderText(vec_trend[1])
				output$variant_covid <- renderText(vec_variant[1])
			
			} else if (target_index == 2) {
				output$abundance_flua <- renderText(vec_abundance[1])
				output$trend_flua <- renderText(vec_trend[1])
			
			} else if (target_index == 3) {
				output$abundance_flub <- renderText(vec_abundance[1])
				output$trend_flub <- renderText(vec_trend[1])
			
			} else if (target_index == 4) {
				output$abundance_rsv <- renderText(vec_abundance[1])
				output$trend_rsv <- renderText(vec_trend[1])
			
			}
		}) # end lapply
		
#	  print("4 updateStatus called!")
	}
	
	
	#
	# Update content in the selection info block (usually on reaction to map click or 
	# site selection).
	#
	updateSelectionInfo <- function() {
	  print("##### updateSelectionInfo called!")
	  
		loc_id <- controlRV$mapClick[controlRV$mapIndex]
		
		# Rollup location ids into a vector, even if it is a single facility.
		# Faciliates calculation of completenes and freshness.
		# Default is a vector containing the facility location id.
		loc_ids <- c(loc_id)
		
		if ("WV" %in% loc_ids) {
			# State specific elements
			loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			county_ids <- unique((df_active_loc %>% filter(location_id %in% loc_ids))$location_counties_served)

			title_text <- "The State of West Virginia"
			selection_text <- "The State of West Virginia"

		} else if ((df_regions %>% filter(region_name == loc_id))$region_geolevel == "county") {
			# County specific elements
			loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
			county_ids <- unique((df_active_loc %>% filter(location_id %in% loc_ids))$location_counties_served)
			
			county_name <- loc_id

			title_text <- paste0(loc_id, " county, WV", sep="")
			selection_text <- paste0(loc_id, " county", sep="")

		} else {
			# Facility specific elements
			loc_ids <- unique((df_active_loc %>% filter(location_id == loc_id))$location_id)
			county_ids <- unique((df_active_loc %>% filter(location_id %in% loc_ids))$location_counties_served)
			county_name <- county_ids[[1]]
			
			loc_name <- unique((df_active_loc %>% filter(location_id == loc_id))$location_common_name)
			
			title_text <- loc_name
			selection_text <- paste0("The ", loc_name, " facility ", sep="")

		}

		popserved_total <- sum(distinct(df_active_loc %>% filter(location_id %in% loc_ids), location_id, location_population_served)$location_population_served)
		#popserved_pct <- 100 * popserved_total / sum((resources$county %>% filter(county_id %in% county_ids))$county_population)
		
		if (popserved_total == -1) {
			popserved_total <- "Unknown"
		}

		num_counties <- length(county_ids)
		counties_text <- "counties"
		if (num_counties == 1) {
			counties_text <- paste0(county_name, " county, WV.", sep="")
		} else {
			counties_text <- paste0(num_counties, " counties.", sep="")
		}

		num_facilities <- length(loc_ids)
		facility_text <- "facilities"
		if (num_facilities == 1) {
			facility_text <- "facility"
		}
		
		selection_details_text <- paste0(
			selection_text, " supports ", num_facilities, " active WaTCH ", facility_text, " serving ", 
			prettyNum(popserved_total, big.mark=","), " residents in ", counties_text, sep=""
		)

		lapply(1:length(controlRV$targetVec), function(i) {
		#for (target_index in controlRV$targetVec) {
			target_index <- controlRV$targetVec[i]
		
			plotSuffix <- tolower(DISEASES[target_index])
			
			# Calculate the freshness of the abundance data.
			#	
			df_fresh <- dflist_rs[[target_index]] %>% filter(location_id %in% loc_ids)
			
			most_recent_sample_date <- max(df_fresh$date_to_plot, na.rm=TRUE)
			most_recent_sample_win <- c(most_recent_sample_date-7, most_recent_sample_date)
	
			sample_freshness_text <- paste0(
				printy_dates(most_recent_sample_win[1]), 
				" - ", printy_dates(most_recent_sample_win[2]), 
				sep="")
	
			# Calculate the completeness of the abundance data.
			#		
			most_recent_contributors <- df_fresh %>% filter(date_to_plot == most_recent_sample_date)
			most_recent_contributor_count <- length(unique(most_recent_contributors$location_id))
			if (most_recent_contributor_count > 1) {
				site_suffix <- "sites"
			} else {
				site_suffix <- "site"
			}
			sample_completeness_text <- paste0(
				"and includes ", 
				most_recent_contributor_count, " reporting ", site_suffix, " (", 
				prettyNum(100*(most_recent_contributor_count/num_facilities), digits=1), "%).", 
				sep="")
	
			# Calculate the freshness and completeness of the variant data (COVID only).
			#
			if (target_index == 1) {
				df_fresh_sq <- df_seqr %>% filter(location_id %in% loc_ids)

				if (nrow(df_fresh_sq) > 0) {
					most_recent_sample_date_sq <- max(df_fresh_sq$date_to_plot, na.rm=TRUE)
					most_recent_sample_win_sq <- c(most_recent_sample_date_sq-7, most_recent_sample_date_sq)
					sample_freshness_text_sq <- paste0(
						printy_dates(most_recent_sample_win_sq[1]), " - ", 
						printy_dates(most_recent_sample_win_sq[2]), sep="")

					most_recent_contributors_sq <- df_fresh_sq %>% filter(date_to_plot == most_recent_sample_date_sq)
					most_recent_contributor_count_sq <- length(unique(most_recent_contributors_sq$location_id))
					if (most_recent_contributor_count_sq > 1) {
						site_suffix_sq <- "sites"
					} else {
						site_suffix_sq <- "site"
					}
					sample_completeness_text_sq <- paste0(
						"and includes ", 
						most_recent_contributor_count_sq, " reporting ", site_suffix_sq, ".", 
						sep="")

				} else {
					sample_freshness_text_sq <- paste0("not reported for this region.", sep="")
					sample_completeness_text_sq <- ""
				}
			}
	
			# Print the title, selection, and last_sample strings to the UI.
			#
			if (target_index == 1) {
				output$selection_title_covid <- renderText(title_text)
				output$selection_details_covid <- renderText(selection_details_text)
				output$selection_freshness_covid <- renderText(sample_freshness_text)
				output$selection_completeness_covid <- renderText(sample_completeness_text)
				output$selectionsq_freshness_covid <- renderText(sample_freshness_text_sq)
				output$selectionsq_completeness_covid <- renderText(sample_completeness_text_sq)
			} else if (target_index == 2) {
				output$selection_title_flu <- renderText(title_text)
				output$selection_details_flu <- renderText(selection_details_text)
				output$selection_freshness_flua <- renderText(sample_freshness_text)
				output$selection_completeness_flua <- renderText(sample_completeness_text)
			} else if (target_index == 3) {
				output$selection_title_flu <- renderText(title_text)
				output$selection_details_flu <- renderText(selection_details_text)
				output$selection_freshness_flub <- renderText(sample_freshness_text)
				output$selection_completeness_flub <- renderText(sample_completeness_text)
			} else if (target_index == 4) {
				output$selection_title_rsv <- renderText(title_text)
				output$selection_details_rsv <- renderText(selection_details_text)
				output$selection_freshness_rsv <- renderText(sample_freshness_text)
				output$selection_completeness_rsv <- renderText(sample_completeness_text)
			}
			
		}) # end lapply
		
	}

	
	getDownloadTableRS <- function(targetIndex) {
		return(dflist_rs[targetIndex])
	}
		

	output$download_data_covid <- downloadHandler(
		filename = "watch-wv_COVID.csv",
		content = function(file) {
			write.csv(getDownloadTableRS(1), file, row.names = FALSE)
		}
	)

	output$download_data_flua <- downloadHandler(
		filename = "watch-wv_FluA.csv",
		content = function(file) {
			write.csv(getDownloadTableRS(2), file, row.names = FALSE)
		}
	)

	output$download_data_flub <- downloadHandler(
		filename = "watch-wv_FluB.csv",
		content = function(file) {
			write.csv(getDownloadTableRS(3), file, row.names = FALSE)
		}
	)

	output$download_data_rsv <- downloadHandler(
		filename = "watch-wv_RSV.csv",
		content = function(file) {
			write.csv(getDownloadTableRS(4), file, row.names = FALSE)
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
	

	# Accepts a menu option.
	#
	# Change the geolevel to the given option in the given map.
	#
	changeGeolevel <- function(clicked) {

		mapProxy <- getMapProxy(controlRV$mapIndex)

		# Update some reactive elements
		controlRV$activeGeoLevel[controlRV$mapIndex] <- clicked
		
		if (clicked == "County") {
			mapProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircles(data = merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[controlRV$mapIndex]], by.x="location_id", by.y="region_name"),
										 layerId = ~location_id, 
										 lat = ~location_lat, 
										 lng = ~location_lng, 
										 radius = ~dotsize, 
#										 radius = 5, 
										 stroke = TRUE,
										 weight = 4, 
										 opacity = 0.5,
										 color = ~abundance_color,
										 fill = TRUE,
										 fillColor = ~trend_color,
										 group = "facility", 
										 label = ~as.character(paste0(location_common_name, " (" , abundance_level, " & ", trend, ")")), 
										 fillOpacity = 0.6) %>%
					addPolygons( 
						data = merge(county_spdf, dflist_alerts[[controlRV$mapIndex]], by.x="NAME", by.y="region_name"), 
						layerId = ~NAME, 
						fillColor = ~trend_color, 
						stroke=TRUE, 
						fillOpacity = 0.7, 
						color="#000000", 
						weight=0.5, 
						group="county",
						label = ~as.character(paste0(NAME, " (", abundance_level, " & ", trend,")")), 
						highlightOptions = highlightOptions(
							weight = 1,
							color = "#00F900",
							fillOpacity = 1.0,
							bringToFront = TRUE)
					)
		} else {
			if (clicked == "Facility") {
				mapProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = merge(county_spdf, dflist_alerts[[controlRV$mapIndex]], by.x="NAME", by.y="region_name"), 
							layerId = ~NAME, 
							#fillColor = ~mypalette(colorby), 
							stroke=TRUE,
							fillOpacity = 0, 
							color="#666666", 
							weight=1, 
							group="county"
						) %>% 
						addCircles(
							data = merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[controlRV$mapIndex]], by.x="location_id", by.y="region_name"),
							layerId = ~location_id, 
							lat = ~location_lat, 
							lng = ~location_lng, 
							radius = ~dotsize, 
							stroke = TRUE,
							weight = 2, 
							opacity = 0.9,
							color = ~abundance_color,
							fill = TRUE,
							fillColor = ~trend_color,
							group = "facility", 
							label = ~as.character(paste0(location_common_name, " (" , abundance_level, " & ", trend, ")")), , 
							highlightOptions = highlightOptions(
								weight = 3,
								color = "#00f900",
								fillOpacity = 0.9,
								bringToFront = TRUE),
							fillOpacity = 0.6)
			}
		}

		# Update the plots only if necessary
		if (controlRV$mapClick[controlRV$mapIndex] != "WV") {
			controlRV$mapClick[controlRV$mapIndex] <- "WV"
			updateAllPlots()
			updateSelectionInfo()
			updateStatus()
		}
		

	}


	# Accepts a menu option.
	#
	# Change the geolevel to the given option in the given map.
	#
	changeMapColor <- function(clicked) {

		mapProxy <- getMapProxy(controlRV$mapIndex)

		# Update some reactive elements
		controlRV$activeMapColor[controlRV$mapIndex] <- clicked
		
		
		df_map_loc <- merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[controlRV$mapIndex]], by.x="location_id", by.y="region_name", all.x = TRUE) 
		df_map_loc <- df_map_loc %>% rename("Lab" = location_primary_lab, "Trend" = trend, "Abundance" = abundance)
				
		df_map_county <- merge(county_spdf, dflist_alerts[[controlRV$mapIndex]], by.x="NAME", by.y="region_name", all.x = TRUE)
		df_map_county <- df_map_county %>% rename("Trend" = trend, "Abundance" = abundance)
		
		if (clicked == "County") {
			mapProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircles(data = df_map_loc,
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
										 fillColor = ~watchPal(clicked),
										 group = "facility", 
										 label = ~as.character(paste0(location_common_name, " (" , prettyNum(location_population_served, big.mark=","), ")")), 
										 fillOpacity = 0.6) %>%
					addPolygons( 
						data = df_map_county, 
						layerId = ~NAME, 
						fillColor = ~watchPal(clicked), 
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
			if (clicked == "Facility") {
				mapProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = df_map_county, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(colorby), 
							stroke=TRUE,
							fillOpacity = 0, 
							color="#666666", 
							weight=1, 
							group="county"
						) %>% 
						addCircles(data = df_map_loc,
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
											fillColor = ~watchPal(clicked),
											group = "facility", 
											label = ~as.character(paste0(location_common_name, " (" , prettyNum(location_population_served, big.mark=","), ")")), , 
											highlightOptions = highlightOptions(
												weight = 3,
												color = "#fff",
												fillOpacity = 0.9,
												bringToFront = TRUE),
											fillOpacity = 0.6)
			}
		}
	}
	

	# Accepts a map marker.
	#
	# Respond to a click on the given map marker in the given map.
	#
	clickMapMarker <- function(clicked) {
    if (length(clicked) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- clicked$id
    }
		
		loc_id <- unique((df_active_loc %>% filter(location_common_name == clickedLocation))$location_id)
    
		#print(loc_id)
		
		if (loc_id %in% df_rs$location_id | clickedLocation == "WV") {
		
			# Update the reactive element
			controlRV$mapClick[controlRV$mapIndex] <- clickedLocation
		
		  updateAllPlots()
			updateSelectionInfo()
			updateStatus()

		} else {
			print(paste0("No data for ", clickedLocation))
		}
	}


	# Accepts a map shape.
	#
	# Respond to a click on the given map shape in the given map.
	#
	clickMapShape <- function(clicked) {
#    print("##### clickMapShape called!")
		if (length(clicked) == 0) {
			clickedLocation <- "WV"
		} else {
    	clickedLocation <- clicked$id
			controlRV$mapClickLat[controlRV$mapIndex] <- clicked$lat
			controlRV$mapClickLng[controlRV$mapIndex] <- clicked$lng
    }

    
		if (clickedLocation %in% df_active_loc$location_id | clickedLocation %in% df_active_loc$location_counties_served | clickedLocation == "WV") {
			controlRV$mapClick[controlRV$mapIndex] <- clickedLocation
		} else {
			clickedLocation <- "WV"
			controlRV$mapClick[controlRV$mapIndex] <- clickedLocation
		}
		
    updateAllPlots()
		updateSelectionInfo()
		updateStatus()

	}
	

	# Accepts a map click.
	#
	# Respond to a generic off-marker click in the given map.
	#
	clickMapOffMarker <- function(clicked) {
#    print("##### clickMapOffMarker called!")
		# only respond if this click is in a new position on the map
		if (clicked$lat != controlRV$mapClickLat[controlRV$mapIndex] | clicked$lng != controlRV$mapClickLng[controlRV$mapIndex]) {
			#print("New position detected!")
			controlRV$mapClickLat[controlRV$mapIndex] <- 0
			controlRV$mapClickLng[controlRV$mapIndex] <- 0

			# Update the reactive element
			controlRV$mapClick[controlRV$mapIndex] <- "WV"
		
			updateAllPlots()
			updateSelectionInfo()
			updateStatus()
		}
	}
	
	
	
	#
	# Render the maps.
	#
	output$map_covid <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[controlRV$mapIndex]], by.x="location_id", by.y="region_name"),
												 layerId = ~location_id, 
												 lat = ~location_lat, 
												 lng = ~location_lng, 
												 radius = ~dotsize, 
												 stroke = TRUE,
												 weight = 4, 
												 opacity = 0.5,
												 color = ~abundance_color,
												 fill = TRUE,
												 fillColor = ~trend_color,
												 group = "facility", 
												 label = ~as.character(paste0(location_common_name, " (" , abundance_level, " & ", trend, ")")), 
												 fillOpacity = 0.6) %>%
			addPolygons( 
				data = merge(county_spdf, dflist_alerts[[controlRV$mapIndex]], by.x="NAME", by.y="region_name"), 
				layerId = ~NAME, 
				fillColor = ~trend_color, 
				stroke = TRUE, 
				fillOpacity = 0.7, 
				color = "#000000", 
				weight = 0.5, 
				group="county",
				label = ~as.character(paste0(NAME, " (", abundance_level, " & ", trend,")")), 
				highlightOptions = highlightOptions(
					weight = 0.5,
					color = "#00F900",
					#dashArray = "",
					fillOpacity = 1.0,
					bringToFront = TRUE)
			)		
	})


	output$map_flu <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[controlRV$mapIndex]], by.x="location_id", by.y="region_name"),
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
												 fillColor = ~trend_color,
												 group = "facility", 
												 label = ~as.character(paste0(location_common_name, " (" , abundance_level, " & ", trend, ")")), 
												 fillOpacity = 0.6) %>%
			addPolygons( 
				data = merge(county_spdf, dflist_alerts[[controlRV$mapIndex]], by.x="NAME", by.y="region_name"), 
				layerId = ~NAME, 
				fillColor = ~trend_color, 
				stroke=TRUE, 
				fillOpacity = 0.7, 
				color="#000000", 
				weight=0.5, 
				group="county",
				label = ~as.character(paste0(NAME, " (", abundance_level, " & ", trend,")")), 
				highlightOptions = highlightOptions(
					weight = 1,
					color = "#00F900",
					#dashArray = "",
					fillOpacity = 1.0,
					bringToFront = TRUE)
			)		
	})


	output$map_rsv <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = merge(df_active_loc %>% filter(location_category == "wwtp"), dflist_alerts[[controlRV$mapIndex]], by.x="location_id", by.y="region_name"),
												 layerId = ~location_id, 
												 lat = ~location_lat, 
												 lng = ~location_lng, 
												 radius = ~dotsize, 
												 stroke = TRUE,
												 weight = 4, 
												 opacity = 0.5,
												 color = ~abundance_color,
												 fill = TRUE,
												 fillColor = ~trend_color,
												 group = "facility", 
												 label = ~as.character(paste0(location_common_name, " (" , abundance_level, " & ", trend, ")")), 
												 fillOpacity = 0.6) %>%
			addPolygons( 
				data = merge(county_spdf, dflist_alerts[[controlRV$mapIndex]], by.x="NAME", by.y="region_name"), 
				layerId = ~NAME, 
				fillColor = ~trend_color, 
				stroke = TRUE, 
				fillOpacity = 0.7, 
				color = "#000000", 
				weight = 0.5, 
				group="county",
				label = ~as.character(paste0(NAME, " (", abundance_level, " & ", trend,")")), 
				highlightOptions = highlightOptions(
					weight = 0.5,
					color = "#00F900",
					#dashArray = "",
					fillOpacity = 1.0,
					bringToFront = TRUE)
			)		
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
    clickMapMarker(input$map_covid_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flu_marker_click, { 
    #print("##### Map MARKER click top")
    clickMapMarker(input$map_flu_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_rsv_marker_click, { 
    #print("##### Map MARKER click top")
    clickMapMarker(input$map_rsv_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)


	# 
	# React to map shape click
	#
  observeEvent(input$map_covid_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_covid_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  observeEvent(input$map_flu_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_flu_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  observeEvent(input$map_rsv_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_rsv_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)


	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	#
  observeEvent(input$map_covid_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_covid_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_flu_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_flu_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$map_rsv_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_rsv_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)



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
	# Change the map active geolayer.
	#
  observeEvent(input$geo_level, {
  	#print("geo level event fired")
		changeGeolevel(input$geo_level)
	}, ignoreInit = TRUE)

	#
	# Change the map color scheme.
	#
  observeEvent(input$map_color, {
  	#print("geo level event fired")
		#changeMapColor(input$map_color)
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



	#
	# Open the COVID color key popups.
	#
	onevent("click", "abundance_level_key_covid", togglePanels(on=c("abundance_level_key_covid_popup")))
	onevent("click", "trend_key_covid", togglePanels(on=c("trend_key_covid_popup")))

	#
	# Close the COVID color key popups.
	#
	observeEvent(input$abundance_level_key_covid_popup_close,{
    togglePanels(off=c("abundance_level_key_covid_popup"))
	}, ignoreInit = TRUE)

	observeEvent(input$trend_key_covid_popup_close,{
    togglePanels(off=c("trend_key_covid_popup"))
	}, ignoreInit = TRUE)
	

	#
	# Open the FLU color key popups.
	#
	onevent("click", "abundance_level_key_flu", togglePanels(on=c("abundance_level_key_flu_popup")))
	onevent("click", "trend_key_flu", togglePanels(on=c("trend_key_flu_popup")))

	#
	# Close the FLU color key popups.
	#
	observeEvent(input$abundance_level_key_flu_popup_close,{
    togglePanels(off=c("abundance_level_key_flu_popup"))
	}, ignoreInit = TRUE)

	observeEvent(input$trend_key_flu_popup_close,{
    togglePanels(off=c("trend_key_flu_popup"))
	}, ignoreInit = TRUE)
	

	#
	# Open the RSV color key popups.
	#
	onevent("click", "abundance_level_key_rsv", togglePanels(on=c("abundance_level_key_rsv_popup")))
	onevent("click", "trend_key_rsv", togglePanels(on=c("trend_key_rsv_popup")))

	#
	# Close the RSV color key popups.
	#
	observeEvent(input$abundance_level_key_rsv_popup_close,{
    togglePanels(off=c("abundance_level_key_rsv_popup"))
	}, ignoreInit = TRUE)

	observeEvent(input$trend_key_rsv_popup_close,{
    togglePanels(off=c("trend_key_rsv_popup"))
	}, ignoreInit = TRUE)
	

	# 
	# React to plot click.
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_rs_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)

})







#