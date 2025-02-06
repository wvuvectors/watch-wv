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
		activeColorBy = c("lab", "lab", "lab"), 
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
      print("tab = covid")
      controlRV$mapIndex <- 1
      controlRV$targetVec <- c(1)

			updateAllPlots()
			updateSelectionInfo()
			updateStatus()

    } else if(input$nav == "Influenza"){
      print("tab = flu")
      controlRV$mapIndex <- 2
      controlRV$targetVec <- c(2,3)

			updateAllPlots()
			updateSelectionInfo()
			updateStatus()
    }
  })
	
	# The following three observers keep the map and controls synched across tabs.
	
	# 
	# Update the map colorBy when the map content changes
	#
	observe({
		updateSelectInput(
			session, 
			inputId = "map_color", 
			selected = controlRV$activeColorBy[controlRV$mapIndex]
		)
	}) 

	# 
	# Update the map geoLevel when the map content changes
	#
	observe({
		updateSelectInput(
			session, 
			inputId = "geo_level", 
			selected = controlRV$activeGeoLevel[controlRV$mapIndex]
		)
	}) 

	# 
	# Update the plot view range when the map content changes
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
	# 	i <- 2
	# 	txt <- ALERT_LEVEL_STRINGS[i]
	# 	color <- ALERT_LEVEL_COLORS[i]
		raw_val <- (dflist_alerts[[target_index]] %>% filter(region_name == loc_id))$latest_abundance
		color <- ALERT_LEVEL_COLORS[2]
		
		val <- prettyNum(raw_val, big.mark=",", digits=1)
		
		return(c(val, color))
	}
	
	
	getTrend <- function(loc_id, target_index) {
	# 	i <- 3
	# 	txt <- ALERT_LEVEL_STRINGS[i]
	# 	color <- ALERT_LEVEL_COLORS[i]
		print((dflist_alerts[[target_index]] %>% filter(region_name == loc_id))$trend)
		
		txt <- (dflist_alerts[[target_index]] %>% filter(region_name == loc_id))$trend
		color <- ALERT_LEVEL_COLORS[2]
	
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
	# Returns a dataframe of abundance data for the current bio target at the given location.
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
	# Returns a dataframe of variant data for the current bio target at the given location.
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
	# Generates a ggplotly object of the abundance data within the date window.
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
		
		# Calculate the mean signal for each target over the last month and last year.
		trend01 <- calcTrend(df_abundance, 1)
		trend06 <- calcTrend(df_abundance, 6)
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
											geom_hline(aes(yintercept=trend01, text=paste0("Most recent month: ", prettyNum(trend01, big.mark=",", digits=1), " virus/person", sep="")), color=TRENDL_01_COLOR, linetype="dotted", alpha=1.0, linewidth=0.5) + 
											geom_hline(aes(yintercept=trend06, text=paste0("Most recent 6 months:", prettyNum(trend06, big.mark=",", digits=1), " virus/person", sep="")), color=TRENDL_06_COLOR, linetype="dotted", alpha=0.8, linewidth=0.5) + 
											geom_hline(aes(yintercept=trend12, text=paste0("Most recent year:", prettyNum(trend12, big.mark=",", digits=1), " virus/person", sep="")), color=TRENDL_12_COLOR, linetype="dotted", alpha=0.6, linewidth=0.5) + 
											geom_point(aes(x = date_to_plot, y = val, text=paste0("Week of ", printy_dates(date_to_plot-7), " - ", printy_dates(date_to_plot), "\n", prettyNum(val, big.mark=",", digits=1), " virus/person (on average).", sep="")), shape = 1, size = 2, alpha=0.9) + 
#											geom_point(aes(x = date_to_plot, y = val, text=paste0("Week that ends ", printy_dates(date_to_plot), "\n", prettyNum(val, big.mark=",", digits=1), " virus/person (on average).", sep="")), shape = 1, size = 2, alpha=0.9) + 
											geom_line(aes(x = date_to_plot, y = val), alpha=0.4, na.rm = TRUE)
		
		#gplot$x$data[[1]]$hoverinfo <- "none"	# Supposed to get rid of the popup on the rect annotation but doesn't work

		ggplotly(gplot, tooltip="text") %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}
	

	# Accepts a dataframe and an integer representing # of months to plot.
	#
	# Generates a ggplotly object of the variant data within the date window.
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
		
		for (target_index in controlRV$targetVec) {
			
			plotSuffix <- tolower(DISEASES[target_index])
			
	    df_abundance <- getAbundanceData(loc_id, target_index)

			if (target_index == 1) {
		    df_variants <- getVariantData(loc_id)
				output$plot_title_covid <- renderText(paste0(plotSuffix, " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_covid <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", TARGETS[target_index], "!"))
				}

				output$plotsq_title_covid <- renderText(paste0(plotSuffix, " variants in ", plot_location, " wastewater"))
				if (nrow(df_variants) > 0) {
					output$plotsq_covid <- renderPlotly({
						plotVariants(df_variants, controlRV$viewMonths[target_index]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plotsq_covid <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No variant data for ", TARGETS[target_index], "!"))
				}

			} else if (target_index == 2) {
				output$plot_title_flua <- renderText(paste0(TARGETS[target_index], " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_flua <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_flua <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", TARGETS[target_index], "!"))
				}
			} else if (target_index == 3) {
				output$plot_title_flub <- renderText(paste0(TARGETS[target_index], " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_flub <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_flub <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", TARGETS[target_index], "!"))
				}
			} else if (target_index == 4) {
				output$plot_title_rsv <- renderText(paste0(TARGETS[target_index], " abundance in ", plot_location, " wastewater"))
				if (nrow(df_abundance) > 0) {
					output$plot_rsv <- renderPlotly({
						plotAbundance(df_abundance, controlRV$viewMonths[controlRV$mapIndex]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
					})
				} else {
					output$plot_rsv <- renderPlotly({ggplotly(ggplot()) %>% config(displayModeBar = FALSE)})
					print(paste0("No data for ", TARGETS[target_index], "!"))
				}
			}
			
		} # end of for loop
  }

  
	#
	# Updates the alert status elements (usually on reaction to map click or site 
	# selection).
	#
	updateStatus <- function() {
	  print("##### updateStatus called!")
	  # Update abundance and trend text blocks, and the dominant variant text block (if applicable).

#	  print("0 updateStatus called!")
		
		# Get the current map click
		loc_id <- controlRV$mapClick[controlRV$mapIndex]
		
#	  print("1 updateStatus called!")
		
		for (target_index in controlRV$targetVec) {
			
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
		}
#	  print("4 updateStatus called!")
	}
	
	
	#
	# Updates content in the selection info block (usually on reaction to map click or 
	# site selection).
	#
	updateSelectionInfo <- function() {
	  print("##### updateSelectionInfo called!")
	  
		loc_id = controlRV$mapClick[controlRV$mapIndex]
		
		if (loc_id == "WV") {
			
			loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			df_this <- df_rs
			df_this_sq <- df_seqr
			
			total_popserved <- sum(distinct(df_active_loc %>% filter(location_category == "wwtp"), location_id, location_population_served)$location_population_served)
			if (total_popserved == -1) {
				total_popserved = "Unknown"
			}

			num_facilities <- n_distinct(df_this$location_id)
			num_counties <- n_distinct((df_active_loc %>% filter(location_category == "wwtp"))$location_counties_served)
			
			selection_text <- "West Virginia"
			title_text <- "State of West Virginia"
			selection_details_text <- paste0(
					"The State of West Virginia represents ", num_facilities, " active treatment facilities serving ", 
					prettyNum(total_popserved, big.mark=","), " residents across ", 
					num_counties, " counties.", sep="")
		} else {
			
			if ((df_regions %>% filter(region_name == loc_id))$region_geolevel == "county") {
				
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_id & location_category == "wwtp"))$location_id
				num_facilities <- n_distinct(loc_ids)
				facility_post <- "facilities"
				if (num_facilities == 1) {
					facility_post = "facility"
				}
				
				df_this <- df_rs %>% filter(location_id %in% loc_ids)
				df_this_sq <- df_seqr %>% filter(location_id %in% loc_ids)
				
				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id %in% loc_ids), location_id, location_population_served)$location_population_served)
				pct_served <- 100 * total_popserved / (resources$county %>% filter(county_id == loc_id))$county_population

				selection_text <- paste0(loc_id, " county")
				title_text <- paste0(loc_id, " county, WV", sep="")
				selection_details_text <- paste0(loc_id, " county in WV supports ", 
					 num_facilities, " WaTCH ", facility_post, ", serving a total of ", 
					 prettyNum(total_popserved, big.mark=","), " residents (", 
					 prettyNum(pct_served, digits=1), "% of the county).",
					 sep="")
			} else {
				loc_ids <- unique((df_active_loc %>% filter(location_id == loc_id))$location_id)
				this_wwtp_id <- unique((df_active_loc %>% filter(location_id == loc_id))$location_primary_wwtp_id)
				loc_name <- unique((df_active_loc %>% filter(location_id == loc_id))$location_common_name)
				
				df_this <- df_rs %>% filter(location_id == loc_ids)
				df_this_sq <- df_seqr %>% filter(location_id == loc_ids)

				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id == loc_ids), location_id, location_population_served)$location_population_served)
				if (total_popserved == -1) {
					total_popserved = "Unknown"
				}

				num_facilities <- 1
				this_county <- (df_active_loc %>% filter(location_id == loc_ids))$location_counties_served
			
				selection_text <- loc_name
				title_text <- loc_name
				selection_details_text <- paste0("The ", loc_name, " facility serves ", 
					 prettyNum(total_popserved, big.mark=","), " residents in ", 
					 this_county, " county, WV.",
					 sep="")
			}
			
		}
		
		for (target_index in controlRV$targetVec) {
		
			plotSuffix <- tolower(DISEASES[target_index])
			
			# Calculate the freshness of the abundance data.
			#	
			df_fresh <- df_this %>% filter(target == TARGETS[target_index])
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
	
			# Calculate the freshness of the variant data.
			if (target_index == 1) {
				if (nrow(df_this_sq) > 0) {
					most_recent_sample_date_sq <- max(df_this_sq$date_to_plot, na.rm=TRUE)
					most_recent_sample_win_sq <- c(most_recent_sample_date_sq-7, most_recent_sample_date_sq)
					sample_freshness_text_sq <- paste0(
						printy_dates(most_recent_sample_win_sq[1]), " - ", 
						printy_dates(most_recent_sample_win_sq[2]), ".", sep="")
				} else {
					sample_freshness_text_sq <- paste0("not reported for this region.", sep="")
				}
			}
	

			# Print the title, selection, and last_sample strings to the UI
			if (target_index == 1) {
				output$selection_covid <- renderText(selection_text)
				output$selection_details_covid <- renderText(selection_details_text)
				output$selection_title_covid <- renderText(title_text)
				output$selection_freshness_covid <- renderText(sample_freshness_text)
				output$selection_completeness_covid <- renderText(sample_completeness_text)
				output$selectionsq_freshness_covid <- renderText(sample_freshness_text_sq)
			} else {
				if (target_index == 2) {
					output$selection_flu <- renderText(selection_text)
					output$selection_details_flu <- renderText(selection_details_text)
					output$selection_title_flu <- renderText(title_text)
					output$selection_completeness_flub <- renderText(sample_completeness_text)
					output$selection_completeness_flub <- renderText(sample_completeness_text)
				} else {
					if (target_index == 3) {
						output$selection_flu <- renderText(selection_text)
						output$selection_details_flu <- renderText(selection_details_text)
						output$selection_title_flu <- renderText(title_text)
						output$selection_completeness_flub <- renderText(sample_completeness_text)
						output$selection_completeness_flub <- renderText(sample_completeness_text)
					} else {
						if (target_index == 4) {
							output$selection_rsv <- renderText(selection_text)
							output$selection_title_rsv <- renderText(title_text)
							output$selection_freshness_rsv <- renderText(sample_freshness_text)
							output$selection_completeness_rsv <- renderText(sample_completeness_text)
						}
					}
				}
			}
		} # end for loop
		
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


	# Accepts a map (target) index.
	#
	#	Resets the zoom and position of the given map.
	#
	resetMap <- function() {
		#print("resetMap called!")
		mapProxy <- getMapProxy(controlRV$mapIndex)
		mapProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
	}
	

	# Accepts a menu option.
	#
	# Changes the geolevel to the given option in the given map.
	#
	changeGeolevel <- function(clicked) {

		mapProxy <- getMapProxy(controlRV$mapIndex)

		# Update some reactive elements
		controlRV$activeGeoLevel[controlRV$mapIndex] <- clicked
		
		if (clicked == "County") {
			mapProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
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
#										 fillColor = ~alertPal(current_fold_change_smoothed), 
										 fillColor = ~colorby,
										 group = "facility", 
										 label = ~as.character(paste0(location_common_name, " (" , prettyNum(location_population_served, big.mark=","), ")")), 
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
			if (clicked == "Facility") {
				mapProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = county_spdf, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(colorby), 
							stroke=TRUE,
							fillOpacity = 0, 
							color="#666666", 
							weight=1, 
							group="county"
						) %>% 
						addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
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
	#										fillColor = ~alertPal(current_fold_change_smoothed), 
											fillColor = ~colorby,
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

		# Update the plots only if necessary
		if (controlRV$mapClick[controlRV$mapIndex] != "WV") {
			controlRV$mapClick[controlRV$mapIndex] <- "WV"
			updateAllPlots()
			updateSelectionInfo()
			updateStatus()
		}
		

	}


	# Accepts a map marker.
	#
	# Responds to a click on the given map marker in the given map.
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
	# Responds to a click on the given map shape in the given map.
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
	# Responds to a generic off-marker click in the given map.
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
	# Renders the maps.
	#
	output$map_covid <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
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
												 fillColor = ~colorby,
												 group = "facility", 
												 label = ~as.character(paste0(location_common_name, " (" , prettyNum(location_population_served, big.mark=","), ")")), 
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


	output$map_flu <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircles(data = df_active_loc %>% filter(location_category == "wwtp"),
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
												 fillColor = ~colorby,
												 group = "facility", 
												 label = ~as.character(paste0(location_common_name, " (" , prettyNum(location_population_served, big.mark=","), ")")), 
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
	# OBSERVER FUNCTIONS - COVID TAB
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_covid_marker_click, { 
    #print("##### Map MARKER click top")
    clickMapMarker(input$map_covid_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_covid_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_covid_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	#
  observeEvent(input$map_covid_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_covid_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	#
	# Reset the map
	#
	observeEvent(input$map_reset, {
		#print("Reset covid map!")
		resetMap()
	}, ignoreInit = TRUE)

	#
	# Change the map active geolayer
	#
  observeEvent(input$geo_level, {
  	#print("geo level event fired")
		changeGeolevel(input$geo_level)
	}, ignoreInit = TRUE)

	#
	# Change the date range to view
	#
  observeEvent(input$view_range, {

		# Update some reactive elements
		controlRV$viewMonths[controlRV$mapIndex] <- as.numeric(input$view_range)
		
		# Update the plots
    updateAllPlots()

	}, ignoreInit = TRUE)

	#
	# Open the color key popup
	#
	onevent("click", "risk_level_key_covid", togglePanels(on=c("risk_level_key_popup")))

	#
	# Close the color key popup
	#
	observeEvent(input$risk_level_key_popup_close,{
    togglePanels(off=c("risk_level_key_popup"))
	}, ignoreInit = TRUE)
	

	# 
	# React to plot click
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_rs_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)


	###########################
	#
	# OBSERVER FUNCTIONS - FLU TAB
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_flu_marker_click, { 
#    print("##### Map MARKER click top")
    clickMapMarker(input$map_flu_marker_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_flu_shape_click, {
#    print("##### Map Shape Click!")
    clickMapShape(input$map_flu_shape_click)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

	# 
	# React to map click (off-marker). This is fired even if a marker or shape has been clicked.
	# In this case, the marker/shape listener fires first.
	#
  observeEvent(input$map_flu_click, { 
#		print("##### Map off-target click top")
    clickMapOffMarker(input$map_flu_click)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

# 	#
# 	# Reset the map
# 	#
# 	observeEvent(input$map_reset_flu, {
# 		#print("Reset flu map!")
# 		resetMap()
# 	}, ignoreInit = TRUE)
# 
# 	#
# 	# Change the map active geolayer
# 	#
#   observeEvent(input$geo_level_flu, {
# 		changeGeolevel(input$geo_level_flu)
# 	}, ignoreInit = TRUE)
# 
# 	#
# 	# Change the date range to view
# 	#
#   observeEvent(input$view_range_flu, {
# 
# 		# Update some reactive elements
# 		controlRV$viewMonths[2] <- as.numeric(input$view_range_flu)
# 		
# 		# Update the plots
#     updateAllPlots()
# 
# 	}, ignoreInit = TRUE)

	#
	# Open the color key popup
	#
	onevent("click", "risk_level_key_flu", togglePanels(on=c("risk_level_key_popup")))

	#
	# Close the color key popup
	#
	observeEvent(input$risk_level_key_popup_close,{
    togglePanels(off=c("risk_level_key_popup"))
	}, ignoreInit = TRUE)
	

	# 
	# React to plot click
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_rs_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)

})







#