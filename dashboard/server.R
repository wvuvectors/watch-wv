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
	rsLeafletProxy <- leafletProxy(mapId="map_rs", session)
	mgLeafletProxy <- leafletProxy(mapId="map_mg", session)
	stLeafletProxy <- leafletProxy(mapId="map_st", session)

	#
	# Initialize reactive values with some defaults.
	# Each element in a controlRV array corresponds to a UI tab.
	#
	controlRV <- reactiveValues(
								activeGeoLevel = c(GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT, GEOLEVELS_DEFAULT), 
								mapClick = c("WV", "WV", "WV"), 
								viewMonths = c(VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY, VIEW_RANGE_PRIMARY), 
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
	

	# Generate a seqr plot on the given data frame.
	#
	plotSeqr <- function(df_plot, date_win) {
		#print("plotSeqr called!")
    #View(df_plot)
    
		gplot <- ggplot(df_plot, aes(fill=lineage_group, y=val, x=date_to_plot)) + labs(y = "", x = "") + 
							geom_bar(position="stack", stat="identity", aes(fill=factor(color_group))) + 
							scale_fill_brewer(type="qual", palette = "Dark2") + 
							labs(x="", y="") + 
							scale_x_date(limits = c(today %m-% months(date_win), today)) + 
							plot_theme() +
							theme(legend.position = "right", legend.title=element_blank())

		ggplotly(gplot) %>% layout(clickmode = list("event"), xaxis = list(showspikes = TRUE, showline = TRUE, spikemode = "across", hovermode = "x"))
	}


	# Generate a dataframe for a basic plot.
	#
	getBasicDF <- function(inputTarget, inputLocus, loc_name, date_win) {
	  #print("getBasicDF called!")
	  
		df_targ <- df_rs %>% filter(target == inputTarget & target_genetic_locus == inputLocus)
		dates <- c(today %m-% months(date_win), today)
		
		if (loc_name == "WV") {
			df_plot <- df_targ %>% 
								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE))
		} else {
			#print(loc_name)
			if (controlRV$activeGeoLevel[1] == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_common_name == loc_name))$location_id
			}
		
			df_plot <- df_targ %>% 
								 filter(location_id %in% loc_ids & date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot) %>% 
								 arrange(date_to_plot) %>%
								 summarize(val := mean(target_copies_fn_per_cap, na.rm = TRUE))
		}

		return(df_plot)
	}
	
	
	getSeqrDF <- function(inputTarget, loc_name, date_win) {
	  #print("getBasicDF called!")
	  
		dates <- c(today %m-% months(date_win), today)
		
		if (loc_name == "WV") {
			df_plot <- df_seqr %>% 
								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot, lineage_group) %>% 
								 summarize(val := mean(percent, na.rm = TRUE),
								 					 color_group := color_group) %>% 
								 arrange(date_to_plot, lineage_group)
		} else {
		
			if (controlRV$activeGeoLevel[1] == "County") {
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
			} else {
				loc_ids <- (df_active_loc %>% filter(location_common_name == loc_name))$location_id
			}

			df_plot <- df_seqr %>% 
								 filter(location %in% loc_ids & date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
								 group_by(date_to_plot, lineage_group) %>% 
								 summarize(val := mean(percent, na.rm = TRUE),
								 					 color_group := color_group) %>% 
								 arrange(date_to_plot, lineage_group)
		}
		#View(df_plot)
		return(df_plot)
	}


  updatePlotsRS <- function() {
    df_plot1 <- getBasicDF(TARGETS_RS[1], GENLOCI_RS[1], controlRV$mapClick[1], controlRV$viewMonths[1])
    output$plot1_rs <- renderPlotly({plotBasic(df_plot1, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})
    
#    df_plot2 <- getBasicDF(TARGETS_RS[2], GENLOCI_RS[2], controlRV$mapClick[1], controlRV$viewMonths[1])
#    output$plot2_rs <- renderPlotly({plotBasic(df_plot2, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})
    
    df_plot3 <- getBasicDF(TARGETS_RS[3], GENLOCI_RS[3], controlRV$mapClick[1], controlRV$viewMonths[1])
    output$plot3_rs <- renderPlotly({plotBasic(df_plot3, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})
    
    df_plot4 <- getBasicDF(TARGETS_RS[4], GENLOCI_RS[4], controlRV$mapClick[1], controlRV$viewMonths[1])
    output$plot4_rs <- renderPlotly({plotBasic(df_plot4, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")})

    # Update the plot titles
    output$plot1_rs_title = renderText(paste0(TARGETS_RS[1], sep=""))
    output$plot2_rs_title = renderText(paste0(TARGETS_RS[2], sep=""))
    output$plot3_rs_title = renderText(paste0(TARGETS_RS[3], sep=""))
    output$plot4_rs_title = renderText(paste0(TARGETS_RS[4], sep=""))
  }

  
	#
	# Generate an alert color string based on the target levels at the given location(s).
	#
	getAlertColor <- function(inputTarget, inputLocus, locations) {
		
		df_targ <- df_rs %>% filter(
			target == inputTarget & 
			target_genetic_locus == inputLocus & 
			location_id %in% locations)
		
		last_date <- max(df_targ$date_primary, na.rm = TRUE)
		
		this_color <- case_when(
			last_date == -Inf ~ STALE_DATA_COLORS[4],
			ymd(today)-ymd(last_date) < STALE_DATA_THRESHOLDS[1] ~ ALERT_LEVEL_COLORS[1],
			ymd(today)-ymd(last_date) < STALE_DATA_THRESHOLDS[2] ~ ALERT_LEVEL_COLORS[2],
			ymd(today)-ymd(last_date) < STALE_DATA_THRESHOLDS[3] ~ ALERT_LEVEL_COLORS[3],
			ymd(today)-ymd(last_date) >= STALE_DATA_THRESHOLDS[3] ~ ALERT_LEVEL_COLORS[4]
		)
		
		return(this_color)
	}
	
	#
	# Generate a color string based on the data freshness at the given location(s).
	#
	getFreshnessColor <- function(inputTarget, inputLocus, locations) {
		
		df_targ <- df_rs %>% filter(
			target == inputTarget & 
			target_genetic_locus == inputLocus & 
			location_id %in% locations)
		
		last_date <- max(df_targ$date_primary, na.rm = TRUE)
		
		this_color <- case_when(
			last_date == -Inf ~ STALE_DATA_COLORS[4],
			ymd(today)-ymd(last_date) < STALE_DATA_THRESHOLDS[1] ~ STALE_DATA_COLORS[1],
			ymd(today)-ymd(last_date) < STALE_DATA_THRESHOLDS[2] ~ STALE_DATA_COLORS[2],
			ymd(today)-ymd(last_date) < STALE_DATA_THRESHOLDS[3] ~ STALE_DATA_COLORS[3],
			ymd(today)-ymd(last_date) >= STALE_DATA_THRESHOLDS[3] ~ STALE_DATA_COLORS[4]
		)
		
		return(this_color)
	}
	
	
	#
	# Write the site info block (reaction to map click or site selection).
	#
	updateSiteInfoRS <- function(geolevel, loc_name, date_win) {
		
		if (missing(geolevel)) { layer = controlRV$activeGeoLevel[1] }
		if (missing(loc_name)) { loc_name = controlRV$mapClick[1] }
		if (missing(date_win)) { dates = controlRV$viewMonths[1] }
		
		#print(loc_name)

		dates <- c(today %m-% months(VIEW_RANGE_PRIMARY), today)
		
		if (loc_name == "WV") {
			
			loc_ids <- unique((df_active_loc %>% filter(location_category == "wwtp"))$location_id)
			df_this <- df_rs
			
			title_text <- "State of West Virginia"

			total_popserved <- sum(distinct(df_active_loc %>% filter(location_category == "wwtp"), location_id, location_population_served)$location_population_served)
			if (total_popserved == -1) {
				total_popserved = "Unknown"
			}

			total_cap <- sum(distinct(df_active_wwtp, wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
			num_facilities <- n_distinct(df_this$location_id)
			num_counties <- n_distinct((df_active_loc %>% filter(location_category == "wwtp"))$location_counties_served)

			# Print the info
			#
			counties_post <- "counties"
			if (num_counties == 1) {
				counties_post <- "county"
			}
			facilities_post <- "facilities"
			if (num_facilities == 1) {
				facilities_post <- "facility"
			}
		
			# Print the info
			#
			output$site_rs_info <- renderText(
				paste0("Representing ", num_facilities, " active treatment facilities serving ", 
							 prettyNum(total_popserved, big.mark=","), " residents across ", 
							 num_counties, " counties.",
							 sep="")
			)

		} else {
			
			if (geolevel == "Facility") {
			
				loc_ids <- unique((df_active_loc %>% filter(location_common_name == loc_name))$location_id)
				this_wwtp_id <- unique((df_active_loc %>% filter(location_common_name == loc_name))$location_primary_wwtp_id)
			
				df_this <- df_rs %>% filter(location_id == loc_ids)

				title_text <- loc_name
				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id == loc_ids), location_id, location_population_served)$location_population_served)
				if (total_popserved == -1) {
					total_popserved = "Unknown"
				}

				total_cap <- sum(distinct(df_active_wwtp %>% filter(wwtp_id == this_wwtp_id), wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
				num_facilities <- 1
				this_county <- (df_active_loc %>% filter(location_id == loc_ids))$location_counties_served
			
				# Print the info
				#
				output$site_rs_info <- renderText(
					paste0("The ", loc_name, " facility serves ", 
								 prettyNum(total_popserved, big.mark=","), " residents in ", 
								 this_county, " county, WV.",
								 sep="")
				)
			} else {
				# county click!
				
				loc_ids <- (df_active_loc %>% filter(location_counties_served == loc_name & location_category == "wwtp"))$location_id
				num_facilities <- n_distinct(loc_ids)
				facility_post <- "facilities"
				if (num_facilities == 1) {
					facility_post = "facility"
				}
				
				df_this <- df_rs %>% filter(location_id %in% loc_ids)
				
				title_text <- paste0(loc_name, " county, WV", sep="")

				total_cap <- sum(distinct(df_active_wwtp %>% filter(wwtp_id %in% df_this), wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1
				total_popserved <- sum(distinct(df_active_loc %>% filter(location_id %in% loc_ids), location_id, location_population_served)$location_population_served)
				pct_served <- 100 * total_popserved / (resources$county %>% filter(county_id == loc_name))$county_population

				# Print the info
				#
				output$site_rs_info <- renderText(
					paste0(loc_name, " county in WV supports ", 
								 num_facilities, " WaTCH ", facility_post, ", serving a total of ", 
								 prettyNum(total_popserved, big.mark=","), " residents (", 
								 prettyNum(pct_served, digits=1), "% of the county).",
								 sep="")
				)
			}
			
		}
		
		# Print the title
		output$site_rs_title <- renderText(title_text)

		# Update the alert div color for each target.
		#
		bgchange1 <- paste0("document.getElementById('rs_CSL_1').style.backgroundColor = '", getAlertColor(TARGETS_RS[1], GENLOCI_RS[1], loc_ids), "';", sep="")
		runjs(bgchange1)
		bgchange2 <- paste0("document.getElementById('rs_CSL_2').style.backgroundColor = '", getAlertColor(TARGETS_RS[2], GENLOCI_RS[2], loc_ids), "';", sep="")
		runjs(bgchange2)
		bgchange3 <- paste0("document.getElementById('rs_CSL_3').style.backgroundColor = '", getAlertColor(TARGETS_RS[3], GENLOCI_RS[3], loc_ids), "';", sep="")
		runjs(bgchange3)
		bgchange4 <- paste0("document.getElementById('rs_CSL_4').style.backgroundColor = '", getAlertColor(TARGETS_RS[4], GENLOCI_RS[4], loc_ids), "';", sep="")
		runjs(bgchange4)


		# Calculate the freshness of the data for each target.
		#
		fresh_1 <- max((df_this %>% filter(target == TARGETS_RS[1] & target_genetic_locus == GENLOCI_RS[1]))$date_primary, na.rm = TRUE)
		fresh_2 <- max((df_this %>% filter(target == TARGETS_RS[2] & target_genetic_locus == GENLOCI_RS[2]))$date_primary, na.rm = TRUE)
		fresh_3 <- max((df_this %>% filter(target == TARGETS_RS[3] & target_genetic_locus == GENLOCI_RS[3]))$date_primary, na.rm = TRUE)
		fresh_4 <- max((df_this %>% filter(target == TARGETS_RS[4] & target_genetic_locus == GENLOCI_RS[4]))$date_primary, na.rm = TRUE)

		if (fresh_1 == -Inf) {
			fresh_1 <- "NA"
		} else {
			fresh_1 <- ymd(today)-ymd(fresh_1)
		}
		if (fresh_2 == -Inf) {
			fresh_2 <- "NA"
		} else {
			fresh_2 <- ymd(today)-ymd(fresh_2)
		}
		if (fresh_3 == -Inf) {
			fresh_3 <- "NA"
		} else {
			fresh_3 <- ymd(today)-ymd(fresh_3)
		}
		if (fresh_4 == -Inf) {
			fresh_4 <- "NA"
		} else {
			fresh_4 <- ymd(today)-ymd(fresh_4)
		}
		
		# Update the freshness.
		#
		output$rs_fresh_1 <- renderText(fresh_1)
		output$rs_fresh_2 <- renderText(fresh_2)
		output$rs_fresh_3 <- renderText(fresh_3)
		output$rs_fresh_4 <- renderText(fresh_4)
		
		# Update the freshness panel text colors.
		#
		tchange1 <- paste0("document.getElementById('rs_fresh_1').style.color = '", getFreshnessColor(TARGETS_RS[1], GENLOCI_RS[1], loc_ids), "';", sep="")
		runjs(tchange1)
		tchange2 <- paste0("document.getElementById('rs_fresh_2').style.color = '", getFreshnessColor(TARGETS_RS[2], GENLOCI_RS[2], loc_ids), "';", sep="")
		runjs(tchange2)
		tchange3 <- paste0("document.getElementById('rs_fresh_3').style.color = '", getFreshnessColor(TARGETS_RS[3], GENLOCI_RS[3], loc_ids), "';", sep="")
		runjs(tchange3)
		tchange4 <- paste0("document.getElementById('rs_fresh_4').style.color = '", getFreshnessColor(TARGETS_RS[4], GENLOCI_RS[4], loc_ids), "';", sep="")
		runjs(tchange4)

		# Calculate and print the site mean flow.
		#		
		mean_flow <- mean((df_this %>% filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]))$sample_flow)
		output$site_rs_flow <- renderText(
			paste0("Mean daily flow is ", prettyNum(mean_flow, digits = 2), 
						 " MGD with a total reported capacity of ", prettyNum(total_cap, digits = 2),
						 " MGD.", sep="")
		)

	}


	#
	# Render the rs info panels for the first time. 
	# Essentially a copy of updateSiteInfoRS() but I cant' figure out how to get a function to run 
	# outside of a reactive element, or how to wrap it in a reactive element so it runs on launch. 
	# So here we are.
	#

	output$site_rs_title <- renderText(paste0("State of West Virginia", sep=""))

	output$rs_CSL_1 <- renderText(DISEASE_RS[1])
	output$rs_CSL_2 <- renderText(DISEASE_RS[2])
	output$rs_CSL_3 <- renderText(DISEASE_RS[3])
	output$rs_CSL_4 <- renderText(DISEASE_RS[4])
	
	total_popserved <- sum(distinct(df_active_loc %>% filter(location_category == "wwtp"), location_id, location_population_served)$location_population_served)
 	if (total_popserved == -1) {
 		total_popserved = "Unknown"
 	}

	fresh_1 <- max((df_rs %>% filter(target == TARGETS_RS[1] & target_genetic_locus == GENLOCI_RS[1]))$date_primary, na.rm = TRUE)
	fresh_2 <- max((df_rs %>% filter(target == TARGETS_RS[2] & target_genetic_locus == GENLOCI_RS[2]))$date_primary, na.rm = TRUE)
	fresh_3 <- max((df_rs %>% filter(target == TARGETS_RS[3] & target_genetic_locus == GENLOCI_RS[3]))$date_primary, na.rm = TRUE)
	fresh_4 <- max((df_rs %>% filter(target == TARGETS_RS[4] & target_genetic_locus == GENLOCI_RS[4]))$date_primary, na.rm = TRUE)
	
	if (fresh_1 == -Inf) {
		fresh_1 <- "NA"
	} else {
		fresh_1 <- ymd(today)-ymd(fresh_1)
	}
	if (fresh_2 == -Inf) {
		fresh_2 <- "NA"
	} else {
		fresh_2 <- ymd(today)-ymd(fresh_2)
	}
	if (fresh_3 == -Inf) {
		fresh_3 <- "NA"
	} else {
		fresh_3 <- ymd(today)-ymd(fresh_3)
	}
	if (fresh_4 == -Inf) {
		fresh_4 <- "NA"
	} else {
		fresh_4 <- ymd(today)-ymd(fresh_4)
	}
	
	total_cap <- sum(distinct(df_active_wwtp, wwtp_id, wwtp_capacity_mgd)$wwtp_capacity_mgd)+1

	dates <- c(today %m-% months(VIEW_RANGE_PRIMARY), today)
	mean_flow <- mean((df_rs %>% filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]))$sample_flow)
	
	output$rs_fresh_1 <- renderText(fresh_1)
	output$rs_fresh_2 <- renderText(fresh_2)
	output$rs_fresh_3 <- renderText(fresh_3)
	output$rs_fresh_4 <- renderText(fresh_4)
	
	output$site_rs_info <- renderText(
		paste0("Representing ", n_distinct(df_rs$location_id), " active treatment facilities serving ", 
					 prettyNum(total_popserved, big.mark=","), " residents across ", 
					 n_distinct((df_active_loc %>% filter(location_category == "wwtp"))$location_counties_served), " counties.",
					 sep="")
	)

	output$site_rs_flow <- renderText(
		paste0("Mean daily flow is ", prettyNum(mean_flow, digits = 2), 
					 " MGD with a total reported capacity of ", prettyNum(total_cap, digits = 2),
					 " MGD.", sep="")
	)


	# Render the rs plots
	#
	output$plot1_rs_title <- renderText(paste0(TARGETS_RS[1], sep=""))
	output$plot1_rs <- renderPlotly({
		df_plot <- getBasicDF(TARGETS_RS[1], GENLOCI_RS[1], controlRV$mapClick[1], controlRV$viewMonths[1])
		plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	

	output$plot2_rs_title <- renderText(paste0(TARGETS_RS[2], sep=""))
	output$plot2_rs <- renderPlotly({
	  df_plot <- getBasicDF(TARGETS_RS[2], GENLOCI_RS[2], controlRV$mapClick[1], controlRV$viewMonths[1])
	  plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	
	
	output$plot3_rs_title <- renderText(paste0(TARGETS_RS[3], sep=""))
	output$plot3_rs <- renderPlotly({
	  df_plot <- getBasicDF(TARGETS_RS[3], GENLOCI_RS[3], controlRV$mapClick[1], controlRV$viewMonths[1])
	  plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	

	output$plot4_rs_title <- renderText(paste0(TARGETS_RS[4], sep=""))
	output$plot4_rs <- renderPlotly({
	  df_plot <- getBasicDF(TARGETS_RS[4], GENLOCI_RS[4], controlRV$mapClick[1], controlRV$viewMonths[1])
	  plotBasic(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	

	output$plotsq_rs_title <- renderText(paste0("SARS-CoV-2 Variant Proportions", sep=""))
	output$plotsq_rs <- renderPlotly({
		df_plot <- getSeqrDF(TARGETS_RS[3], controlRV$mapClick[1], controlRV$viewMonths[1])
		plotSeqr(df_plot, controlRV$viewMonths[1]) %>% config(displayModeBar = FALSE)# %>% style(hoverinfo = "skip")
	})	


	#
	# Render the maps
	#
	output$map_rs <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircleMarkers(data = df_active_loc %>% filter(location_category == "wwtp"),
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
				fillColor = ~color_group, 
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


	output$map_mg <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircleMarkers(data = df_active_loc %>% filter(location_category == "wwtp"),
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
				fillColor = ~color_group, 
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


	output$map_st <- renderLeaflet({
		leaflet() %>% 
				addTiles() %>% 
				setView(lng = MAP_CENTER$lng, lat = MAP_CENTER$lat, zoom = MAP_CENTER$zoom) %>% 
				addCircleMarkers(data = df_active_loc %>% filter(location_category == "wwtp"),
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
				fillColor = ~color_group, 
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
	# OBSERVER FUNCTIONS, ROUTINE SURVEILLANCE TAB
	#
	###########################

	# 
	# React to map marker click
	#
  observeEvent(input$map_rs_marker_click, { 
		#print("Map MARKER click top")
		
    clickedLocation <- input$map_rs_marker_click$id
		
		loc_id <- unique((df_active_loc %>% filter(location_common_name == clickedLocation))$location_id)

		if (loc_id %in% df_rs$location_id) {
		
			# Update the reactive element
			controlRV$mapClick[1] <- clickedLocation
		
			# Update the plots to include data for clicked marker
		  updatePlotsRS()
		  
		  # Update the site info panel
			updateSiteInfoRS(controlRV$activeGeoLevel[1], controlRV$mapClick[1], controlRV$viewMonths[1])
		  
		} else {
			print(paste0("No data for ", clickedLocation))
		}
				
  }, ignoreInit = TRUE)

	# 
	# React to map shape click
	#
  observeEvent(input$map_rs_shape_click, { 
    clickedLocation <- input$map_rs_shape_click$id
		
		if (clickedLocation %in% df_active_loc$location_counties_served) {

			controlRV$mapClick[1] <- clickedLocation

			# Update the site info panel
			updateSiteInfoRS(controlRV$activeGeoLevel[1], controlRV$mapClick[1], controlRV$viewMonths[1])
		} else {
			# Make some message box appear.
		}
		
  }, ignoreInit = TRUE)

	# 
	# React to map click (off-marker)
	#
#   observeEvent(input$map_rs_click, { 
# 		#print("Map click top")
# 		
# 		# only respond if this click is in a new position on the map
# 		if (input$map_rs_click$lat != controlRV$clickLat | input$map_rs_click$lng != controlRV$clickLng) {
# 			controlRV$clickLat <- 0
# 			controlRV$clickLng <- 0
# 
# 			# Update the reactive element
# 			controlRV$mapClick[1] <- "WV"
# 		
# 			# Update the plots to show state-wide data
# 			updatePlotsRS()
# 			
# 			# Update the site info panel
# 			updateSiteInfoRS(controlRV$activeGeoLevel[1], controlRV$mapClick[1], controlRV$viewMonths[1])
# 		}
# 				  				
#   }, ignoreInit = TRUE)



	#
	# Re-center the map
	#
	observeEvent(input$center_map_rs, {
		
		rsLeafletProxy %>% setView(MAP_CENTER$lng, MAP_CENTER$lat, zoom = MAP_CENTER$zoom)
#		rsLeafletProxy %>% flyTo(map_center$lng, map_center$lat, zoom = map_center$zoom, options = {animate = TRUE})
	}, ignoreInit = TRUE)

	#
	# Change the map active geolayer
	#
  observeEvent(input$geo_level_rs, {

		# Update some reactive elements
		controlRV$activeGeoLevel[1] <- input$geo_level_rs
		
		# Reconfigure the map markers and shapes
		if (tolower(input$geo_level_rs) == "county") {
			rsLeafletProxy %>% 
					clearMarkers() %>% 
					clearShapes() %>% 
					addCircleMarkers(data = df_active_loc %>% filter(location_category == "wwtp"),
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
						fillColor = ~color_group, 
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
			if (tolower(input$geo_level_rs) == "facility") {
				rsLeafletProxy %>% 
						clearMarkers() %>% 
						clearShapes() %>% 
						addPolygons( 
							data = county_spdf, 
							layerId = ~NAME, 
							#fillColor = ~mypalette(color_group), 
							stroke=TRUE, 
							fillOpacity = 0, 
							color="black", 
							weight=1.0, 
							group="county"
						) %>% 
						addCircleMarkers(data = df_active_loc %>% filter(location_category == "wwtp"),
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
			updatePlotsRS()
		}
		
	}, ignoreInit = TRUE)

	#
	# Change the date range to view
	#
  observeEvent(input$view_range_rs, {

		# Update some reactive elements
		controlRV$viewMonths[1] <- as.numeric(input$view_range_rs)
		
		# Update the plots
    updatePlotsRS()

	}, ignoreInit = TRUE)


	# 
	# React to plot click
	#
#  observeEvent(plotww_wide$plotly_click, { 
    #clickedLocation <- input$map_rs_shape_click$id
		#print(plotww_wide$plotly_click)
		
#	}, ignoreInit = TRUE)

})





# WORKING!

# 	# Generate a dataframe for a plot of hospitalization data
# 	getHospPlotData <- function(date_win) {
# 		
# 		dates <- c(today %m-% months(date_win), today)
# 		#years <- c(lubridate::year(dates[1]), lubridate::year(dates[2]))
# 		#weeks <- c(lubridate::week(dates[1]), lubridate::week(dates[2]))
# 		
# 		df_plot <- df_hospital
# 		df_plot$date_to_plot <- as.Date(with(df_plot, paste(mmr_year, mmr_week, 1, sep="-")),"%Y-%U-%u")
# 		
# 		loc_name <- controlRV$mapClick[1]
# 		if (loc_name == "WV") {
# 			df_plot <- df_plot %>% 
# 								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
# 								 #filter(mmr_year >= years[1] & mmr_year <= years[2] & mmr_week >= weeks[1] & mmr_week <= weeks[2]) %>%
# 								 group_by(date_to_plot) %>% 
# 								 arrange(date_to_plot) %>%
# 								 summarize(val := mean(weekly_sum, na.rm = TRUE))
# 		} else {
# 			loc_id <- unique((df_location %>% filter(location_common_name == loc_name))$location_id)
# 		
# 			df_plot <- left_join(
# 					df_targ %>% filter(location_id == loc_id & date_to_plot >= dates[1] & date_to_plot <= dates[2]), 
# 					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
# 					by="location_id")
# 			df_plot <- df_plot %>% 
# 								 group_by(date_to_plot, target_genetic_locus) %>% 
# 								 arrange(date_to_plot, target_genetic_locus) %>%
# 								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
# 		}
# 
# 		df_plot$target <- 'SARS-CoV-2'
# 		return(df_plot)
# 	}
# 
#
# 	# Generate a dataframe for a class plot
# 	getClassPlotData <- function(date_win) {
# 		
# 		targ <- controlRV$visibleTarget
# 		tclass <- controlRV$visibleClass
# 		
# 		df_related <- resources$target %>% filter(target_class == tclass)
# 		#df_targ <- df_rs %>% filter(target %in% df_related$target_id & target != targ)
# 		df_targ <- df_rs %>% filter(target %in% df_related$target_id)
# 		#print(df_targ)
# 		dates <- c(today %m-% months(date_win), today)
# 		
# 		loc_name <- controlRV$mapClick[1]
# 		if (loc_name == "WV") {
# 			df_plot <- df_targ %>% 
# 								 filter(date_to_plot >= dates[1] & date_to_plot <= dates[2]) %>%
# 								 group_by(date_to_plot, target, target_genetic_locus) %>% 
# 								 arrange(date_to_plot, target, target_genetic_locus) %>%
# 								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
# 		} else {
# 			loc_id <- unique((resources$location %>% filter(location_common_name == loc_name))$location_id)
# 		
# 			df_plot <- left_join(
# 					df_targ %>% filter(location_id == loc_id & date_to_plot >= dates[1] & date_to_plot <= dates[2]), 
# 					resources$location %>% filter(tolower(location_status) == "active") %>% select(location_id, location_common_name), 
# 					by="location_id")
# 			df_plot <- df_plot %>% 
# 								 group_by(date_to_plot, target, target_genetic_locus) %>% 
# 								 arrange(date_to_plot, target, target_genetic_locus) %>%
# 								 summarize(val := mean(target_copies_per_l, na.rm = TRUE))
# 		}
# 		return(df_plot)
# 	}
# 	
