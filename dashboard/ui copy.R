shinyUI(fluidPage(
	useShinyjs(),
	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main_styles.css")
  ),
	
	tags$div(class="nbcontainer", 
	
	navbarPage(
		theme = shinytheme("flatly"), 
		id="nav",
		HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Wastewater Testing for Community Health in WV (WaTCH-WV)</a>'), 
		windowTitle = "WaTCH-WV",
		collapsible = FALSE,

		tabPanel(
			"OLD",
			fluidRow(
				column(5, 
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						column(12,
							leafletOutput("map_covid", width = "100%", height = "443px")
						)
					), # fluidRow (map)
					fluidRow(
						column(12,
							div("MAP & PLOT CONTROLS", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
						)
					), # fluidRow (controls title)
					fluidRow(
						column(8,
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:105px;",
								"Map View",
								selectInput(
									"geo_level_covid",
									label = NULL,
									choices = GEOLEVELS, 
									selected = GEOLEVELS_DEFAULT
								)
							),
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:105px;",
								"Map Color",
								selectInput(
									"map_color_covid",
									label = NULL,
									choices = c("By Lab" = "lab"),
			#							choices = c("By Lab" = "lab", "FluA" = "FluA", "SARS-CoV-2" = "SARS-CoV-2", "RSV" = "RSV"),
									selected = "lab"
								)
							),
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:120px;",
								"Plot View",
								selectInput(
									"view_range_covid",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						),
						column(4,
							div(
								style = "display: inline-block;margin-top: 25px;",
								actionBttn(inputId="center_map_rs", label="Reset Map", style="pill", size="xs", color="success")
							),
							div(
								class = "logo",
								style = "display: inline-block;margin-left: 15px;",
								tags$a(href='https://www.watch-wv.com/', 
								tags$img(src='WaTCH-WV_logo.png',height='50',width='50'))
							)
						)
					), # fluidRow (controls)
					fluidRow(
						column(12,
							style = "margin-top: 5px;",
							# SARS seqr data
							div(textOutput("plotsq_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plotsq_rs", height="325px", width="100%")
						)
					) # fluidRow (variant plot)
				),
				column(7,
					style = "margin-top: 5px; margin-bottom: 5px;",
					fluidRow(
						column(6,
							#style = "padding: 3px;",
							fluidRow(
								#style = "border: 2px solid #941100;",
								column(12,
									# title
									style = "text-align: center;font-size: 18px;font-weight: 800;padding: 3px;",
									div(
										textOutput("site_rs_title"), 
										style="padding: 4px;color: #ffffff;background-color: #303D4E;"
									)
								)
							), # fluidRow (selection title)
							fluidRow(
								# alert blocks
								#style = "border: 2px solid #941100;",
								column(4,
									style = "text-align: right;font-size: 13px;font-weight: 800;padding: 3px;",
									div(id = "rs_CSL_label", "Alert Color Level: ", style = "padding: 2px;height: 25px;"),
									div(id = "rs_delta_label", paste0("% of ", VIEW_RANGE_PRIMARY, " Mo. Mean:"), style = "padding: 2px;height: 25px;"),
									div(id = "rs_fresh_label", "Data Age (days):", style = "padding: 2px;height: 25px;")
								),
								column(2,
									style = "text-align: center;font-size: 14px;font-weight: 800;padding: 3px;",
									div(id = "rs_CSL_1", textOutput("rs_CSL_1"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_delta_1", textOutput("rs_delta_1"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_fresh_1", textOutput("rs_fresh_1"), style = "height: 25px;padding: 2px;")
								),
								column(2,
									style = "text-align: center;font-size: 14px;font-weight: 800;padding: 3px;",
									div(id = "rs_CSL_2", textOutput("rs_CSL_2"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_delta_2", textOutput("rs_delta_2"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_fresh_2", textOutput("rs_fresh_2"), style = "height: 25px;padding: 2px;")
								),
								column(2,
									style = "text-align: center;font-size: 14px;font-weight: 800;padding: 3px;",
									div(id = "rs_CSL_3", textOutput("rs_CSL_3"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_delta_3", textOutput("rs_delta_3"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_fresh_3", textOutput("rs_fresh_3"), style = "height: 25px;padding: 2px;")
								),
								column(2,
									style = "text-align: center;font-size: 14px;font-weight: 800;padding: 3px;",
									div(id = "rs_CSL_4", textOutput("rs_CSL_4"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_delta_4", textOutput("rs_delta_4"), style = "height: 25px;padding: 2px;"),
									div(id = "rs_fresh_4", textOutput("rs_fresh_4"), style = "height: 25px;padding: 2px;")
								)
							), # fluidRow (alert blocks)
							fluidRow(
								# alert key & site info
								style = "background-color: #000000;color: #ffffff;text-align: center;margin-top: 10px;padding: 5px;",
								column(5,
									div(
										id = "alert_scale", 
										#class = "panel panel-default",
										style = "margin 5px;text-align: center;",
										div("Alert Level Key", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;width:125px;"),
										div("LOW", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #3288BD;"),
										div("MODERATE", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #E6F598;"),
										div("HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #FDAE61;"),
										div("VERY HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #D53E4F;"),
										div("UNKNOWN", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #EEEEEE;")
									)
								),
								column(7,
									style = "padding: 3px;text-align: center;font-size: 13px; font-weight: 400;",
									div(
										div(textOutput("site_rs_info"), style="padding: 8px 0px;"),
										div(textOutput("site_rs_flow"), style="padding: 2px 0px;")
									)
								)	# column
							) # fluidRow (color key & info)
						),
						column(6,
							fluidRow(
								column(12,
									# title
									style = "padding: 3px 0px 0px 0px;text-align: center;font-size: 13px;font-weight: 800;",
									div(
										paste0("Percent of ", VIEW_RANGE_PRIMARY, " Month Mean as of ", today, " (By County)"), 
										style = "background-color: #000000;color: #FFFFFF;"
									)
								) # column
							), # fluidRow (county table title)
							fluidRow(
								# alert blocks
								column(4,
									style = "background-color: #000000;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
									div("County", style = "padding: 2px;")
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: right;font-size: 15px;font-weight: 800;padding: 3px;",
									div(paste0(DISEASE_RS[1]), style = "padding: 2px;"),
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: center;font-size: 15px;font-weight: 800;padding: 3px;",
									div(paste0(DISEASE_RS[2]), style = "padding: 2px;"),
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: center;font-size: 15px;font-weight: 800;padding: 3px;",
									div(paste0(DISEASE_RS[3]), style = "padding: 2px;"),
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
									div(paste0(DISEASE_RS[4]), style = "padding: 2px;")
								),
							), # fluidRow (table header)
							fluidRow(
								column(12,
									style = "overflow: auto; padding: 3px;height: 236px;",
									tableHTML(
										df_regions %>% filter(region_geolevel == "county") %>% select(region_name, any_of(DISEASE_RS)) %>% rename(County = region_name), 
										#collapse = "separate_shiny", 
										#spacing = "5px 2px", 
										rownames = FALSE, 
										border = 0
										#widths = c(96, 60, 60, 65, 60)
									) %>% 
									add_css_thead(css = list("background-color", "#000000")) %>% 
									add_css_thead(css = list("color", "#000000")) %>% 
									add_css_thead(css = list("font-size", "1px")) %>% 
									add_css_column(css = list("text-align", "left"), columns=names(df_regions[2:5])) %>% 
									add_css_table(css = list("width", "100%")) %>% 
									add_css_table(css = list("background-color", "#ffffff")) %>% 
									add_css_row(css = list("background-color", "#f2f2f2"), rows = odd(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1)) %>%
									add_css_row(css = list("background-color", "#e6f0ff"), rows = even(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1))
								) # column
							) # fluidRow (table)
						)
					), # fluidRow (selection info)
					fluidRow(
						style = "margin-top: 5px;",
						column(6,
							# target 1 plot (fluA)
							style = "padding: 3px;",
							div(textOutput("plot1_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plot1_rs", height="210px")
						),
						column(6,
							# target 2 plot (fluB)
							style = "padding: 3px;",
							div(textOutput("plot2_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plot2_rs", height="210px")
						)
					), # fluidRow (T1 & T2)
					fluidRow(
						style = "margin-top: 5px;",
						column(6,
							# target 3 plot (COVID)
							style = "padding: 3px;",
							div(textOutput("plot3_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plot3_rs", height="210px")
						),
						column(6,
							# target 4 plot (RSV)
							style = "padding: 3px;",
							div(textOutput("plot4_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plot4_rs", height="210px")
						)
					), # fluidRow (T3 & T4)
					fluidRow(
						style = "background-color: #000000;color: #ffffff;",
						column(7,
							# explanation of plots
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
								style = "padding: 2px 20px 2px 2px;text-align: left;",
								"Plot values are reported as the copies of viral particles per person after correction for daily flow, averaged across all selected sites. SARS-CoV-2 variant proportions are shown as the percent of total identified variants, averaged across all selected sites."
							)
						),
						column(5,
							# explanation of trend lines
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
								style = "padding: 2px;text-align: left;",
								span("The "),
								span(style="color: #EAAA00;", "solid gold line"),
								span(paste0("on each plot represents the ", VIEW_RANGE_PRIMARY, " month average level of each target, used to generate the percent change. The ")),
								span(style="color: #00B140;", "green dashed line"),
								span("represents the average over the most recent month.")
							)
						)
					) # fluidRow (footnotes)
				)
			),
			
			hidden(
				absolutePanel(
					id = "alert_scale_info",
					class = "mdinfo",
					top = 365, left = 610, width = 580, height = 350,
					div("Alert Color Explanations", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "alertinfo",
						id = "level_1",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
						span(paste0("CODE BLUE. The latest amount of this disease agent is less than 50% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be low."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_2",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
						span(paste0("CODE YELLOW. The latest amount of this disease agent is between 50% and 100% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be moderate."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_3",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
						span(paste0("CODE ORANGE. The latest amount of this disease agent is between 100% and 150% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_4",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
						span(paste0("CODE RED. The latest amount of this disease agent is greater than 150% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be very high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_5",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
						span("The data for this disease agent is either missing or too old to make an accurate determination.", style="font-size: 13px;")
					),
					div(
						style="padding-top: 15px; padding-right: 5px; float: right;",
						actionBttn(inputId="alert_scale_info_close", label="Close", style="pill", size="xs", color="success")
					) # button div
				)
			) # hidden (alert level info)

		), # tabPanel (routine surveillance)

		tabPanel(
			"COVID",
			fluidRow(
				column(5, 
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						column(12,
							leafletOutput("map_covid", width = "100%", height = "443px")
						)
					), # fluidRow (map)
					fluidRow(
						column(12,
							div("MAP & PLOT CONTROLS", 
							style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
						)
					), # fluidRow (controls title)
					fluidRow(
						column(8,
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:105px;",
								"Map View",
								selectInput(
									"geo_level_covid",
									label = NULL,
									choices = GEOLEVELS, 
									selected = GEOLEVELS_DEFAULT
								)
							),
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:105px;",
								"Map Color",
								selectInput(
									"map_color_covid",
									label = NULL,
									choices = c("Lab" = "lab"),
			#							choices = c("Status" = "Status", "Risk" = "Risk", "Freshness" = "Freshness"),
									selected = "lab"
								)
							),
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:120px;",
								"Plot View",
								selectInput(
									"view_range_covid",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						),
						column(4,
							div(
								style = "display: inline-block;margin-top: 25px;",
								actionBttn(inputId="map_center_covid", label="Reset Map", style="pill", size="xs", color="success")
							),
							div(
								class = "logo",
								style = "display: inline-block;margin-left: 15px;",
								tags$a(href='https://www.watch-wv.com/', 
								tags$img(src='WaTCH-WV_logo.png',height='50',width='50'))
							)
						)
					), # fluidRow (controls)
					fluidRow(
						#style = "border: 2px solid #941100;",
						column(12,
							# title
							div(
								"State of West Virginia",
								#textOutput("site_title"), 
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E"
							)
						)
					), # fluidRow (selection title)
					fluidRow(
						# alert key & site info
						style = "background-color: #000000;color: #fffff;",
						column(5,
							div(
								id = "alert_scale", 
								#class = "panel panel-default",
								style = "margin 5px;text-align: center;",
								div("Alert Level Key", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;width:125px;"),
								div("LOW", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #3288BD;"),
								div("MODERATE", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #E6F598;"),
								div("HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #FDAE61;"),
								div("VERY HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #D53E4F;"),
								div("UNKNOWN", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:125px;color: #000000;background-color: #EEEEEE;")
							)
						),
						column(7,
							style = "padding: 3px;text-align: center;font-size: 13px; font-weight: 400;",
							div(
								div(
									"Info on the currently selected site or county appears here.",
									#textOutput("site_info"), 
									style="padding: 8px 0px;"
								),
								div(
									"Flow for the currently selected site(s) appears here.",
									#textOutput("site_flow"), 
									style="padding: 2px 0px;"
								)
							)
						)	# column
					) # fluidRow (color key & selection info)
# 					fluidRow(
# 						style = "background-color: #000000;color: #ffffff;",
# 						column(7,
# 							# explanation of plots
# 							style = "font-size: 14px;padding: 3px;font-weight: 400;",
# 							div(
# 								style = "padding: 2px 20px 2px 2px;text-align: left;",
# 								"Plot values are reported as the copies of viral particles per person after correction for daily flow, averaged across all selected sites. SARS-CoV-2 variant proportions are shown as the percent of total identified variants, averaged across all selected sites."
# 							)
# 						),
# 						column(5,
# 							# explanation of trend lines
# 							style = "font-size: 14px;padding: 3px;font-weight: 400;",
# 							div(
# 								style = "padding: 2px;text-align: left;",
# 								span("The "),
# 								span(style="color: #EAAA00;", "solid gold line"),
# 								span(paste0("on each plot represents the ", VIEW_RANGE_PRIMARY, " month average level of each target, used to generate the percent change. The ")),
# 								span(style="color: #00B140;", "green dashed line"),
# 								span("represents the average over the most recent month.")
# 							)
# 						)
# 					) # fluidRow (footnotes)
				),
				column(5,
					style = "margin-top: 5px;padding: 0px;",
					fluidRow(
						# alert blocks
						column(3,
							fluidRow(
								column(12,
									div("TREND", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
								)
							) # fluidRow (trend title)
						),
						column(3,
							fluidRow(
								column(12,
									div("LEVEL", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
								)
							) # fluidRow (level title)
						),
						column(1,
							fluidRow(
								column(12,
									div("FRESH", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
								)
							) # fluidRow (freshness title)
						),
						column(4,
							fluidRow(
								column(12,
									div("VARIANT", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
								)
							) # fluidRow (variant title)
						),
						column(1,
							fluidRow(
								column(12,
									div("FRESH", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
								)
							) # fluidRow (sequence freshness title)
						),
					), # fluidRow (alert blocks)
					fluidRow(
						style = "margin-top: 5px;",
						column(12,
							# Plot of COVID change over time
							style = "padding: 3px;",
							div("SARS-CoV-2 in Wastewater", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plot_covid", height="210px")
						)
					), # fluidRow (COVID plot)
					fluidRow(
						column(12,
							style = "margin-top: 5px;",
							# SARS seqr data
							div(textOutput("plotsq_covid_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
							plotlyOutput("plotsq_covid", height="325px", width="100%")
						)
					) # fluidRow (variant plot)
				),
				column(2,
					style = "margin-top: 5px;padding: 3px;",
							fluidRow(
								column(12,
									# title
									style = "padding: 3px 0px 0px 0px;text-align: center;font-size: 13px;font-weight: 800;",
									div(
										paste0("Percent of ", VIEW_RANGE_PRIMARY, " Month Mean as of ", today, " (By County)"), 
										style = "background-color: #000000;color: #FFFFFF;"
									)
								) # column
							), # fluidRow (county table title)
							fluidRow(
								column(6,
									style = "background-color: #000000;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
									div("County", style = "padding: 2px;")
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: right;font-size: 15px;font-weight: 800;padding: 3px;",
									div("Level", style = "padding: 2px;"),
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: center;font-size: 15px;font-weight: 800;padding: 3px;",
									div("Trend", style = "padding: 2px;"),
								),
								column(2,
									style = "background-color: #000000;color: #FFFFFF;text-align: center;font-size: 15px;font-weight: 800;padding: 3px;",
									div("Freshness", style = "padding: 2px;"),
								)
							), # fluidRow (table header)
							fluidRow(
								column(12,
									style = "overflow: auto; padding: 3px;height: 236px;",
									tableHTML(
										df_regions %>% filter(region_geolevel == "county") %>% select(region_name, any_of(DISEASE_RS)) %>% rename(County = region_name), 
										#collapse = "separate_shiny", 
										#spacing = "5px 2px", 
										rownames = FALSE, 
										border = 0
										#widths = c(96, 60, 60, 65, 60)
									) %>% 
									add_css_thead(css = list("background-color", "#000000")) %>% 
									add_css_thead(css = list("color", "#000000")) %>% 
									add_css_thead(css = list("font-size", "1px")) %>% 
									add_css_column(css = list("text-align", "left"), columns=names(df_regions[2:5])) %>% 
									add_css_table(css = list("width", "100%")) %>% 
									add_css_table(css = list("background-color", "#ffffff")) %>% 
									add_css_row(css = list("background-color", "#f2f2f2"), rows = odd(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1)) %>%
									add_css_row(css = list("background-color", "#e6f0ff"), rows = even(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1))
								) # column
							) # fluidRow (table)
				)
			)
		), # tabPanel (COVID)

		tabPanel(
			"Outbreaks",
			div(
				class="outer",
				leafletOutput("map_st", width=645, height=457),

				absolutePanel(
					class = "logo", 
					top = 380, left = 535, width = 80, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					tags$a(href='https://www.watch-wv.com/', 
					tags$img(src='WaTCH-WV_logo.png',height='80',width='80'))
				), # absolutePanel
				
				absolutePanel(
					class = "map-recenter", 
					top = 460, left = 530, width = 125, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					actionBttn(inputId="center_map_st", label="Recenter Map", style="pill", size="xs", color="success")
				), # absolutePanel
				
				absolutePanel(
					#class = "panel panel-default",
					class = "controls",
					top = 70, left = 260, height = 160, width = 450,
					fixed=TRUE, draggable=FALSE,
					div(
						class = "map_embed",
						selectInput(
							"geo_level_st",
							label = NULL,
							#width = "200px",
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						),
						selectInput(
							"view_range_st",
							label = NULL,
							#width = "60px",
							choices = VIEW_RANGES, 
							selected = VIEW_RANGE_PRIMARY
						)
					)
				), # absolutePanel

			) # tabPanel div
		) # tabPanel (outbreaks)

	)) # navbarPage and container div
))

