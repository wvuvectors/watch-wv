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
			"COVID",
			fluidRow(
				column(12,
					# title
					style = "border: 0px solid #000000;margin-left: 0px; margin-right: 0px;",
					div(
						textOutput("selection_title_covid"), 
						style="font-size: 24px;padding-top: 2px;padding-bottom: 2px;font-weight: 800;margin: 0px;text-align: center;color: #ffffff; background-color: #008F00"
					)
				)
			), # fluidRow (selection title)
			fluidRow(
				column(5, 
					style = "margin-top: 5px;",
					fluidRow(
						# alert blocks
						style = "margin-top: 0px;", 
						column(4,
							div(
								id = "abundance_title_covid",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_covid",
								textOutput("abundance_covid"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(4,
							div(
								id = "trend_title_covid",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_covid",
								textOutput("trend_covid"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(4,
							div(
								id = "variant_title_covid",
								"Variant",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "variant_text_covid",
								textOutput("variant_covid"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						)
					), # fluidRow (alert blocks)
					fluidRow(
						# Alert details
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFE7E5;color: #000000;",
						column(12,
							div(
								textOutput("alert_details_covid"), 
								style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert details)
					fluidRow(
						# Selection details
						style = "margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_covid"), 
								style="font-size: 16px;font-weight: 800;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
# 					fluidRow(
# 						style = "background-color: #fbfbfb; margin-left: 0px;margin-right: 0px;margin-bottom: 2px;", 
# 						column(12,
# 							style = "text-align: right",
# 							div(
# 								style = "display: inline-block;margin-top: 0px;",
# 								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success")
# 							)
# 						)
#						column(6,
# 							div(
# 								class = "map_embed",
# 								style = "display: inline-block;font-size: 12px;font-weight: 400;text-align: center; width:120px;",
# 								"Color By Amount of:",
# 								selectInput(
# 									"map_color",
# 									label = NULL,
# 									choices = c("COVID"),
# 									selected = "COVID"
# 								)
# 							),
# 							div(
# 								class = "map_embed",
# 								style = "margin-left: 5px; display: inline-block;font-size: 12px;font-weight: 400;text-align: center; width:120px;",
# 								"Plot Most Recent:",
# 								selectInput(
# 									"view_range",
# 									label = NULL,
# 									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
# 									selected = VIEW_RANGE_PRIMARY
# 								)
# 							)
# 						)
#					), # fluidRow (controls)
					fluidRow(
						column(12,
							leafletOutput("map_covid", width = "100%", height = "470px")
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								"Click any county on the map below to see specific data for that region. Click anywhere else on the map to return to statewide results. ",
								"Map counties are colored according to the amount (abundance) of the disease present in the latest reporting period. ", 
								"IMPORTANT: Data is subject to change as additional sites report results.", 
								style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
							)
						)	# column
					) # fluidRow (footnotes)
				),
				column(7,
					style = "margin-top: 5px;",
					fluidRow(
						column(12,
							div(
								textOutput("plot_title_covid"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (COVID plot title)
					fluidRow(
						column(12,
							div(
								textOutput("selection_freshness_covid"), 
								style="color: #000000;background-color: #fbfbfb;font-size: 16px;padding-top: 5px;font-weight: 800;text-align: center;"
							),
							div(
								textOutput("selection_completeness_covid"), 
								style="color: #000000;background-color: #fbfbfb;font-size: 16px;padding-bottom: 5px;font-weight: 800;text-align: center;"
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						column(12,
							# Plot of COVID change over time
							plotlyOutput("plot_covid", height="350px", width="100%")
						)
					), # fluidRow (COVID plot)
					fluidRow(
						#style = "margin-left: 0px; margin-right: 0px;", 
						style = "margin-top: 4px;", 
						column(12,
							div(
								textOutput("plotsq_title_covid"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (COVID variant plot title)
					fluidRow(
						column(12,
							# Plot of COVID variant proportions over time
							plotlyOutput("plotsq_covid", height="320px", width="100%")
						)
					) # fluidRow (variant plot)
				)
# 				column(2,
# 					style = "margin-top: 8px;",
# 					fluidRow(
# 						# Data info
# 						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
# 						column(12,
# 							fluidRow(
# 								column(12,
# 									div(
# 									"Latest variant data for this region is for ",
# 										style="font-size: 18px;padding-top: 10px;font-weight: 800;text-align: center;color: #000000;"
# 									),
# 									div(
# 										textOutput("selectionsq_freshness_covid"), 
# 										style="font-size: 24px;padding-bottom: 0px;font-weight: 800;text-align: center;color: #9437FF;"
# 									),
# 									div(
# 										textOutput("selectionsq_completeness_covid"), 
# 										style="font-size: 16px;padding-bottom: 10px;font-weight: 400;text-align: center;color: #000000;"
# 									)
# 								)
# 							), # fluidRow (variant data freshness)
# 							fluidRow(
# 								column(12,
# 									div(
# 									"NOTE: Recent data is subject to change as more sites report results.",
# 										style="font-size: 14px;font-style: italic;padding-top: 10px; padding-bottom: 5px;font-weight: 400;text-align: center;color: #333333"
# 									)
# 								)
# 							) # fluidRow (about the data)
# 						)	# column
# 					), # fluidRow (data info)
# 					fluidRow(
# 						style = "margin-top: 12px; background-color: #f3f3e1; color: #000000; padding: 5px;",
# 						column(12,
# 							# explanation of plots
# 							style = "font-size: 14px;padding: 3px;font-weight: 400;",
# 							div(
# #								style = "padding: 2px 20px 2px 2px;text-align: left;",
# 								"Plot values are reported as the relative abundance of pathogen after correction for daily flow and population served, averaged across all selected sites."
# 							)
# 						)
# 					), # fluidRow (abundance plot explanation)
# 					fluidRow(
# 						style = "margin-top: 12px; background-color: #000000; color: #ffffff; padding: 5px;",
# 						column(12,
# 							# explanation of trend lines
# 							style = "font-size: 14px;padding: 3px;font-weight: 400;",
# 							div(
# #								style = "padding: 2px;text-align: left;",
# 								span("The "),
# 								span(style="color: #00B140;", "green "),
# 								span("and "),
# 								span(style="color: #EAAA00;", "gold "),
# 								span("lines represent the average abundance over the most recent "),
# 								span(style="color: #00B140;", "3 months "),
# 								span("and "),
# 								span(style="color: #EAAA00;", "12 months "),
# 								span("respectively."),
# 							)
# 						)
# 					), # fluidRow (trend line descriptions)
# 					fluidRow(
# 						style = "margin-top: 12px; background-color: #f3f3e1; color: #000000; padding: 5px;",
# 						column(12,
# 							# explanation of plots
# 							style = "font-size: 14px;padding: 3px;font-weight: 400;",
# 							div(
# #								style = "padding: 2px 20px 2px 2px;text-align: left;",
# 								"Variant proportions are shown as the percent of total identified variants and averaged across all selected sites."
# 							)
# 						)
# 					), # fluidRow (variant plot explanation)
# 					fluidRow(
# 						style = "margin-top: 12px; background-color: #ffffff; color: #ffffff; padding: 5px;",
# 						column(12,
# 							#style = "font-size: 14px;padding: 3px;font-weight: 400;",
# 							#div("WaTCH-WV is supported by CDC-sponsored grants from the West Virginia Department of Health to West Virginia University & Marshall University.")
# 							div(
# 								style = "margin-top: 20px;text-align: center;",
# 								downloadBttn(outputId="download_data_covid", label="Download COVID Data", style="pill", size="s", color="royal")
# 							)
# 						)
# 					) # fluidRow (data download)
# 				)
			),

			hidden(
				absolutePanel(
					id = "missing_data_p1_popup",
					class = "mdinfo",
					top = 320, left = 710, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			), # hidden (missing data p1 popup)
			hidden(
				absolutePanel(
					id = "missing_data_p0_popup",
					class = "mdinfo",
					top = 650, left = 690, width = 450, height = 100,
					div("There is no variant data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			) # hidden (missing data p0 popup)

		), # tabPanel (COVID)

		tabPanel(
			"Influenza",
			fluidRow(
				column(5, 
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						#style = "border: 2px solid #941100;",
						column(12,
							# title
							div(
								textOutput("selection_title_flu"), 
								style="font-size: 24px;padding: 6px;margin-bottom: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #008F00"
							)
						)
					), # fluidRow (selection title)
					fluidRow(
						# Selection details
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							fluidRow(
								column(12,
									div(
										textOutput("selection_details_flu"), 
										style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
									),
									div(
									"Click any county on the map below to see specific data for that region. Click anywhere else on the map to return to statewide results.",
										style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
									),
								)
							) # fluidRow (selection details)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						style = "background-color: #F3EFEA; margin-left: 0px;margin-right: 0px;",
						column(6,
							style = "text-align: left",
							div(
								style = "display: inline-block;margin-top: 20px;",
								actionBttn(inputId="map_reset", label="Reset Zoom", style="pill", size="xs", color="success")
							)
						), 
						column(6,
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 12px;font-weight: 400;text-align: center; width:120px;",
								"Color By Amount of:",
								selectInput(
									"map_color",
									label = NULL,
									choices = c("FluA", "FluB"),
									selected = "FluA"
								)
							),
							div(
								class = "map_embed",
								style = "margin-left: 5px; display: inline-block;font-size: 12px;font-weight: 400;text-align: center; width:120px;",
								"Plot Most Recent:",
								selectInput(
									"view_range",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (controls)
					fluidRow(
						column(12,
							leafletOutput("map_flu", width = "100%", height = "420px")
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							fluidRow(
								column(12,
									div(
										"It is now ", printy_dates(today), ".", 
										style="font-size: 18px;padding-top: 10px;font-weight: 400;text-align: center;color: #333333;"
									)
								)
							),
							fluidRow(
								column(12,
									div(
									"Latest Flu A abundance data for this region is from ",
										style="font-size: 18px;padding-top: 10px;font-weight: 800;text-align: center;color: #000000;"
									),
									div(
										textOutput("selection_freshness_flua"), 
										style="font-size: 24px;padding: 0px;font-weight: 800;text-align: center;color: #ee3333;"
									),
									div(
										textOutput("selection_completeness_flua"), 
										style="font-size: 16px;padding-bottom: 10px;font-weight: 400;text-align: center;color: #000000;"
									)
								)
							), # fluidRow (fluA abundance data freshness)
							fluidRow(
								column(12,
									div(
									"Latest Flu B abundance data for this region is from ",
										style="font-size: 18px;padding-top: 10px;font-weight: 800;text-align: center;color: #000000;"
									),
									div(
										textOutput("selection_freshness_flub"), 
										style="font-size: 24px;padding-bottom: 0px;font-weight: 800;text-align: center;color: #9437FF;"
									),
									div(
										textOutput("selection_completeness_flub"), 
										style="font-size: 16px;padding-bottom: 10px;font-weight: 400;text-align: center;color: #000000;"
									)
								)
							), # fluidRow (fluB abundance data freshness)
							fluidRow(
								column(12,
									div(
									"NOTE: Recent data is subject to change as more sites report results.",
										style="font-size: 14px;font-style: italic;padding-top: 10px; padding-bottom: 5px;font-weight: 400;text-align: center;color: #333333"
									)
								)
							) # fluidRow (about the data)
						)	# column
					), # fluidRow (data info)
				),
				column(5,
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						#style = "margin-left: 0px; margin-right: 0px;", 
						column(12,
							div(
								textOutput("plot_title_flua"), 
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (fluA plot title)
					fluidRow(
						style = "margin-top: 5px;", 
						column(6,
							div(
								id = "abundance_title_flua",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_flua",
								textOutput("abundance_flua"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(6,
							div(
								id = "trend_title_flua",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_flua",
								textOutput("trend_flua"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						)
					), # fluidRow (fluA status blocks)
					fluidRow(
						column(12,
							# Plot of fluA change over time
							plotlyOutput("plot_flua", height="350px", width="100%")
						)
					), # fluidRow (fluA plot)
					fluidRow(
						#style = "margin-left: 0px; margin-right: 0px;", 
						column(12,
							div(
								textOutput("plot_title_flub"), 
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (fluB plot title)
					fluidRow(
						style = "margin-top: 5px;", 
						column(6,
							div(
								id = "abundance_title_flub",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_flub",
								textOutput("abundance_flub"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(6,
							div(
								id = "trend_title_flub",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_flub",
								textOutput("trend_flub"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						)
					), # fluidRow (fluB status blocks)
					fluidRow(
						column(12,
							# Plot of fluB change over time
							plotlyOutput("plot_flub", height="350px", width="100%")
						)
					) # fluidRow (fluB plot)
				),
				column(2,
					style = "margin-top: 8px;",
					fluidRow(
						column(12,
						style = "padding-top: 4px;padding-bottom: 4px;margin-left: 0px; margin-right: 0px;background-color: #000000;color: #ffffff;",
							div(
								class = "logo",
								id = "abundance_level_key_flu", 
								#class = "panel panel-default",
								style = "margin 5px;padding: 3px;text-align: center;",
								div("Abundance Level Key", style="font-size: 14px;font-weight: 800;text-align: center;padding: 2px;"),
								div("VERY HIGH", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #D53E4F;"),
								div("HIGH", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #FDAE61;"),
								div("MODERATE", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #E6F598;"),
								div("LOW", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #3288BD;"),
								div("UNKNOWN", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #EEEEEE;")
							)
						) # column
					), # fluidRow (abundance level color key)
					fluidRow(
						column(12,
						style = "padding-top: 4px;padding-bottom: 4px;margin-left: 0px; margin-right: 0px;margin-top: 10px;margin-bottom: 20px;background-color: #000000;color: #ffffff;",
							div(
								class = "logo",
								id = "trend_key_flu", 
								#class = "panel panel-default",
								style = "margin 5px;padding: 3px;text-align: center;",
								div("Trend Level Key", style="font-size: 14px;font-weight: 800;text-align: center;padding: 2px;"),
								div("SPIKING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #FF85FF;"),
								div("DESPIKING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #96F786;"),
								div("INCREASING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #D53E4F;"),
								div("VARIABLE", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #FDAE61;"),
								div("STABLE", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #E6F598;"),
								div("DECREASING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #3288BD;"),
								div("UNKNOWN", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #EEEEEE;")
							)
						) #  column
					), # fluidRow (trend color key)
					fluidRow(
						style = "margin-top: 12px; background-color: #f3f3e1; color: #000000; padding: 5px;",
						column(12,
							# explanation of plots
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px 20px 2px 2px;text-align: left;",
								"Plot values are reported as the relative abundance of pathogen after correction for daily flow and population served, averaged across all selected sites."
							)
						)
					), # fluidRow (abundance plot explanation)
					fluidRow(
						style = "margin-top: 12px; background-color: #000000; color: #ffffff; padding: 5px;",
						column(12,
							# explanation of trend lines
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px;text-align: left;",
								span("The "),
								span(style="color: #00B140;", "green "),
								span("and "),
								span(style="color: #EAAA00;", "gold "),
								span("lines represent the average abundance over the most recent "),
								span(style="color: #00B140;", "3 months "),
								span("and "),
								span(style="color: #EAAA00;", "12 months "),
								span("respectively."),
							)
						)
					), # fluidRow (trend line descriptions)
					fluidRow(
						style = "margin-top: 12px; background-color: #ffffff; color: #ffffff; padding: 5px;",
						column(12,
							#style = "font-size: 14px;padding: 3px;font-weight: 400;",
							#div("WaTCH-WV is supported by CDC-sponsored grants from the West Virginia Department of Health to West Virginia University & Marshall University.")
							div(
								style = "margin-top: 20px;text-align: center;",
								downloadBttn(outputId="download_data_flua", label="Download FluA Data", style="pill", size="s", color="royal")
							),
							div(
								style = "margin-top: 20px;text-align: center;",
								downloadBttn(outputId="download_data_flub", label="Download FluB Data", style="pill", size="s", color="primary")
							)
						)
					) # fluidRow (data download)
				)
			),

			hidden(
				absolutePanel(
					id = "abundance_level_key_flu_popup",
					class = "mdinfo",
					top = 365, left = 710, width = 580, height = 350,
					div("Abundance Level Color Key", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "alertinfo",
						id = "level_4",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
						span(paste0("VERY HIGH. The latest amount of this disease agent is greater than 150% of the 3 month average. Community transmission is estimated to be very high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_3",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
						span(paste0("HIGH. The latest amount of this disease agent is between 100% and 150% of the 3 month average. Community transmission is estimated to be high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_2",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
						span(paste0("MODERATE. The latest amount of this disease agent is between 50% and 100% of the 3 month average. Community transmission is estimated to be moderate."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_1",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
						span(paste0("LOW. The latest amount of this disease agent is less than 50% of the 3 month average. Community transmission is estimated to be low."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_5",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
						span("UNKNOWN. The data for this disease agent is either missing or too old to make an accurate determination.", style="font-size: 13px;")
					),
					div(
						style="padding-top: 15px; padding-right: 5px; float: right;",
						actionBttn(inputId="abundance_level_key_flu_popup_close", label="Close", style="pill", size="xs", color="success")
					) # button div
				)
			), # hidden (abundance level key popup)

			hidden(
				absolutePanel(
					id = "trend_key_flu_popup",
					class = "mdinfo",
					top = 265, left = 710, width = 580, height = 480,
					div("Trend Color Key", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "alertinfo",
						id = "level_5",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FF85FF; background-color: #FF85FF; border: 1px solid black; border-radius: 3px;"),
						span("SPIKING. The latest abundance of this disease agent is more than 500% higher than the previous value.", style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_6",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #96F786; background-color: #96F786; border: 1px solid black; border-radius: 3px;"),
						span("DESPIKING. The latest abundance of this disease agent has dropped more than 500% from than the previous value.", style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_4",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
						span(paste0("INCREASING. The abundance of this disease agent has risen more than it has declined over the most recent period.."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_3",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
						span(paste0("VARIABLE. The abundance of this disease agent has been variable over the most recent period."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_2",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
						span(paste0("STABLE. The abundance of this disease agent has been relatively stable over the most recent period."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_1",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
						span(paste0("DECREASING. The abundance of this disease agent has generally declined over the most recent period (usually 4 weeks)."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_7",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
						span("INDETERMINATE. There is not enough data for this disease agent to make an accurate determination.", style="font-size: 13px;")
					),
					div(
						style="padding-top: 15px; padding-right: 5px; float: right;",
						actionBttn(inputId="trend_key_flu_popup_close", label="Close", style="pill", size="xs", color="success")
					) # button div
				)
			), # hidden (trend key popup)

			hidden(
				absolutePanel(
					id = "missing_data_p2_popup",
					class = "mdinfo",
					top = 320, left = 710, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			), # hidden (missing data p3 popup)

			hidden(
				absolutePanel(
					id = "missing_data_p3_popup",
					class = "mdinfo",
					top = 650, left = 690, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			) # hidden (missing data p4 popup)

		), # tabPanel (FLU)

		tabPanel(
			"RSV",
			fluidRow(
				column(5, 
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						#style = "border: 2px solid #941100;",
						column(12,
							# title
							div(
								textOutput("selection_title_rsv"), 
								style="font-size: 24px;padding: 6px;margin-bottom: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #008F00"
							)
						)
					), # fluidRow (selection title)
					fluidRow(
						# Selection details
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							fluidRow(
								column(12,
									div(
										textOutput("selection_details_rsv"), 
										style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
									),
									div(
									"Click any county on the map below to see specific data for that region. Click anywhere else on the map to return to statewide results.",
										style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
									),
								)
							) # fluidRow (selection details)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						style = "background-color: #F3EFEA; margin-left: 0px;margin-right: 0px;",
						column(6,
							style = "text-align: left",
							div(
								style = "display: inline-block;margin-top: 20px;",
								actionBttn(inputId="map_reset", label="Reset Zoom", style="pill", size="xs", color="success")
							)
						), 
						column(6,
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 12px;font-weight: 400;text-align: center; width:120px;",
								"Color By Amount of:",
								selectInput(
									"map_color",
									label = NULL,
									choices = c("RSV"),
									selected = "RSV"
								)
							),
							div(
								class = "map_embed",
								style = "margin-left: 5px; display: inline-block;font-size: 12px;font-weight: 400;text-align: center; width:120px;",
								"Plot Most Recent:",
								selectInput(
									"view_range",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (controls)
					fluidRow(
						column(12,
							leafletOutput("map_rsv", width = "100%", height = "420px")
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							fluidRow(
								column(12,
									div(
										"It is now ", printy_dates(today), ".", 
										style="font-size: 18px;padding-top: 10px;font-weight: 400;text-align: center;color: #333333;"
									)
								)
							),
							fluidRow(
								column(12,
									div(
									"Latest RSV abundance data for this region is from ",
										style="font-size: 18px;padding-top: 10px;font-weight: 800;text-align: center;color: #000000;"
									),
									div(
										textOutput("selection_freshness_rsv"), 
										style="font-size: 24px;padding-bottom: 0px;font-weight: 800;text-align: center;color: #9437FF;"
									),
									div(
										textOutput("selection_completeness_rsv"), 
										style="font-size: 16px;padding-bottom: 10px;font-weight: 400;text-align: center;color: #000000;"
									)
								)
							), # fluidRow (abundance data freshness)
							fluidRow(
								column(12,
									div(
									"NOTE: Recent data is subject to change as more sites report results.",
										style="font-size: 14px;font-style: italic;padding-top: 10px; padding-bottom: 5px;font-weight: 400;text-align: center;color: #333333"
									)
								)
							) # fluidRow (about the data)
						)	# column
					), # fluidRow (data info)
				),
				column(5,
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						#style = "margin-left: 0px; margin-right: 0px;", 
						column(12,
							div(
								textOutput("plot_title_rsv"), 
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (RSV plot title)
					fluidRow(
						style = "margin-top: 5px;", 
						column(6,
							div(
								id = "abundance_title_rsv",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_rsv",
								textOutput("abundance_rsv"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(6,
							div(
								id = "trend_title_rsv",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_rsv",
								textOutput("trend_rsv"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						)
					), # fluidRow (RSV status blocks)
					fluidRow(
						column(12,
							# Plot of RSV change over time
							plotlyOutput("plot_rsv", height="550px", width="100%")
						)
					) # fluidRow (RSV plot)
				),
				column(2,
					style = "margin-top: 8px;",
					fluidRow(
						column(12,
						style = "padding-top: 4px;padding-bottom: 4px;margin-left: 0px; margin-right: 0px;background-color: #000000;color: #ffffff;",
							div(
								class = "logo",
								id = "abundance_level_key_rsv", 
								#class = "panel panel-default",
								style = "margin 5px;padding: 3px;text-align: center;",
								div("Abundance Level Key", style="font-size: 14px;font-weight: 800;text-align: center;padding: 2px;"),
								div("VERY HIGH", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #D53E4F;"),
								div("HIGH", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #FDAE61;"),
								div("MODERATE", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #E6F598;"),
								div("LOW", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #3288BD;"),
								div("UNKNOWN", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #EEEEEE;")
							)
						) # column
					), # fluidRow (abundance level color key)
					fluidRow(
						column(12,
						style = "padding-top: 4px;padding-bottom: 4px;margin-left: 0px; margin-right: 0px;margin-top: 10px;margin-bottom: 20px;background-color: #000000;color: #ffffff;",
							div(
								class = "logo",
								id = "trend_key_rsv", 
								#class = "panel panel-default",
								style = "margin 5px;padding: 3px;text-align: center;",
								div("Trend Level Key", style="font-size: 14px;font-weight: 800;text-align: center;padding: 2px;"),
								div("SPIKING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #FF85FF;"),
								div("DESPIKING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #96F786;"),
								div("INCREASING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #D53E4F;"),
								div("VARIABLE", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #FDAE61;"),
								div("STABLE", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #E6F598;"),
								div("DECREASING", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #3288BD;"),
								div("UNKNOWN", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;color: #000000;background-color: #EEEEEE;")
							)
						) #  column
					), # fluidRow (trend color key)
					fluidRow(
						style = "margin-top: 12px; background-color: #f3f3e1; color: #000000; padding: 5px;",
						column(12,
							# explanation of plots
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px 20px 2px 2px;text-align: left;",
								"Plot values are reported as the relative abundance of pathogen after correction for daily flow and population served, averaged across all selected sites."
							)
						)
					), # fluidRow (abundance plot explanation)
					fluidRow(
						style = "margin-top: 12px; background-color: #000000; color: #ffffff; padding: 5px;",
						column(12,
							# explanation of trend lines
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px;text-align: left;",
								span("The "),
								span(style="color: #00B140;", "green "),
								span("and "),
								span(style="color: #EAAA00;", "gold "),
								span("lines represent the average abundance over the most recent "),
								span(style="color: #00B140;", "3 months "),
								span("and "),
								span(style="color: #EAAA00;", "12 months "),
								span("respectively."),
							)
						)
					), # fluidRow (trend line descriptions)
					fluidRow(
						style = "margin-top: 12px; background-color: #ffffff; color: #ffffff; padding: 5px;",
						column(12,
							#style = "font-size: 14px;padding: 3px;font-weight: 400;",
							#div("WaTCH-WV is supported by CDC-sponsored grants from the West Virginia Department of Health to West Virginia University & Marshall University.")
							div(
								style = "margin-top: 20px;text-align: center;",
								downloadBttn(outputId="download_data_rsv", label="Download RSV Data", style="pill", size="s", color="royal")
							)
						)
					) # fluidRow (data download)
				)
			),

			hidden(
				absolutePanel(
					id = "abundance_level_key_flu_popup",
					class = "mdinfo",
					top = 365, left = 710, width = 580, height = 350,
					div("Abundance Level Color Key", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "alertinfo",
						id = "level_4",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
						span(paste0("VERY HIGH. The latest amount of this disease agent is greater than 150% of the 3 month average. Community transmission is estimated to be very high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_3",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
						span(paste0("HIGH. The latest amount of this disease agent is between 100% and 150% of the 3 month average. Community transmission is estimated to be high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_2",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
						span(paste0("MODERATE. The latest amount of this disease agent is between 50% and 100% of the 3 month average. Community transmission is estimated to be moderate."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_1",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
						span(paste0("LOW. The latest amount of this disease agent is less than 50% of the 3 month average. Community transmission is estimated to be low."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_5",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
						span("UNKNOWN. The data for this disease agent is either missing or too old to make an accurate determination.", style="font-size: 13px;")
					),
					div(
						style="padding-top: 15px; padding-right: 5px; float: right;",
						actionBttn(inputId="abundance_level_key_rsv_popup_close", label="Close", style="pill", size="xs", color="success")
					) # button div
				)
			), # hidden (abundance level key popup)

			hidden(
				absolutePanel(
					id = "trend_key_rsv_popup",
					class = "mdinfo",
					top = 265, left = 710, width = 580, height = 480,
					div("Trend Color Key", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "alertinfo",
						id = "level_5",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FF85FF; background-color: #FF85FF; border: 1px solid black; border-radius: 3px;"),
						span("SPIKING. The latest abundance of this disease agent is more than 500% higher than the previous value.", style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_6",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #96F786; background-color: #96F786; border: 1px solid black; border-radius: 3px;"),
						span("DESPIKING. The latest abundance of this disease agent has dropped more than 500% from than the previous value.", style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_4",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
						span(paste0("INCREASING. The abundance of this disease agent has risen more than it has declined over the most recent period.."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_3",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
						span(paste0("VARIABLE. The abundance of this disease agent has been variable over the most recent period."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_2",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
						span(paste0("STABLE. The abundance of this disease agent has been relatively stable over the most recent period."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_1",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
						span(paste0("DECREASING. The abundance of this disease agent has generally declined over the most recent period (usually 4 weeks)."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_7",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
						span("INDETERMINATE. There is not enough data for this disease agent to make an accurate determination.", style="font-size: 13px;")
					),
					div(
						style="padding-top: 15px; padding-right: 5px; float: right;",
						actionBttn(inputId="trend_key_rsv_popup_close", label="Close", style="pill", size="xs", color="success")
					) # button div
				)
			), # hidden (trend key popup)

			hidden(
				absolutePanel(
					id = "missing_data_p4_popup",
					class = "mdinfo",
					top = 320, left = 710, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			), # hidden (missing data p4 popup)

		) # tabPanel (RSV)

	)) # navbarPage and container div
))

