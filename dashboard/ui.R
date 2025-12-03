shinyUI(fluidPage(
	useShinyjs(),
	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main_styles.css")
  ),
	
	tags$div(class="nbcontainer", 
	
	navbarPage(
		theme = shinytheme("flatly"), 
		id="nav",
		HTML('<a style="text-decoration:none;cursor:default;color:#FFEBCF;" class="active" href="#">Wastewater Testing for Community Health in WV (WaTCH-WV)</a>'), 
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
					fluidRow(
						column(12,
							leafletOutput("map_covid", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 80,
								top = 320,
								style = "z-index: 1000;"
							)
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
								textOutput("aplot_title_covid"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						column(10,
							div(
								textOutput("selection_freshness_covid"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-top: 5px;font-weight: 800;text-align: center;"
							),
							div(
								textOutput("selection_completeness_covid"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-bottom: 5px;font-weight: 800;text-align: center;"
							)
						),
						column(2,
							div(
								#class = "map_embed",
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px;",
								"Plot Most Recent:",
								selectInput(
									"view_range",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_covid", height="350px", width="100%")
						)
					), # fluidRow (abundance plot)
					fluidRow(
						column(8,
							fluidRow(
								style = "margin-top: 4px;margin-left: 9px;margin-right: 3px;", 
								div(
									textOutput("vplot_title_covid"), 
									style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
								)
							), # fluidRow (variant plot title)
							fluidRow(
								column(12,
									# Plot of variant proportions over time
									plotlyOutput("vplot_covid", height="300px", width="100%")
								)
							) # fluidRow (variant plot)
						),
						column(4,
							#style = "border: 2px solid #000000;padding: 5px;",
							fluidRow(
								style = "margin-top: 4px;margin-left: 3px;margin-right: 9px;", 
								div(
									"ESSENCE Data", 
									style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
								)
							), # fluidRow (ESSENCE plot title)
							div(
								style = "margin-top: 20px;text-align: center;",
								downloadBttn(outputId="download_data_covid", label="Download COVID Data", style="pill", size="s", color="royal")
							)
						)
					) # fluidRow (variant and ESSENCE plots)
				)
			),

			hidden(
				absolutePanel(
					id = "missing_data_p1_popup",
					class = "mdinfo",
					top = 320, left = 880, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			), # hidden (missing data p1 popup)
			hidden(
				absolutePanel(
					id = "missing_data_p0_popup",
					class = "mdinfo",
					top = 660, left = 705, width = 450, height = 100,
					div("There is no variant data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			) # hidden (missing data p0 popup)

		), # tabPanel (COVID)

		tabPanel(
			"FluA",
			fluidRow(
				column(12,
					# title
					style = "border: 0px solid #000000;margin-left: 0px; margin-right: 0px;",
					div(
						textOutput("selection_title_flua"), 
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
								id = "abundance_title_flua",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_flua",
								textOutput("abundance_flua"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(4,
							div(
								id = "trend_title_flua",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_flua",
								textOutput("trend_flua"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(4,
							div(
								id = "variant_title_flua",
								"Variant",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "variant_text_flua",
								textOutput("variant_flua"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						)
					), # fluidRow (alert blocks)
					fluidRow(
						# Alert details
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFE7E5;color: #000000;",
						column(12,
							div(
								textOutput("alert_details_flua"), 
								style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert details)
					fluidRow(
						# Selection details
						style = "margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_flua"), 
								style="font-size: 16px;font-weight: 800;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_flua", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 80,
								top = 320,
								style = "z-index: 1000;"
							)
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
								textOutput("aplot_title_flua"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						column(10,
							div(
								textOutput("selection_freshness_flua"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-top: 5px;font-weight: 800;text-align: center;"
							),
							div(
								textOutput("selection_completeness_flua"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-bottom: 5px;font-weight: 800;text-align: center;"
							)
						),
						column(2,
							div(
								#class = "map_embed",
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px;",
								"Plot Most Recent:",
								selectInput(
									"view_range",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_flua", height="350px", width="100%")
						)
					), # fluidRow (abundance plot)
					fluidRow(
						column(12,
							div(
								#textOutput("eplot_title_flua"), 
								"ESSENCE Data Here",
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (ESSENCE plot title)
					fluidRow(
						column(12,
							# Plot of ESSENCE data over time
							#plotlyOutput("eplot_flua", height="350px", width="100%")
							"ESSENCE Plot Here!"
						)
					), # fluidRow (ESSENCE plot)
				)
			),

			hidden(
				absolutePanel(
					id = "missing_data_p2_popup",
					class = "mdinfo",
					top = 320, left = 880, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			) # hidden (missing data p2 popup)
		), # tabPanel (FluA)

		tabPanel(
			"FluB",
			fluidRow(
				column(12,
					# title
					style = "border: 0px solid #000000;margin-left: 0px; margin-right: 0px;",
					div(
						textOutput("selection_title_flub"), 
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
						column(2),
						column(4,
							div(
								id = "abundance_title_flub",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_flub",
								textOutput("abundance_flub"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(4,
							div(
								id = "trend_title_flub",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_flub",
								textOutput("trend_flub"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(2)
					), # fluidRow (alert blocks)
					fluidRow(
						# Alert details
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFE7E5;color: #000000;",
						column(12,
							div(
								textOutput("alert_details_flub"), 
								style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert details)
					fluidRow(
						# Selection details
						style = "margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_flub"), 
								style="font-size: 16px;font-weight: 800;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_flub", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 80,
								top = 320,
								style = "z-index: 1000;"
							)
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
								textOutput("aplot_title_flub"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						column(10,
							div(
								textOutput("selection_freshness_flub"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-top: 5px;font-weight: 800;text-align: center;"
							),
							div(
								textOutput("selection_completeness_flub"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-bottom: 5px;font-weight: 800;text-align: center;"
							)
						),
						column(2,
							div(
								#class = "map_embed",
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px;",
								"Plot Most Recent:",
								selectInput(
									"view_range",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_flub", height="350px", width="100%")
						)
					), # fluidRow (abundance plot)
					fluidRow(
						column(12,
							div(
								#textOutput("eplot_title_flub"), 
								"ESSENCE Data Here",
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (ESSENCE plot title)
					fluidRow(
						column(12,
							# Plot of ESSENCE data over time
							#plotlyOutput("eplot_flub", height="350px", width="100%")
							"ESSENCE Plot Here!"
						)
					), # fluidRow (ESSENCE plot)
				)
			),

			hidden(
				absolutePanel(
					id = "missing_data_p3_popup",
					class = "mdinfo",
					top = 320, left = 880, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			) # hidden (missing data p2 popup)
		), # tabPanel (FluB)

		tabPanel(
			"RSV",
			fluidRow(
				column(12,
					# title
					style = "border: 0px solid #000000;margin-left: 0px; margin-right: 0px;",
					div(
						textOutput("selection_title_rsv"), 
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
						column(2),
						column(4,
							div(
								id = "abundance_title_rsv",
								"Abundance",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_rsv",
								textOutput("abundance_rsv"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(4,
							div(
								id = "trend_title_rsv",
								"Trend",
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_rsv",
								textOutput("trend_rsv"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(2)
					), # fluidRow (alert blocks)
					fluidRow(
						# Alert details
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFE7E5;color: #000000;",
						column(12,
							div(
								textOutput("alert_details_rsv"), 
								style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert details)
					fluidRow(
						# Selection details
						style = "margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_rsv"), 
								style="font-size: 16px;font-weight: 800;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_rsv", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 80,
								top = 320,
								style = "z-index: 1000;"
							)
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
								textOutput("aplot_title_rsv"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						column(10,
							div(
								textOutput("selection_freshness_rsv"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-top: 5px;font-weight: 800;text-align: center;"
							),
							div(
								textOutput("selection_completeness_rsv"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 16px;padding-bottom: 5px;font-weight: 800;text-align: center;"
							)
						),
						column(2,
							div(
								#class = "map_embed",
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px;",
								"Plot Most Recent:",
								selectInput(
									"view_range",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_rsv", height="350px", width="100%")
						)
					), # fluidRow (abundance plot)
					fluidRow(
						column(12,
							div(
								#textOutput("eplot_title_rsv"), 
								"ESSENCE Data Here",
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (ESSENCE plot title)
					fluidRow(
						column(12,
							# Plot of ESSENCE data over time
							#plotlyOutput("eplot_rsv", height="350px", width="100%")
							"ESSENCE Plot Here!"
						)
					), # fluidRow (ESSENCE plot)
				)
			),

			hidden(
				absolutePanel(
					id = "missing_data_p4_popup",
					class = "mdinfo",
					top = 320, left = 880, width = 450, height = 100,
					div("There is no data for this target at the selected region during the requested time period. Please try a different region or time period.", 
							style="font-size: 16px;padding: 4px;font-weight: 800;text-align: center;color: #000000;")
				)
			) # hidden (missing data p4 popup)
		), # tabPanel (RSV)
		
		tabPanel(
			HTML("<span style=\"color: #FFC749; font-style: italic;\">Get The Data</span>")
		) # tabPanel (Get The Data)
		
	)) # navbarPage and container div
))

