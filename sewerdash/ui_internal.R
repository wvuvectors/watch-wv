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
			style = "border: 0px solid #000000; background-color: #ffffff;padding: 0px;margin: 0px;", 
			fluidRow(
				column(5, 
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-right: 5px;",
					fluidRow(
						# Alert blocks
						style = "border: 1px solid #DDC0AD;padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFEBCF;color: #000000;",
						column(4,
							div(
								id = "abundance_title_covid",
								textOutput("abundance_head_covid"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_covid",
								textOutput("abundance_covid"),
								style="margin-bottom: 8px;font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;"),
							div(
								id = "trend_title_covid",
								textOutput("trend_head_covid"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_covid",
								textOutput("trend_covid"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						# Alert details
						column(8,
							div(
								textOutput("alert_details_covid"), 
								style="font-size: 18px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert blocks and details)
					fluidRow(
						# Selection details
						style = "border: 1px solid #000000;margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_covid"), 
								style="font-size: 13px;font-weight: 400;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_covid", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 75,
								top = 290,
								style = "z-index: 1000;"
							)
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								"Click any county on the map above to see specific data for that region. ", 
								"Click anywhere else on the map to return to statewide results. Each map ",
								"county is colored according to the amount (abundance) of the disease present, ", 
								"with the county border colored according to the trend of the disease. ", 
								"IMPORTANT: Data is subject to change as additional sites report results. ",
								"The dashboard updates weekly except during major holidays.", 
								style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
							)
						)	# column
					) # fluidRow (footnotes)
				),
				column(7,
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-left: 5px;",
					fluidRow(
						column(12,
							div(
								textOutput("aplot_title_covid"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						style = "border: 0px solid #C6BFDD;", 
						column(10,
							div(
								textOutput("selection_freshness_covid"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 15px;padding: 5px;font-weight: 400;text-align: center;"
							)
						),
						column(2,
							div(
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px; height: 65px;",
								#"Plot Most Recent:",
								selectInput(
									"view_range",
									label = "Plot Most Recent:",
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_covid", height="350px", width="100%")
						)
					), # fluidRow (plot of abundance data over time)
					fluidRow(
						style = "border: 0px solid #000000;margin-top: 7px;",
						column(12,
							div(
								"ESSENCE data for COVID",
								#textOutput("vplot_title_covid"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (ESSENCE plot title)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of ESSENCE data
							plotlyOutput("eplot_covid", height="285px", width="100%")
						)
					) # fluidRow (ESSENCE plot)
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
			style = "border: 0px solid #000000; background-color: #ffffff;padding: 0px;margin: 0px;", 
			fluidRow(
				column(5, 
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-right: 5px;",
					fluidRow(
						# Alert blocks
						style = "border: 1px solid #DDC0AD;padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFEBCF;color: #000000;",
						column(4,
							div(
								id = "abundance_title_flua",
								textOutput("abundance_head_flua"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_flua",
								textOutput("abundance_flua"),
								style="margin-bottom: 8px;font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;"),
							div(
								id = "trend_title_flua",
								textOutput("trend_head_flua"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_flua",
								textOutput("trend_flua"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						# Alert details
						column(8,
							div(
								textOutput("alert_details_flua"), 
								style="font-size: 18px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert blocks and details)
					fluidRow(
						# Selection details
						style = "border: 1px solid #000000;margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_flua"), 
								style="font-size: 13px;font-weight: 400;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_flua", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 75,
								top = 290,
								style = "z-index: 1000;"
							)
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								"Click any county on the map above to see specific data for that region. ", 
								"Click anywhere else on the map to return to statewide results. Each map ",
								"county is colored according to the amount (abundance) of the disease present, ", 
								"with the county border colored according to the trend of the disease. ", 
								"IMPORTANT: Data is subject to change as additional sites report results. ",
								"The dashboard updates weekly except during major holidays.", 
								style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
							)
						)	# column
					) # fluidRow (footnotes)
				),
				column(7,
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-left: 5px;",
					fluidRow(
						column(12,
							div(
								textOutput("aplot_title_flua"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						style = "border: 0px solid #C6BFDD;", 
						column(10,
							div(
								textOutput("selection_freshness_flua"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 15px;padding: 5px;font-weight: 400;text-align: center;"
							)
						),
						column(2,
							div(
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px; height: 65px;",
								#"Plot Most Recent:",
								selectInput(
									"view_range",
									label = "Plot Most Recent:",
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_flua", height="350px", width="100%")
						)
					), # fluidRow (plot of abundance data over time)
					fluidRow(
						style = "border: 0px solid #000000;margin-top: 7px;",
						column(12,
							div(
								"ESSENCE data for Flu A",
								#textOutput("eplot_title_flua"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (ESSENCE plot title)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of ESSENCE data
							plotlyOutput("eplot_flua", height="285px", width="100%")
						)
					) # fluidRow (ESSENCE plot)
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

		), # tabPanel (FluA)

		tabPanel(
			"FluB",
			style = "border: 0px solid #000000; background-color: #ffffff;padding: 0px;margin: 0px;", 
			fluidRow(
				column(5, 
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-right: 5px;",
					fluidRow(
						# Alert blocks
						style = "border: 1px solid #DDC0AD;padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFEBCF;color: #000000;",
						column(4,
							div(
								id = "abundance_title_flub",
								textOutput("abundance_head_flub"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_flub",
								textOutput("abundance_flub"),
								style="margin-bottom: 8px;font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;"),
							div(
								id = "trend_title_flub",
								textOutput("trend_head_flub"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_flub",
								textOutput("trend_flub"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						# Alert details
						column(8,
							div(
								textOutput("alert_details_flub"), 
								style="font-size: 18px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert blocks and details)
					fluidRow(
						# Selection details
						style = "border: 1px solid #000000;margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_flub"), 
								style="font-size: 13px;font-weight: 400;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_flub", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 75,
								top = 290,
								style = "z-index: 1000;"
							)
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								"Click any county on the map above to see specific data for that region. ", 
								"Click anywhere else on the map to return to statewide results. Each map ",
								"county is colored according to the amount (abundance) of the disease present, ", 
								"with the county border colored according to the trend of the disease. ", 
								"IMPORTANT: Data is subject to change as additional sites report results. ",
								"The dashboard updates weekly except during major holidays.", 
								style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
							)
						)	# column
					) # fluidRow (footnotes)
				),
				column(7,
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-left: 5px;",
					fluidRow(
						column(12,
							div(
								textOutput("aplot_title_flub"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						style = "border: 0px solid #C6BFDD;", 
						column(10,
							div(
								textOutput("selection_freshness_flub"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 15px;padding: 5px;font-weight: 400;text-align: center;"
							)
						),
						column(2,
							div(
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px; height: 65px;",
								#"Plot Most Recent:",
								selectInput(
									"view_range",
									label = "Plot Most Recent:",
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_flub", height="350px", width="100%")
						)
					), # fluidRow (plot of abundance data over time)
					fluidRow(
						style = "border: 0px solid #000000;margin-top: 7px;",
						column(12,
							div(
								"ESSENCE data for Flu B",
								#textOutput("vplot_title_flub"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (ESSENCE plot title)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of ESSENCE data
							plotlyOutput("eplot_flub", height="285px", width="100%")
						)
					) # fluidRow (ESSENCE plot)
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

		), # tabPanel (FluB)

		tabPanel(
			"RSV",
			style = "border: 0px solid #000000; background-color: #ffffff;padding: 0px;margin: 0px;", 
			fluidRow(
				column(5, 
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-right: 5px;",
					fluidRow(
						# Alert blocks
						style = "border: 1px solid #DDC0AD;padding-top: 4px;padding-bottom: 4px;margin-bottom: 5px; margin-left: 0px; margin-right: 0px;background-color: #FFEBCF;color: #000000;",
						column(4,
							div(
								id = "abundance_title_rsv",
								textOutput("abundance_head_rsv"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "abundance_text_rsv",
								textOutput("abundance_rsv"),
								style="margin-bottom: 8px;font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;"),
							div(
								id = "trend_title_rsv",
								textOutput("trend_head_rsv"),
								style="font-size: 14px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								id = "trend_text_rsv",
								textOutput("trend_rsv"),
								style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						# Alert details
						column(8,
							div(
								textOutput("alert_details_rsv"), 
								style="font-size: 18px;font-weight: 800;text-align: center;padding-bottom: 5px;"
							)
						)	# column
					), # fluidRow (alert blocks and details)
					fluidRow(
						# Selection details
						style = "border: 1px solid #000000;margin-bottom: 5px;margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								textOutput("selection_details_rsv"), 
								style="font-size: 13px;font-weight: 400;text-align: center;padding-bottom: 3px;"
							)
						)	# column
					), # fluidRow (selection info)
					fluidRow(
						column(12,
							leafletOutput("map_rsv", width = "100%", height = "470px"),
							fixedPanel(
								actionBttn(inputId="map_reset", label="Reset Map", style="pill", size="xs", color="success"),
								left = 75,
								top = 290,
								style = "z-index: 1000;"
							)
						)
					), # fluidRow (map)
					fluidRow(
						style = "padding-top: 4px;padding-bottom: 4px;margin-bottom: 12px; margin-left: 0px; margin-right: 0px;background-color: #fbfbfb;color: #000000;",
						column(12,
							div(
								"Click any county on the map above to see specific data for that region. ", 
								"Click anywhere else on the map to return to statewide results. Each map ",
								"county is colored according to the amount (abundance) of the disease present, ", 
								"with the county border colored according to the trend of the disease. ", 
								"IMPORTANT: Data is subject to change as additional sites report results. ",
								"The dashboard updates weekly except during major holidays.", 
								style="font-size: 14px;font-style: italic;font-weight: 400;text-align: center;color: #000000"
							)
						)	# column
					) # fluidRow (footnotes)
				),
				column(7,
					style = "margin-top: 5px;border: 0px solid #000000;margin-left: 0px; margin-right: 0px;padding: 0px;padding-left: 5px;",
					fluidRow(
						column(12,
							div(
								textOutput("aplot_title_rsv"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (abundance plot title)
					fluidRow(
						style = "border: 0px solid #C6BFDD;", 
						column(10,
							div(
								textOutput("selection_freshness_rsv"), 
								style="color: #000000;background-color: #E8E1FF;font-size: 15px;padding: 5px;font-weight: 400;text-align: center;"
							)
						),
						column(2,
							div(
								style = "display: inline-block;font-size: 12px;font-weight: 800;text-align: center; width:110px; height: 65px;",
								#"Plot Most Recent:",
								selectInput(
									"view_range",
									label = "Plot Most Recent:",
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						)
					), # fluidRow (abundance data freshness)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of abundance over time
							plotlyOutput("aplot_rsv", height="350px", width="100%")
						)
					), # fluidRow (plot of abundance data over time)
					fluidRow(
						style = "border: 0px solid #000000;margin-top: 7px;",
						column(12,
							div(
								"Alert history for RSV",
								#textOutput("vplot_title_rsv"), 
								style="font-size: 16px;padding: 5px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"
							)
						)
					), # fluidRow (alert history plot title)
					fluidRow(
						style = "border: 1px solid #000000;margin-top: 0px;margin-left: 0px; margin-right: 0px;padding: 0px;",
						column(12,
							# Plot of ESSENCE data
							plotlyOutput("eplot_rsv", height="285px", width="100%")
						)
					) # fluidRow (alert history plot)
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

		), # tabPanel (RSV)
		
		tabPanel(
			HTML("<span style=\"color: #FFC749; font-style: italic;\">Get The Data</span>")
		) # tabPanel (Get The Data)
		
	)) # navbarPage and container div
))

