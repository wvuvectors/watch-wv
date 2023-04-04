shinyUI(bootstrapPage(
	useShinyjs(),
	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main_styles.css")
  ),
	
	navbarPage(
		theme = shinytheme("flatly"), 
		id="nav",
		HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Wastewater Testing for Community Health in WV (WaTCH)</a>'), 
		windowTitle = "WaTCH-WV",
		collapsible = FALSE,

		tabPanel(
			"SARS-CoV-2",
			div(
				class="outer",
				
				leafletOutput("watch_map", width=567, height=425),

				absolutePanel(
					id = "geolevels",
					top = 65, left = 462, height = 47, width = 100,
					fixed=TRUE, draggable=FALSE,
					div(
						class = "map_embed",
						selectInput(
							"geo_level",
							label = NULL,
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						)
					)
				), # absolutePanel
				
				absolutePanel(
					id = "location-info",
					class = "mdblock",
					top = 60, left = 570, height = 120, width = 334,
					fixed=TRUE, draggable = FALSE,
					span("Mean daily flow", style="font-size: 13px;"),
					#span(textOutput("mean_flow"), style="font-size: 16px;"),
					span("(million gallons per day)", style="font-size: 11px;")
				), # absolutePanel
				
				absolutePanel(
					id = "hospital-plot",
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 180, left = 570, width = 334, height = 295, 
					div(textOutput("plothosp_title"), style="padding-top: 2px; font-size: 14px; font-style: italic; color:#045a8d; text-align: center;"),
					plotlyOutput("plothosp", height="270px", width="100%")
#					div(style="font-size: 9px;", dataTableOutput("tableWW"))
#					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
				), # absolutePanel

				absolutePanel(
					id = "other-data",
					class = "panel panel-default",
					top = 60, left = 910, height = 885, width = 400,
					fixed=TRUE, draggable = FALSE,
					div(
						style="font-size: 13px;font-weight: 800;margin-bottom: 6px;text-align: center;background-color: #f3f3e1;",
						span("Other Related Data")
					),
					fluidRow(
						column(
							width = 12,
							div(
								#class = "control_group",
								div(textOutput("plotclass_title"), style="padding-top: 2px; font-size: 14px; font-style: italic; color:#045a8d; text-align: center;"),
								plotlyOutput("plotclass", height="300px", width="100%")
								
							)
						),
					) # fluidRow
				), # absolutePanel
				
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 485, left = 4, width = 900, height = 400, 
					div(textOutput("plotww_title"), style="padding-top: 5px; font-size: 16px; font-style: bold; color:#000000; text-align: center;"),
					fluidRow(
						column(
							width = 12,
							div(
								div(textOutput("plotww_wide_title"), style="padding-top: 2px; font-size: 14px; font-style: italic; color:#045a8d; text-align: center;"),
								plotlyOutput("plotww_wide", height="300px", width="100%")
			#					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
							) # div
#						), # column
#						column(
#							width = 6,
#							div(
#								div(textOutput("plotww_narrow_title"), style="padding-top: 2px; font-size: 14px; font-style: italic; color:#045a8d; text-align: center;"),
#								plotlyOutput("plotww_narrow", height="300px", width="100%")
#			#					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
#							) # div
						) # column
					) # fluidRow
				) # absolutePanel
				
			) # div
		), # tabPanel
		
		tabPanel(
			"Flu",
			div(
				class="outer",
				leafletOutput("map_flu", width=567, height=425)
			) # div

		), # tabPanel

		tabPanel(
			"RSV",
			div(
				class="outer",
				leafletOutput("map_rsv", width=567, height=425)
			) # div

		), # tabPanel
		
	)	# navbarPage
))

