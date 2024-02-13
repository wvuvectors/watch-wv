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
			"Respiratory",
			div(
				class="outer",

				leafletOutput("watch_map", width = 690, height = 450),

				absolutePanel(
					class = "panel panel-default",
					id = "viewranges",
					top = 505, left = 0, height = 50, width = 100,
					fixed=TRUE, draggable=FALSE,
					div(
						#class = "map_embed",
						selectInput(
							"view_range",
							label = NULL,
							choices = VIEW_RANGES, 
							selected = VIEW_RANGE_PRIMARY
						)
					)
				), # absolutePanel
				
				absolutePanel(
					class = "panel panel-default",
					id = "geolevels",
					top = 505, left = 150, height = 50, width = 100,
					fixed=TRUE, draggable=FALSE,
					div(
						#class = "map_embed",
						selectInput(
							"geo_level",
							label = NULL,
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						)
					)
				), # absolutePanel

				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 75, left = 700, width = 700, height = 280, 
					div(textOutput("plot1_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot1", height="250px", width="100%")
				), # absolutePanel

				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 360, left = 700, width = 700, height = 280, 
					div(textOutput("plot2_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot2", height="250px", width="100%")
				), # absolutePanel

				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 645, left = 700, width = 700, height = 280, 
					div(textOutput("plot3_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot3", height="250px", width="100%")
				) # absolutePanel
				
			) # tabPanel div
		), # tabPanel

		tabPanel(
			"Emerging",
			div(
				class="outer",
				leafletOutput("map_flu", width=567, height=425)
			) # tabPanel div
		) # tabPanel

	) # navbarPage
))

