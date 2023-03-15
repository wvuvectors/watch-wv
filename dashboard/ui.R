shinyUI(bootstrapPage(
	useShinyjs(),
	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main_styles.css")
  ),
	
	navbarPage(
		theme = shinytheme("flatly"), 
		id="nav",
		HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Wastewater Testing for Community Health in WV</a>'), 
		windowTitle = "WaTCH-WV",
		collapsible = FALSE,

		tabPanel(
			"Viral WaTCH",
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
					id = "controls",
					#class = "panel panel-default",
					top = 60, left = 635, height = 120, width = 320,
					fixed=TRUE, draggable = FALSE,
					div(
						style="font-size: 13px;font-weight: 800;margin-bottom: 6px;text-align: center;background-color: #f3f3e1;",
						span("Pathogen Targets and Test Loci")
					),
					fluidRow(
						column(
							width = 6,
							div(
								class = "control_group",
								prettyCheckboxGroup(
									inputId = "target_control",
									label = NULL,
									choices = TARGETS,
									icon = icon("check-square"), 
									status = "primary",
									selected = TARGETS_DEFAULT,
									outline = TRUE
								)
							)
						),
						column(
							width = 6,
							div(
								class = "control_group",
								prettyCheckboxGroup(
									inputId = "locus_control",
									label = NULL,
									choices = LOCI,
									icon = icon("check-square"), 
									status = "primary",
									selected = LOCI_DEFAULT,
									outline = TRUE
								)
							)
						)
					) # fluidRow
				), # absolutePanel
				
				absolutePanel(
					id = "tableview", 
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 170, left = 570, width = 290, height = 310, 
#					div(textOutput("plot_title"), style="padding-top: 10px; font-size: 17px; font-style: bold; color:#045a8d; text-align: center;"),
					div(style="font-size: 9px;", dataTableOutput("tableWW"))
#					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
				), # absolutePanel

				absolutePanel(
					id = "plotview", 
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 540, left = 0, width = 900, height = 350, 
#					div(textOutput("plot_title"), style="padding-top: 10px; font-size: 17px; font-style: bold; color:#045a8d; text-align: center;"),
					plotlyOutput("plotWW", height="300px", width="100%")
#					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
				) # absolutePanel
				
			) # div
		) # tabPanel

	)	# navbarPage
))

