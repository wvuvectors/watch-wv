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

				leafletOutput("map_resp", width = 645, height = 457),

				absolutePanel(
					class = "logo", 
					top = 74, left = 60, width = 80, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					tags$a(href='https://www.wvuvectors.com/', 
					tags$img(src='WaTCH-WV_logo.png',height='80',width='80'))
				), # absolutePanel
				
				absolutePanel(
					class = "map-recenter", 
					top = 460, left = 530, width = 125, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					actionBttn(inputId="center_map_resp", label="Recenter Map", style="pill", size="xs", color="success")
				), # absolutePanel
				

				absolutePanel(
					#class = "panel panel-default",
					class = "controls",
					top = 70, left = 260, height = 160, width = 450,
					fixed=TRUE, draggable=FALSE,
					div(
						class = "map_embed",
						selectInput(
							"geo_level_resp",
							label = NULL,
							#width = "200px",
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						),
						selectInput(
							"view_range_resp",
							label = NULL,
							#width = "60px",
							choices = VIEW_RANGES, 
							selected = VIEW_RANGE_PRIMARY
						)
					)
				), # absolutePanel

				# SARS seqr data
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 520, left = 5, width = 637, height = 300, 
					"SEQR PLOT",
#					div(textOutput("plot1_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot1", height="290px", width="100%")
				), # absolutePanel
				
				# SARS trend data
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 520, left = 647, width = 850, height = 300, 
					div(textOutput("plot1_resp_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot1_resp", height="290px", width="100%")
				), # absolutePanel


				# Top left
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 648, width = 425, height = 225, 
					"SITE SPECIFIC INFO",
#					div(textOutput("plot2_resp_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot2_resp", height="210px", width="100%")
				), # absolutePanel

				# Top right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 1075, width = 425, height = 225, 
					div(textOutput("plot2_resp_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot2_resp", height="210px", width="100%")
				), # absolutePanel


				# Bottom left
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 290, left = 648, width = 425, height = 225, 
					div(textOutput("plot4_resp_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot4_resp", height="210px", width="100%")
				), # absolutePanel

				# Bottom right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 290, left = 1075, width = 425, height = 225, 
					div(textOutput("plot3_resp_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot3_resp", height="210px", width="100%")
				) # absolutePanel
			) # tabPanel div
		), # tabPanel

		tabPanel(
			"Emerging",
			div(
				class="outer",
				leafletOutput("map_emerg", width=645, height=457),

				absolutePanel(
					class = "logo", 
					top = 74, left = 60, width = 80, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					tags$a(href='https://www.wvuvectors.com/', 
					tags$img(src='WaTCH-WV_logo.png',height='80',width='80'))
				), # absolutePanel
				
				absolutePanel(
					class = "map-recenter", 
					top = 460, left = 530, width = 125, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					actionBttn(inputId="center_map_emerg", label="Recenter Map", style="pill", size="xs", color="success")
				), # absolutePanel
				
				absolutePanel(
					#class = "panel panel-default",
					class = "controls",
					top = 70, left = 260, height = 160, width = 450,
					fixed=TRUE, draggable=FALSE,
					div(
						class = "map_embed",
						selectInput(
							"geo_level_emerg",
							label = NULL,
							#width = "200px",
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						),
						selectInput(
							"view_range_emerg",
							label = NULL,
							#width = "60px",
							choices = VIEW_RANGES, 
							selected = VIEW_RANGE_PRIMARY
						)
					)
				), # absolutePanel

			) # tabPanel div
		) # tabPanel

	) # navbarPage
))

