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
			"Routine Surveillance",
			div(
				class="outer",

				leafletOutput("map_rs", width = 645, height = 527),

				absolutePanel(
					class = "logo", 
#					top = 74, left = 60, width = 80, 
					top = 445, left = 535, width = 80, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					tags$a(href='https://www.watch-wv.com/', 
					tags$img(src='WaTCH-WV_logo.png',height='80',width='80'))
				), # absolutePanel
				
				absolutePanel(
					class = "map-recenter", 
					top = 530, left = 530, width = 125, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					actionBttn(inputId="center_map_rs", label="Recenter Map", style="pill", size="xs", color="success")
				), # absolutePanel
				

				absolutePanel(
					class = "controls",
					top = 70, left = 255, height = 230, width = 450,
					fixed=TRUE, draggable=FALSE,
					div(
						class = "map_embed",
						selectInput(
							"geo_level_rs",
							label = NULL,
							#width = "200px",
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						),
						selectInput(
							"view_range_rs",
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
					top = 590, left = 5, width = 637, height = 290, 
					div(textOutput("plotsq_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plotsq_rs", height="270px", width="100%")
				), # absolutePanel

				
				# Top left (info panels)
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 648, width = 505, height = 37, 
					div(textOutput("site_rs_title"), style="font-size: 18px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff;background-color: #000000;")
				), # absolutePanel
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 104, left = 648, width = 505, height = 32, 
					div(id="rs_CSL_label", "Alert Level: ", style="display: inline-block;font-size: 13px;font-weight: 800;text-align: right;padding: 2px;width:85px;"),
					div(id="rs_CSL_1", textOutput("rs_CSL_1"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin-right: 1px;margin-left: 3px;padding: 2px;width:98px;background-color: #ee4444;"),
					div(id="rs_CSL_2", textOutput("rs_CSL_2"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:98px;background-color: #44ee44;"),
					div(id="rs_CSL_3", textOutput("rs_CSL_3"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:98px;background-color: #4444ee;"),
					div(id="rs_CSL_4", textOutput("rs_CSL_4"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:98px;background-color: #eeee44;")
				), # absolutePanel
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 138, left = 648, width = 505, height = 30, 
					div(id="rs_fresh_label", "Data Age (d):", style="display: inline-block;font-size: 13px;font-weight: 800;text-align: right;padding: 2px;width:85px;"),
					div(id="rs_fresh_1", textOutput("rs_fresh_1"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin-right: 1px;margin-left: 3px;padding: 2px;width:98px;background-color: #fafafa;"),
					div(id="rs_fresh_2", textOutput("rs_fresh_2"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:98px;background-color: #fafafa;"),
					div(id="rs_fresh_3", textOutput("rs_fresh_3"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:98px;background-color: #fafafa;"),
					div(id="rs_fresh_4", textOutput("rs_fresh_4"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:98px;background-color: #fafafa;")
				), # absolutePanel
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 170, left = 648, width = 505, height = 183, 
					span(textOutput("site_rs_info"), style="font-size: 14px; font-weight: 400;text-align: center;padding-top: 1px;"),
					span(textOutput("site_rs_flow"), style="font-size: 14px; font-weight: 400;text-align: center;padding: 1px;")
				), # absolutePanel

				# Top right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 1154, width = 345, height = 290, 
					"Table in here"
				), # absolutePanel


				# Middle left
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 358, left = 648, width = 425, height = 225, 
					div(textOutput("plot1_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot1_rs", height="210px", width="100%")
				), # absolutePanel

				# Middle right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 358, left = 1075, width = 425, height = 225, 
					div(textOutput("plot2_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot2_rs", height="210px", width="100%")
				), # absolutePanel

				# Bottom left
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 590, left = 648, width = 425, height = 225, 
					div(textOutput("plot3_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot3_rs", height="210px", width="100%")
				), # absolutePanel

				# Bottom right
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 590, left = 1075, width = 425, height = 225, 
					div(textOutput("plot4_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
					plotlyOutput("plot4_rs", height="210px", width="100%")
				), # absolutePanel

				# Tag line (below the plot grid)
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 815, left = 648, width = 850, height = 66, 
					"Info about testing in here"
				), # absolutePanel
			) # tabPanel div
		), # tabPanel

		tabPanel(
			"Mass Gatherings",
			div(
				class="outer",
				leafletOutput("map_mg", width=645, height=457),

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
					actionBttn(inputId="center_map_mg", label="Recenter Map", style="pill", size="xs", color="success") #002855
				), # absolutePanel
				
				absolutePanel(
					#class = "panel panel-default",
					class = "controls",
					top = 70, left = 260, height = 160, width = 450,
					fixed=TRUE, draggable=FALSE,
					div(
						class = "map_embed",
						selectInput(
							"geo_level_mg",
							label = NULL,
							#width = "200px",
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						),
						selectInput(
							"view_range_mg",
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
					"SEQR PLOTS",
#					div(textOutput("plot1_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot1", height="290px", width="100%")
				), # absolutePanel
				
				# Top left
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 648, width = 425, height = 225, 
					"SITE SPECIFIC INFO",
#					div(textOutput("plot2_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot2_rs", height="210px", width="100%")
				), # absolutePanel

				# Top right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 1075, width = 425, height = 225, 
					"PLOT TARGET 1",
#					div(textOutput("plot1_mg_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot1_mg", height="210px", width="100%")
				), # absolutePanel

				# Middle left
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 290, left = 648, width = 425, height = 225, 
					"PLOT TARGET 2",
#					div(textOutput("plot2_mg_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot2_mg", height="210px", width="100%")
				), # absolutePanel

				# Middle right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 290, left = 1075, width = 425, height = 225, 
					"PLOT TARGET 3",
#					div(textOutput("plot3_mg_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot3_mg", height="210px", width="100%")
				), # absolutePanel

				# Bottom left
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 520, left = 648, width = 425, height = 300, 
					"PLOT TARGET 4",
#					div(textOutput("plot4_mg_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot4_mg", height="290px", width="100%")
				), # absolutePanel

				# Bottom right
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 520, left = 1075, width = 425, height = 300, 
					"PLOT TARGET 5",
#					div(textOutput("plot5_mg_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;background-color: #f3f3e1;"),
#					plotlyOutput("plot5_mg", height="290px", width="100%")
				), # absolutePanel

			) # tabPanel div
		), # tabPanel

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
		) # tabPanel

	) # navbarPage
))

