#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

shinyUI(bootstrapPage(
	useShinyjs(),
	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main_styles.css")
  ),
	
	navbarPage(
		theme = shinytheme("flatly"), 
		id="nav",
		HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">WaTCH-WV</a>'), 
		windowTitle = "WaTCH-WV",
		collapsible = TRUE,

		tabPanel(
			"Wastewater Nowcast",
			div(
				class="outer",

				leafletOutput("watch_map", width="100%", height="100%"),

				absolutePanel(
					id = "details", 
					class = "panel panel-default",
					top = 162, left = 53, width = 600, fixed=TRUE,
					draggable = FALSE, height = 430,
					plotlyOutput("watch_plot", height="380px", width="100%"),
					span(textOutput("plot_title"), style="font-size: 14px; font-style: bold; color:#045a8d; text-align: center;"),
					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
				), # absolutePanel

				absolutePanel(
					id = "controls",
					class = "control_panel",
					top = 593, left = 53, height = "auto", width = 600,
					fixed=TRUE, draggable = FALSE,
					fluidRow(
						column(
							width = 4,
							#style = "border-right: 1px solid #666666;",
							div(
								class = "control_group",
								selectInput(
									"infections_control",
									label = NULL,
									choices = INFECTIONS, 
									selected = INFECTIONS_DEFAULT
								)
							),
							div(
								style = "margin-left: 20px;font-size: 12px;",
								prettyCheckboxGroup(
									inputId = "targets_control",
									label = NULL,
									choiceNames = TARGETS,
									choiceValues = TARGET_VALUES,
									icon = icon("check-square"), 
									status = "primary",
									selected = TARGETS_DEFAULT,
									outline = TRUE
								)
							)
						),
						column(
							width = 4,
							div(
								class = "control_group",
								div(
									style="font-size: 13px;font-weight: 800;margin-bottom: 22px;text-align: center;",
									span("Choose a range of dates")
								),
								sliderTextInput(
									inputId = "dates_control",
									force_edges = TRUE,
									#width = "90%",
									label=NULL,
									choices = sort(unique(df_watch$week_starting)),
									selected = c(min(df_watch$week_starting), max(df_watch$week_starting)),
				#							animate=animationOptions(interval = 3000, loop = FALSE),
									grid = FALSE
								)
							)
						),
						column(
							width = 4,
							div(
								class = "control_group",
								div(
									style="font-size: 13px;font-weight: 800;margin-bottom: 20px;text-align: center;",
									span("Choose a rolling mean (days)")
								),
								selectInput(
									"roll_control",
									label = NULL,
									choices = SMOOTHER_OPTIONS, 
									selected = SMOOTHER_DEFAULT
								),
								materialSwitch(inputId = "ci_control", label = "Show 90% CI:", value = FALSE, status="primary")
							)
						)
					) # fluidRow
				),					

				absolutePanel(
					class="mdblock",
					id = "site_status_panel", 
					top = 145, left = 590, height=120, width = 130,
					fixed=TRUE, draggable=FALSE,
					span("Signal at this site is", style="font-size: 10px; line-height: 12px;"),
					span(textOutput("site_signal"), style="font-size: 18px;font-weight: 800; line-height: 22px;"),
					span("over the past 2 weeks", style="font-size: 12px; line-height: 18px;"),
				),
				
				hidden(
					absolutePanel(
						id = "site_status_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("What do the various levels of behavior mean, and how are they calculated? Maybe show a plot focus on last several weeks?"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="site_status_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "alert_panel", 
					top = 65, left = 45, height=75, width = 125,
					fixed=TRUE, draggable=FALSE,
					span("Current alert level is", style="font-size: 13px;"),
					span(textOutput("alert_level"), style="font-size: 24px; font-weight: 800; color: #009900;line-height: 28px;")
				),
				
				hidden(
					absolutePanel(
						id = "alert_level_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("What does the alert level mean?"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="alert_level_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "scope_panel", 
					top = 65, left = 170, height=75, width = 240,
					fixed=TRUE, draggable=FALSE,
					span(textOutput("scope"), style="font-size: 19px; font-weight: 800; line-height: 26px;"),
					span(textOutput("scope_count"), style="font-size: 15px; font-weight: 200;"),
					span(textOutput("sample_count"), style="font-size: 11px; font-weight: 200;")
				),
				
				hidden(
					absolutePanel(
						id = "scope_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Plot showing receipt of samples?"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="scope_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "population_panel", 
					top = 65, left = 410, height=75, width = 130,
					fixed=TRUE, draggable=FALSE,
					style = "padding-top: 8px;",
					span(textOutput("counties_served"), style="font-size: 14px;"),
					span("total population", style="font-size: 10px;"),
					span(textOutput("county_population"), style="font-size: 15px;")
				),
				
				hidden(
					absolutePanel(
						id = "population_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("If more than one county, show some metadata about all (name, population, etc.)."),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="population_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "networkpop_panel", 
					top = 65, left = 540, height=75, width = 130,
					fixed=TRUE, draggable=FALSE,
					span("Population on network", style="font-size: 12px; line-height: 12px;"),
					span(textOutput("population_served"), style="font-size: 18px;line-height: 22px;"),
					#span("or", style="font-size: 10px; line-height: 14px;"),
					span(textOutput("population_served_pct"), style="font-size: 12px;line-height: 18px;")
				),
				
				hidden(
					absolutePanel(
						id = "networkpop_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Anything to show/tell here?"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="networkpop_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "daily_flow_panel", 
					top = 65, left = 670, height=75, width = 130,
					fixed=TRUE, draggable=FALSE,
					span("Mean daily flow", style="font-size: 13px;"),
					span(textOutput("mean_flow"), style="font-size: 16px;"),
					span("(million gallons per day)", style="font-size: 11px;")
				),
				
				hidden(
					absolutePanel(
						id = "daily_flow_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Plot showing daily flow over time at this site, with capacity."),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="daily_flow_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "collection_panel", 
					top = 65, left = 800, height=75, width = 130,
					fixed=TRUE, draggable=FALSE,
					span("This facility is providing", style="font-size: 10px;"),
					span(textOutput("collection_frequency"), style="font-size: 13px;line-height: 12px;"),
					span("over the last 4 weeks", style="font-size: 10px; line-height: 11px;")
					#span("as", style="font-size: 10px; line-height: 11px;"),
					#span((textOutput("collection_scheme")), style="font-size: 11px;line-height: 11px;")
				),
				
				hidden(
					absolutePanel(
						id = "collection_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Plot showing samples submitted (maybe number/week) over time for this site."),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="collection_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),

				absolutePanel(
					class="mdblock",
					id = "last_date_panel", 
					top = 65, left = 930, height=75, width = 140,
					fixed=TRUE, draggable=FALSE,
					span("Last update was", style="font-size: 16px;"),
					span(textOutput("last_update"), style="font-size: 20px;")
				),
				
				hidden(
					absolutePanel(
						id = "last_date_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Text describing how samples are received by this site (e.g., 'usually on a Wed' kind of thing)."),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="last_date_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),

				absolutePanel(
					id = "recenter_panel", 
					class = "card", 
					bottom = 28, right = 125, width = 80, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					actionBttn(inputId="center_map", label="Recenter", style="pill", size="xs", color="success")
				), # absolutePanel
				
				absolutePanel(
					id = "logo", 
					class = "card", 
					bottom = 5, left = 5, width = 80, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					tags$a(href='https://www.wvuvectors.com/', 
					tags$img(src='WaTCH-WV_logo.png',height='80',width='80'))
				), # absolutePanel
				
			) # div
		), # tabPanel

		tabPanel(
			"Variant analysis",
			tags$div(
				tags$h4("Coming soon..."), 
				#h6(paste0(update)),
#				"This site is updated weekly or biweekly. ", 
#				tags$a(href="https://experience.arcgis.com/experience/685d0ace521648f8a5beeeee1b9125cd", "the WHO,"),
				"This tab will present results of our RNA virus genome sequencing efforts."
			)
		), # tabPanel

		tabPanel(
			"Ethics of testing",
			tags$div(
				tags$h4("Coming soon..."), 
				#h6(paste0(update)),
#				"This site is updated weekly or biweekly. ", 
#				tags$a(href="https://experience.arcgis.com/experience/685d0ace521648f8a5beeeee1b9125cd", "the WHO,"),
				"This tab will present results of our RNA virus genome sequencing efforts."
			)
		), # tabPanel

		tabPanel(
			"About this site",
			div(
				class="info",
				tags$i(paste0("This is WaTCH dashboard version ", WATCH_VERSION, ", released on ", format(ymd(WATCH_RELEASE_DATE), format="%d %b %Y"), ".",  sep=""))
			),
			div(
				class="info",
				"Wastewater monitoring is a valuable and cost-effective method for detecting certain infections in the community, including COVID-19, influenza, and RSV. ",
				"A robust wastewater testing program represents an early-warning system, particularly in communities with little or no individual testing, and can inform effective public health interventions. "
			),
			div(
				class="info",
				"In contrast to individual testing, wastewater surveillance does not require individuals to seek care and can passively track the change in community infection levels over time. ",
				"It is useful in detecting the beginning of outbreaks, may also be valuable in estimating disease burden in a community, and can be adapted readily to monitor future epidemics of novel diseases. ",
				"This is a key approach to directing appropriate medical care – particularly for populations that are vulnerable or hesitant about seeking care."
			),
			div(
				class="info",
				"West Virginia is one of the most vulnerable states to infections such as COVID-19 based on socioeconomic and health status, household composition/disability, and epidemiological factors. ",
				"Almost all the state is rural, 20% of WV population is 65 and older, and WV ranks at or near the bottom in most US chronic disease categories. This is especially concerning due to the high ",
				"percentage of grandparents living with school age children across WV (up to over 70% in many counties). Moreover, almost the entire state is defined as medically underserved by the Health Resources and Services Administration."
			)
	#		tags$a(href="https://experience.arcgis.com/experience/685d0ace521648f8a5beeeee1b9125cd", "the WHO,"),
		) # tabPanel
	)	# navbarPage
))

