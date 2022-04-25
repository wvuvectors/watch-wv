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
			"Treatment Facilities",
			div(
				class="outer",

				leafletOutput("map_wwtp", width="100%", height="100%"),

				absolutePanel(
					id = "details", 
					class = "panel panel-default",
					top = 162, left = 53, width = 600, fixed=TRUE,
					draggable = FALSE, height = 420,

					plotlyOutput("plot_wwtp", height="380px", width="100%"),
					#span(tags$i(h6("All reported values are 3-day rolling means, and adjusted for average daily flow rate.")), style="color:#045a8d; text-align: right; padding-right: 20px;"),

#					div(
#						style="background-color:#ffffff; margins: 1px 0px 5px 0px; padding: 3px; float: right;",
#						actionBttn(inputId="embiggen_open_wwtp", label="Enlarge This Plot", style="pill", size="xs", color="success")
#					)
				), # absolutePanel

				absolutePanel(
					id = "controls_wwtp",
					class = "control_panel",
					top = 583, left = 53, height = "auto", width = 600,
					fixed=TRUE, draggable = FALSE,
					fluidRow(
						column(
							width = 4,
							#style = "border-right: 1px solid #666666;",
							div(
								class = "control_group",
								selectInput(
									"infections_wwtp",
									label = NULL,
									choices = INFECTIONS, 
									selected = INFECTIONS_DEFAULT
								)
							),
							div(
								style = "margin-left: 20px;font-size: 12px;",
								prettyCheckboxGroup(
									inputId = "targets_wwtp",
									label = NULL,
									choiceNames = TARGETS,
									choiceValues = TARGET_VALUES,
									icon = icon("check-square"), 
									status = "success",
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
									inputId = "dates_wwtp",
									force_edges = TRUE,
									#width = "90%",
									label=NULL,
									choices = sort(unique(as.Date(ymd(df_watch$week_starting)))),
									selected = c(min(as.Date(ymd(df_watch$week_starting))), max(as.Date(ymd(df_watch$week_starting)))),
				#							animate=animationOptions(interval = 3000, loop = FALSE),
									grid = FALSE
								)
							)
						),
						column(
							width = 4,
							div(
								class = "control_group",
								#style = "border-left: 1px solid #666666;",
								div(
									style="font-size: 13px;font-weight: 800;margin-bottom: 20px;text-align: center;",
									span("Choose a rolling mean (days)")
								),
								sliderInput(
									inputId = "roll_wwtp",
									label=NULL,
									min=2, max=10, value=SMOOTHER_DEFAULT
								)
							)
						)
					) # fluidRow
				),					

				absolutePanel(
					class="mdblock",
					id = "site_status_panel", 
					top = 145, left = 590, height=120, width = 120,
					fixed=TRUE, draggable=FALSE,
					span("Signal at this site is", style="font-size: 10px; line-height: 12px;"),
					span(textOutput("site_signal"), style="font-size: 16px;font-weight: 800; line-height: 22px;"),
					span("over the past week", style="font-size: 12px; line-height: 18px;"),
				),
				
				hidden(
					absolutePanel(
						id = "site_status_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="site_status_info_close", label="Close", style="pill", size="xs", color="warning")
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
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="alert_level_info_close", label="Close", style="pill", size="xs", color="warning")
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
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="scope_info_close", label="Close", style="pill", size="xs", color="warning")
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
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="population_info_close", label="Close", style="pill", size="xs", color="warning")
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
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="networkpop_info_close", label="Close", style="pill", size="xs", color="warning")
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
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="daily_flow_info_close", label="Close", style="pill", size="xs", color="warning")
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
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="collection_info_close", label="Close", style="pill", size="xs", color="warning")
						)
					)
				),

				absolutePanel(
					class="mdblock",
					id = "last_date_panel", 
					top = 65, left = 930, height=75, width = 140,
					fixed=TRUE, draggable=FALSE,
					span("Last update from this site", style="font-size: 11px;"),
					span(textOutput("last_update"), style="font-size: 20px;")
				),
				
				hidden(
					absolutePanel(
						id = "last_date_info",
						class = "mdinfo",
						top = 265, left = 590, height="auto", width = 300,
						fixed=TRUE, draggable=TRUE,
						div("Hi there!"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="last_date_info_close", label="Close", style="pill", size="xs", color="warning")
						)
					)
				),

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

#			hidden(
#				div(
#					id = "conditionalPanelWWTP",
#					absolutePanel(
#						id = "controls_big_wwtp",
#						class = "panel panel-default",
#						style="background-color:#ffffff; padding: 4px; margins: 0;border: 2px solid #909090; border-radius: 5px;",
#						top = 67, left = 50, width = "90%", height = "85%", fixed=TRUE,
##						draggable = TRUE,
#						fluidRow(
#							column(
#								width = 3,
#								div(
#									style="padding-top: 15px; padding-left: 5px; float: left;",
#									materialSwitch(inputId = "plot_ci_wwtp", label = "90% CI:", value = FALSE, status="warning")
#								)
#							),
#							column(
#								width=6,
#								h4(textOutput("facility_name_big"), align = "center", style="background-color: #ffffff;color: #000000;padding: 5px; font-weight: bold;"),
#							),
#							column(
#								width=3,
#								div(
#									style="padding-top: 15px; padding-right: 5px; float: right;",
#									actionBttn(inputId="embiggen_close_wwtp", label="Close", style="pill", size="xs", color="warning")
#								)
#							)
#						),
#						span(tags$i(h6("All reported values are adjusted for average daily flow rate.")), style="color:#045a8d; text-align: left"),
#						plotlyOutput("plot_wwtp_big", height="auto", width="auto")
#					) # absolutePanel
#				) # hidden panel div
#			) # hidden panel
		), # tabPanel

		tabPanel(
			"Sewer Network Sites",
			div(
				class="outer",
				leafletOutput("map_upstream", width="100%", height="100%"),
	
				absolutePanel(
					id = "controls", 
					class = "panel panel-default",
					top = 75, left = 55, width = 525, fixed=TRUE,
					draggable = TRUE, height = "auto",

					span(tags$i(h4(textOutput("last_update_upstream"))), style="color:#000000"),
					span(tags$i(h6("Sewer network sites include any collections between source and treatment facility, including buildings, complexes, and pump stations.")), style="color:#045a8d; text-align: left"),
					h4(textOutput("all_facilities_upstream"), align = "left"),
					h5(textOutput("all_samples_upstream"), align = "left"),
					h4(textOutput("facility_name_upstream"), align = "center", style="background-color: #FFEFCE;color: #000000;padding: 5px;margin-top: 20px"),

					span(tags$i(h6("All reported values are 3-day rolling means. They are not adjusted for flow or population due to a lack of data.")), style="color:#045a8d; text-align: center"),

					plotlyOutput("plot_upstream", height="250px", width="100%"),

					div(
						style="background-color:#f0f0e3; margins: 1px 0px 10px 0px; padding: 8px;",
						span(tags$i(h5("Upstream Location Information")), style="color:#000000; text-align: center;"),
						h6(textOutput("facility_samples_upstream"), align = "left")
#						h6(textOutput("facility_capacity"), align = "left"),
#						h6(textOutput("facility_popserved"), align = "left"),
#						h6(textOutput("facility_counties"), align = "left")
					),
		
					div(
						style="background-color:#ffffff; margins: 1px 0px 5px 0px; padding: 3px; float: right;",
						actionBttn(inputId="embiggen_open_upstream", label="Enlarge This Plot", style="pill", size="xs", color="success")
					)					
				), # absolutePanel
	
				absolutePanel(
					id = "targets_popup_upstream_min",
					style = "text-align:center; background-color: #000000; color: #ffffff; border: 1px solid #dddddd; border-radius: 5px; opacity: 0.9;",
					bottom = 0, left = 5, width = 80, height = 25,
					fixed=TRUE, draggable = FALSE,
					span(tags$i("Targets"))
				),					
				absolutePanel(
					id = "targets_popup_upstream", 
					style="background-color: #000000;color: #ffffff;opacity: 1.0;padding: 3px; border: 1px solid #dddddd; border-radius: 5px; font-size: 11px;", 
					bottom = 25, left = 25, width = 300, height = 130,
					fixed=TRUE, draggable=FALSE, 
					h5("Choose the targets to display:", align = "center", style="background-color: #000000;color: #ffffff;padding: 2px;margin-top: 5px"),
					prettyCheckboxGroup(
						inputId = "targets_upstream",
						label=NULL,
						choiceNames = TARGETS,
						choiceValues = TARGET_VALUES,
						icon = icon("check-square"), 
						status = "success",
						selected = TARGETS_DEFAULT,
						outline = TRUE
	#							animation = "pulse"
					)
				),
				absolutePanel(
					id = "dates_popup_upstream_min",
					style = "text-align:center; background-color: #000000; color: #ffffff; border: 1px solid #dddddd; border-radius: 5px; opacity: 0.9;",
					bottom = 0, left = 90, width = 100, height = 25,
					fixed=TRUE, draggable = FALSE,
					span(tags$i("Date Range"))
				),					
				absolutePanel(
					id = "dates_popup_upstream", 
					style="background-color: #000000;opacity: 0.8;padding: 3px; border: 1px solid #dddddd; border-radius: 5px; font-size: 11px;", 
					bottom = 25, left = 110, width = 300, height = 130,
					fixed=TRUE, draggable=FALSE, 
					h5("Choose a range of dates to display:", align = "center", style="background-color: #000000;color: #ffffff;padding: 2px;margin-top: 5px"),
					sliderTextInput(
						inputId = "dates_upstream",
						force_edges = TRUE,
						#width = "90%",
						label=NULL,
						choices = sort(unique(as.Date(ymd(df_watch$week_starting)))),
						selected = c(min(as.Date(ymd(df_watch$week_starting))), max(as.Date(ymd(df_watch$week_starting)))),
	#							animate=animationOptions(interval = 3000, loop = FALSE),
						grid = FALSE
					)
				),
				absolutePanel(
					id = "roll_popup_upstream_min",
					style = "text-align:center; background-color: #000000; color: #ffffff; border: 1px solid #dddddd; border-radius: 5px; opacity: 0.9;",
					bottom = 0, left = 195, width = 115, height = 25,
					fixed=TRUE, draggable = FALSE,
					span(tags$i("Data Smoothing"))
				),					
				absolutePanel(
					id = "roll_popup_upstream", 
					style="background-color: #000000;color: #ffffff;opacity: 1.0;padding: 3px; border: 1px solid #dddddd; border-radius: 5px; font-size: 11px;", 
					bottom = 25, left = 215, width = 300, height = 130,
					fixed=TRUE, draggable=FALSE, 
					h5("Choose a rolling mean window (days):", align = "center", style="background-color: #000000;color: #ffffff;padding: 2px;margin-top: 5px"),
					sliderInput(
						inputId = "roll_upstream",
						label=NULL,
						min=2, max=10, value=SMOOTHER_DEFAULT
					)
				), # absolutePanel

				absolutePanel(
					id = "logo", 
					class = "card", 
					bottom = 10, right = 10, width = 80, 
					fixed=TRUE, draggable = FALSE, 
					height = "auto",
					tags$a(href='https://www.wvuvectors.com/', 
					tags$img(src='WaTCH-WV_logo.png',height='80',width='80'))
				)

#				absolutePanel(
#					id = "logo", 
#					class = "card", 
#					bottom = 10, left = 70, width = 50, 
#					fixed=TRUE, draggable = FALSE, 
#					height = "auto",
#					actionButton(
#						"twitter_share", 
#						label = "", 
#						icon = icon("twitter"), 
#						style='padding:5px',
#						onclick = sprintf("window.open('%s')","https://twitter.com/intent/tweet?text=%20@WVUvectors%20wastewater%20testing&url=https://wvuvectors.shinyapps.io/WaTCH-WV&hashtags=COVID,WaTCH-WV")
#					)
#				)
			), # tabPanel div
	
			hidden(
				div(
					id = "conditionalPanelUpstream",
					absolutePanel(
						id = "controls_big_upstream",
						class = "panel panel-default",
						style="background-color:#ffffff; padding: 4px; margins: 0;border: 2px solid #909090; border-radius: 5px;",
						top = 67, left = 50, width = "90%", height = "85%", fixed=TRUE,
#						draggable = TRUE,
						fluidRow(
							column(
								width = 3,
								div(
									style="padding-top: 15px; padding-left: 5px; float: left;",
									materialSwitch(inputId = "plot_ci_upstream", label = "90% CI:", value = FALSE, status="warning")
								)
							),
							column(
								width=6,
								h4(textOutput("facility_name_big_upstream"), align = "center", style="background-color: #ffffff;color: #000000;padding: 5px; font-weight: bold;"),
							),
							column(
								width=3,
								div(
									style="padding-top: 15px; padding-right: 5px; float: right;",
									actionBttn(inputId="embiggen_close_upstream", label="Close", style="pill", size="xs", color="warning")
								)
							)
						),
						span(tags$i(h6("All values are plotted as rolling means of RNA copies/L wastewater.")), style="color:#045a8d; text-align: left"),

						plotlyOutput("plot_upstream_big", height="auto", width="auto")
					) # hidden absolutePanel
				) # hidden panel div
			) # hidden panel

		), # tabPanel

		tabPanel(
			"Genome Sequencing",
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

