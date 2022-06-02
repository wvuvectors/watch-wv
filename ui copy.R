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
					id = "metadata", 
					class = "panel panel-default",
					bottom = 0, left = 0, height=100, width = "100%",
					fixed=TRUE, draggable = FALSE,
					div(h4("METADATA HERE"))
				),
				
				absolutePanel(
					id = "details", 
					class = "panel panel-default",
					top = 75, right = 10, width = 600, fixed=TRUE,
					draggable = TRUE, height = "auto",

					span(tags$i(h4(textOutput("last_update"))), style="color:#000000"),
					h4(textOutput("all_facilities"), align = "left"),
					h5(textOutput("all_samples"), align = "left"),
					h4(textOutput("facility_name"), align = "center", style="background-color: #FFEFCE;color: #000000;padding: 5px;margin-top: 20px"),

					span(tags$i(h6("All reported values are 3-day rolling means, and adjusted for average daily flow rate.")), style="color:#045a8d; text-align: center"),

					plotlyOutput("plot_wwtp", height="250px", width="100%"),

					div(
						style="background-color:#f0f0e3; margins: 1px 0px 10px 0px; padding: 8px;",
						span(tags$i(h5("Facility Information")), style="color:#000000; text-align: center;"),
						h6(textOutput("facility_samples"), align = "left"),
						h6(textOutput("facility_capacity"), align = "left"),
						h6(textOutput("facility_popserved"), align = "left"),
						h6(textOutput("facility_counties"), align = "left")
					),
		
					div(
						style="background-color:#ffffff; margins: 1px 0px 5px 0px; padding: 3px; float: right;",
						actionBttn(inputId="embiggen_open_wwtp", label="Enlarge This Plot", style="pill", size="xs", color="success")
					)
				), # absolutePanel
				
				absolutePanel(
					id = "targets_popup_wwtp_min",
					style = "text-align:center; background-color: #000000; color: #ffffff; border: 1px solid #dddddd; border-radius: 5px; opacity: 0.9;",
					bottom = 0, left = 5, width = 80, height = 25,
					fixed=TRUE, draggable = FALSE,
					span(tags$i("Targets"))
				),					
				absolutePanel(
					id = "targets_popup_wwtp", 
					style="background-color: #000000;color: #ffffff;opacity: 1.0;padding: 3px; border: 1px solid #dddddd; border-radius: 5px; font-size: 11px;", 
					bottom = 25, left = 25, width = 300, height = 130,
					fixed=TRUE, draggable=FALSE, 
					h5("Choose the targets to display:", align = "center", style="background-color: #000000;color: #ffffff;padding: 2px;margin-top: 5px"),
					prettyCheckboxGroup(
						inputId = "targets_wwtp",
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
					id = "dates_popup_wwtp_min",
					style = "text-align:center; background-color: #000000; color: #ffffff; border: 1px solid #dddddd; border-radius: 5px; opacity: 0.9;",
					bottom = 0, left = 90, width = 100, height = 25,
					fixed=TRUE, draggable = FALSE,
					span(tags$i("Date Range"))
				),					
				absolutePanel(
					id = "dates_popup_wwtp", 
					style="background-color: #000000;opacity: 0.8;padding: 3px; border: 1px solid #dddddd; border-radius: 5px; font-size: 11px;", 
					bottom = 25, left = 110, width = 300, height = 130,
					fixed=TRUE, draggable=FALSE, 
					h5("Choose a range of dates to display:", align = "center", style="background-color: #000000;color: #ffffff;padding: 2px;margin-top: 5px"),
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
				),
				absolutePanel(
					id = "roll_popup_wwtp_min",
					style = "text-align:center; background-color: #000000; color: #ffffff; border: 1px solid #dddddd; border-radius: 5px; opacity: 0.9;",
					bottom = 0, left = 195, width = 115, height = 25,
					fixed=TRUE, draggable = FALSE,
					span(tags$i("Data Smoothing"))
				),					
				absolutePanel(
					id = "roll_popup_wwtp", 
					style="background-color: #000000;color: #ffffff;opacity: 1.0;padding: 3px; border: 1px solid #dddddd; border-radius: 5px; font-size: 11px;", 
					bottom = 25, left = 215, width = 300, height = 130,
					fixed=TRUE, draggable=FALSE, 
					h5("Choose a rolling mean window (days):", align = "center", style="background-color: #000000;color: #ffffff;padding: 2px;margin-top: 5px"),
					sliderInput(
						inputId = "roll_wwtp",
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
				) # absolutePanel

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
			), # div
			hidden(
				div(
					id = "conditionalPanelWWTP",
					absolutePanel(
						id = "controls_big_wwtp",
						class = "panel panel-default",
						style="background-color:#ffffff; padding: 4px; margins: 0;border: 2px solid #909090; border-radius: 5px;",
						top = 67, left = 50, width = "90%", height = "85%", fixed=TRUE,
#						draggable = TRUE,
						fluidRow(
							column(
								width = 3,
								div(
									style="padding-top: 15px; padding-left: 5px; float: left;",
									materialSwitch(inputId = "plot_ci_wwtp", label = "90% CI:", value = FALSE, status="warning")
								)
							),
							column(
								width=6,
								h4(textOutput("facility_name_big"), align = "center", style="background-color: #ffffff;color: #000000;padding: 5px; font-weight: bold;"),
							),
							column(
								width=3,
								div(
									style="padding-top: 15px; padding-right: 5px; float: right;",
									actionBttn(inputId="embiggen_close_wwtp", label="Close", style="pill", size="xs", color="warning")
								)
							)
						),
						span(tags$i(h6("All reported values are adjusted for average daily flow rate.")), style="color:#045a8d; text-align: left"),
						plotlyOutput("plot_wwtp_big", height="auto", width="auto")
					) # absolutePanel
				) # hidden panel div
			) # hidden panel
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

