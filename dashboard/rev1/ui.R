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
					top = 162, left = 6, width = 600, fixed=TRUE,
					draggable = FALSE, height = 470,
					div(textOutput("plot_title"), style="padding-top: 10px; font-size: 17px; font-style: bold; color:#045a8d; text-align: center;"),
					plotlyOutput("watch_plot", height="380px", width="100%"),
					span(textOutput("data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
				), # absolutePanel

				absolutePanel(
					id = "controls",
					class = "control_panel",
					top = 633, left = 6, height = "auto", width = 600,
					fixed=TRUE, draggable = FALSE,
					fluidRow(
						column(
							width = 4,
							#style = "border-right: 1px solid #666666;",
							div(
								class = "control_group",
								div(
									style="font-size: 13px;font-weight: 800;margin-bottom: 22px;text-align: left;",
									span("Choose an infection to track:")
								),
								selectInput(
									"infections_control",
									label = NULL,
									choices = INFECTIONS, 
									selected = INFECTIONS_DEFAULT
								)
							)
#							div(
#								style = "margin-left: 20px;font-size: 12px;",
#								prettyCheckboxGroup(
#									inputId = "targets_control",
#									label = NULL,
#									choiceNames = TARGETS,
#									choiceValues = TARGET_VALUES,
#									icon = icon("check-square"), 
#									status = "primary",
#									selected = TARGETS_DEFAULT,
#									outline = TRUE
#								)
#							)
						),
						column(
							width = 8,
							div(
								class = "control_group",
								div(
									style="font-size: 13px;font-weight: 800;margin-bottom: 22px;text-align: left;",
									span("Choose a range of dates to display:")
								),
								sliderTextInput(
									inputId = "dates_control",
									force_edges = TRUE,
									width = "90%",
									label=NULL,
									choices = sort(unique(df_watch$week_starting)),
									selected = c(min(df_watch$week_starting), max(df_watch$week_starting)),
				#							animate=animationOptions(interval = 3000, loop = FALSE),
									grid = FALSE
								)
							)
#						),
#						column(
#							width = 4,
#							div(
#								class = "control_group",
#								div(
#									style="font-size: 13px;font-weight: 800;margin-bottom: 20px;text-align: center;",
#									span("Choose a rolling mean (days)")
#								),
#								selectInput(
#									"roll_control",
#									label = NULL,
#									choices = SMOOTHER_OPTIONS, 
#									selected = SMOOTHER_DEFAULT
#								),
#								materialSwitch(inputId = "ci_control", label = "Show 90% CI:", value = FALSE, status="primary")
#							)
						)
					) # fluidRow
				),					

				absolutePanel(
					class="mdblock",
					id = "site_status_panel", 
#					top = 145, left = 590, height=120, width = 130,
					top = 145, left = 545, height=70, width = 200,
					fixed=TRUE, draggable=FALSE,
					span(textOutput("site_change_hdr"), style="font-size: 13px; color: #000000; line-height: 22px;"),
					span(textOutput("site_change_above"), style="font-size: 14px; font-weight: 800; color: #000000; line-height: 20px;"),
					span(textOutput("site_change_below"), style="font-size: 12px; color: #000000; line-height: 22px;")
				),
				
				hidden(
					absolutePanel(
						id = "site_change_info",
						class = "mdinfo",
						top = 215, left = 545, height="auto", width = 470,
						fixed=TRUE, draggable=TRUE,
						div(textOutput("change_plot_title"), style="padding-top: 10px; font-size: 13px; font-style: normal; color:#045a8d; text-align: center;"),
						plotlyOutput("change_plot", height="318px", width="100%"),
						div("Columns show the fold-change of the signal for that day, compared to lowest recorded value for this site. The line shows the 5-day rolling mean of fold change. The colored boxes indicate the approximate community transmission, which is calculated from the 5-day rolling mean.", style="font-size: 11px; color:#888888; text-align: center;"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="site_change_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "site_trend_panel", 
					top = 145, left = 745, height=70, width = 200,
					fixed=TRUE, draggable=FALSE,
					span(textOutput("site_trend_hdr"), style="font-size: 13px; line-height: 18px;"),
					span(textOutput("site_trend_txt"), style="font-size: 18px; font-weight: 800; line-height: 28px;"),
					span(textOutput("site_trend_level"), style="font-size: 12px; line-height: 18px;"),
				),
				
				hidden(
					absolutePanel(
						id = "site_focus_info",
						class = "mdinfo",
						top = 215, left = 545, height="auto", width = 470,
						fixed=TRUE, draggable=TRUE,
						div(textOutput("focus_plot_title"), style="padding-top: 10px; font-size: 13px; font-style: normal; color:#045a8d; text-align: center;"),
						plotlyOutput("focus_plot", height="318px", width="100%"),
						span(textOutput("focus_data_format"), style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="site_focus_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock-nopop",
					id = "scope_panel", 
					top = 65, left = 45, height=75, width = 240,
					fixed=TRUE, draggable=FALSE,
					span(textOutput("scope"), style="font-size: 19px; font-weight: 800; line-height: 26px;"),
					span(textOutput("population_served"), style="font-size: 15px;font-weight: 200;"),
					span(textOutput("sample_count"), style="font-size: 11px; font-weight: 200;")
				),
				
				absolutePanel(
					class="mdblock-nopop",
					id = "population_panel", 
					top = 65, left = 285, height=75, width = 130,
					fixed=TRUE, draggable=FALSE,
					style = "padding-top: 8px;",
					span(textOutput("counties_served"), style="font-size: 14px;"),
					span("county population", style="font-size: 10px;"),
					span(textOutput("county_population"), style="font-size: 15px;")
				),
				
				absolutePanel(
					class="mdblock",
					id = "daily_flow_panel", 
					top = 65, left = 415, height=75, width = 130,
					fixed=TRUE, draggable=FALSE,
					span("Mean daily flow", style="font-size: 13px;"),
					span(textOutput("mean_flow"), style="font-size: 16px;"),
					span("(million gallons per day)", style="font-size: 11px;")
				),
				
				hidden(
					absolutePanel(
						id = "daily_flow_info",
						class = "mdinfo",
						top = 215, left = 545, height="auto", width = 470,
						fixed=TRUE, draggable=TRUE,
						plotlyOutput("flow_plot", height="300px", width="100%"),
						span(textOutput("flow_plot_title"), style="font-size: 14px; font-style: bold; color:#045a8d; text-align: center;"),
						div("Daily flow is provided by the facility with each sample, in million gallons per day.", style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="daily_flow_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),
					
				absolutePanel(
					class="mdblock",
					id = "collection_panel", 
					top = 65, left = 545, height=75, width = 160,
					fixed=TRUE, draggable=FALSE,
					span("Last data update was", style="font-size: 13px;"),
					span(textOutput("last_update"), style="font-size: 19px;"),
					span(textOutput("last_update_stamp"), style="font-size: 11px;")
				),
				
				hidden(
					absolutePanel(
						id = "collection_info",
						class = "mdinfo",
						top = 215, left = 545, height="auto", width = 470,
						fixed=TRUE, draggable=TRUE,
						plotlyOutput("collection_plot", height="300px", width="100%"),
						span(textOutput("collection_frequency"), style="font-size: 14px; font-style: bold; color:#045a8d; text-align: center;"),
						#div("Daily flow is provided by the facility with each sample, in million gallons per day.", style="font-size: 12px; font-style: italic; color:#888888; text-align: center;"),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="collection_info_close", label="Close", style="pill", size="xs", color="success")
						)
					)
				),

				absolutePanel(
					class="mdblock",
					id = "alert_panel", 
					top = 65, left = 707, height=75, width = 238,
					fixed=TRUE, draggable=FALSE,
					span(textOutput("alert_hdr"), style="font-size: 13px; color: #000000; line-height: 22px;"),
					span(textOutput("alert_change"), style="font-size: 20px; font-weight: 800; color: #000000; line-height: 20px;"),
					span(textOutput("alert_trend"), style="font-size: 16px; font-weight: 800; color: #000000; line-height: 22px;")
				),
				
				hidden(
					absolutePanel(
						id = "alert_level_info",
						class = "mdinfo",
						top = 215, left = 545, height="auto", width = 450,
						fixed=TRUE, draggable=TRUE,
						div(
							class = "alertinfo",
							span("XXX", style="height: 50px; width: 50px; margin: 5px; color: #c5000b; background-color: #c5000b; border: 1 px solid black; border-radius: 3px;"),
							span("Red. The total amount of infection in the wastewater is more than 500X above the lowest reported value. Community transmission is extremely high.", style="font-size: 14px;")
						),
						div(
							class = "alertinfo",
							span("XXX", style="height: 50px; width: 50px; margin: 5px; color: #ff950e; background-color: #ff950e; border: 1 px solid black; border-radius: 3px;"),
							span("Orange. Wastewater levels are 250-500X above the lowest value for this infection. Community transmission is very high.", style="font-size: 14px;")
						),
						div(
							class = "alertinfo",
							span("XXX", style="height: 50px; width: 50px; margin: 5px; color: #eedc82; background-color: #eedc82; border: 1 px solid black; border-radius: 3px;"),
							span("Gold. Wastewater levels are 100-250X above the lowest value for this infection. Community transmission is high.", style="font-size: 14px;")
						),
						div(
							class = "alertinfo",
							span("XXX", style="height: 50px; width: 50px; margin: 5px; color: #ffd320; background-color: #ffd320; border: 1 px solid black; border-radius: 3px;"),
							span("Yellow. Wastewater levels are 25-100X above the lowest value for this infection. Community transmission is moderate.", style="font-size: 14px;")
						),
						div(
							class = "alertinfo", 
							span("XXX", style="height: 50px; width: 50px; margin: 5px; color: #579d1c; background-color: #579d1c; border: 1 px solid black; border-radius: 3px;"),
							span("Blue. Wastewater levels are 0-25X above the lowest reported value for this infection. Community transmission is limited.", style="font-size: 14px;")
						), # div
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="alert_level_info_close", label="Close", style="pill", size="xs", color="success")
						) # button div
					)
				),
					
				absolutePanel(
					id = "legend_panel", 
					class = "card", 
					bottom = 70, right = 5, width = 125, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					style = "background-color: #ffffff;border: 1px solid #000000;border-radius: 3px;",
					div("Community Spread", style="text-align: center; font-weight: 800; font-size: 13px; margin-bottom: 5px;"),
					div(
						span("XX", style="height: 20px; width: 20px; margin: 5px; color: #c5000b; background-color: #c5000b; border: 1 px solid black; border-radius: 3px;"),
						span("Extremely high", style="font-size: 11px;")
					),
					div(
						span("XX", style="height: 20px; width: 20px; margin: 5px; color: #ff950e; background-color: #ff950e; border: 1 px solid black; border-radius: 3px;"),
						span("Very high", style="font-size: 11px;")
					),
					div(
						span("XX", style="height: 20px; width: 20px; margin: 5px; color: #eedc82; background-color: #eedc82; border: 1 px solid black; border-radius: 3px;"),
						span("High", style="font-size: 11px;")
					),
					div(
						span("XX", style="height: 20px; width: 20px; margin: 5px; color: #ffd320; background-color: #ffd320; border: 1 px solid black; border-radius: 3px;"),
						span("Moderate", style="font-size: 11px;")
					),
					div(
						span("XX", style="height: 20px; width: 20px; margin: 5px; color: #579d1c; background-color: #579d1c; border: 1 px solid black; border-radius: 3px;"),
						span("Low", style="font-size: 11px;")
					)
				), # absolutePanel
				
				absolutePanel(
					id = "recenter_panel", 
					class = "card", 
					bottom = 28, right = 0, width = 125, 
					fixed = TRUE, draggable = FALSE, 
					height = "auto",
					actionBttn(inputId="center_map", label="Recenter Map", style="pill", size="xs", color="success")
				), # absolutePanel
				
				absolutePanel(
					id = "logo", 
					class = "card", 
					bottom = 5, left = 160, width = 80, 
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

#		tabPanel(
#			"Insights",
#			div(
#				class="outer",
#				
#				absolutePanel(
#					id = "insights-2way", 
#					class = "panel panel-default",
#					top = 65, left = 6, width = 400, fixed=TRUE,
#					span(textOutput("plot_2way_title"), style="font-size: 14px; font-style: bold; color:#045a8d; text-align: center;"),
#					plotlyOutput("insights_plot_2way", height="380px", width="100%")
#				), # absolutePanel
#				
##				absolutePanel(
##					id = "insights-2way-qq1", 
##					class = "panel panel-default",
##					top = 65, left = 406, width = 200, fixed=TRUE,
##					plotlyOutput("insights_plot_2wayQQ_1", height="200px", width="100%")
##				), # absolutePanel
##				
##				absolutePanel(
##					id = "insights-2way-qq2", 
##					class = "panel panel-default",
##					top = 265, left = 406, width = 200, fixed=TRUE,
##					plotlyOutput("insights_plot_2wayQQ_2", height="200px", width="100%")
##				), # absolutePanel
#				
#				absolutePanel(
#					id = "insights-2way-cor", 
#					class = "panel panel-default",
#					top = 465, left = 6, width = 400, fixed=TRUE,
#					span(textOutput("insights_2way_cor"), style="font-size: 12px; font-style: bold; color:#045a8d; text-align: center;")
#				) # absolutePanel
#				
#			) # div, outer
#		), # tabPanel
#
#		tabPanel(
#			"Ethics of testing",
#			tags$div(
#				tags$h4("Coming soon..."), 
#				#h6(paste0(update)),
##				"This site is updated weekly or biweekly. ", 
##				tags$a(href="https://experience.arcgis.com/experience/685d0ace521648f8a5beeeee1b9125cd", "the WHO,"),
#				"This tab will present information about the ethics of wastewater testing."
#			)
#		), # tabPanel

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

