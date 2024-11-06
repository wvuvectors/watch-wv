shinyUI(fluidPage(
	useShinyjs(),
	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "main_styles.css")
  ),
	
	tags$div(class="nbcontainer", 
	
	navbarPage(
		theme = shinytheme("flatly"), 
		id="nav",
		HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Wastewater Testing for Community Health in WV (WaTCH-WV)</a>'), 
		windowTitle = "WaTCH-WV",
		collapsible = FALSE,

		tabPanel(
			"COVID",
			fluidRow(
				column(5, 
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						column(12,
							leafletOutput("map_covid", width = "100%", height = "443px")
						)
					), # fluidRow (map)
					fluidRow(
						column(12,
							div("MAP & PLOT CONTROLS", 
							style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
						)
					), # fluidRow (controls title)
					fluidRow(
						column(8,
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center; width:105px;",
								"Map View",
								selectInput(
									"geo_level_covid",
									label = NULL,
									choices = GEOLEVELS, 
									selected = GEOLEVELS_DEFAULT
								)
							),
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center; width:105px;",
								"Map Color",
								selectInput(
									"map_color_covid",
									label = NULL,
									choices = c("Lab" = "lab"),
			#							choices = c("Status" = "Status", "Risk" = "Risk", "Freshness" = "Freshness"),
									selected = "lab"
								)
							),
							div(
								class = "map_embed",
								style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center; width:120px;",
								"Plot View",
								selectInput(
									"view_range_covid",
									label = NULL,
									choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
									selected = VIEW_RANGE_PRIMARY
								)
							)
						),
						column(4,
							style = "text-align: right;",
							div(
								style = "display: inline-block;margin-top: 25px;",
								actionBttn(inputId="map_reset_covid", label="Reset Map", style="pill", size="xs", color="success")
							)
						)
					), # fluidRow (controls)
					fluidRow(
						#style = "border: 2px solid #941100;",
						column(12,
							# title
							div(
								textOutput("selection_title_covid"), 
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E"
							)
						)
					), # fluidRow (selection title)
					fluidRow(
						# alert key & site info
						style = "padding-top: 4px;padding-bottom: 4px;margin-left: 0px; margin-right: 0px;background-color: #000000;color: #ffffff;",
						column(3,
							div(
								class = "logo",
								id = "risk_level_key_covid", 
								#class = "panel panel-default",
								style = "margin 5px;padding: 3px;text-align: center;",
								div("Risk Level Key", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;width:100px;"),
								div("VERY HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:100px;color: #000000;background-color: #D53E4F;"),
								div("HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:100px;color: #000000;background-color: #FDAE61;"),
								div("MODERATE", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:100px;color: #000000;background-color: #E6F598;"),
								div("LOW", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:100px;color: #000000;background-color: #3288BD;"),
								div("UNKNOWN", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:100px;color: #000000;background-color: #EEEEEE;")
							),
							div(
								class = "logo",
								style = "padding-top: 10px;text-align: center;",
								tags$a(href='https://www.watch-wv.com/', 
								tags$img(src='WaTCH-WV_logo.png',height='75',width='75'))
							)
						),
						column(9,
							fluidRow(
								column(12,
									div(
										textOutput("selection_details_covid"), 
										style="font-size: 14px;font-weight: 800;text-align: center;padding-bottom: 5px;"
									)
								)
							), # fluidRow (selection details)
							fluidRow(
								column(12,
									div(
										textOutput("selection_samples_covid"), 
										style="font-size: 16px;padding: 10px 0px;font-weight: 800;text-align: center;color: #FFFC79;"
									)
								)
							), # fluidRow (abundance data freshness)
							fluidRow(
								column(12,
									div(
										textOutput("selectionsq_samples_covid"), 
										style="font-size: 16px;padding: 10px 0px;font-weight: 800;text-align: center;color: #FFFC79;"
									)
								)
							), # fluidRow (variant data freshness)
							fluidRow(
								style = "margin-left: 0px; margin-right: 0px;background-color: #000000; color: #ffffff;",
								column(12,
									style = "font-size: 13px;padding: 15px 0px;font-weight: 400;text-align: center;",
									div("WaTCH-WV is supported by CDC-sponsored grants from the West Virginia Department of Health, to West Virginia University and Marshall University.")
								)
							) # fluidRow (logo and sponsors)
						)	# column
					), # fluidRow (color key & selection info)
				),
				column(5,
					style = "margin-top: 5px;padding: 3px;",
					fluidRow(
						# alert blocks
						column(3,
							div(
								"Risk Level",
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								textOutput("risk_covid"),
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(3,
							div(
								"Abundance",
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								textOutput("abundance_covid"),
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(3,
							div(
								"Trend",
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								textOutput("trend_covid"),
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						),
						column(3,
							div(
								"Variant",
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #000000;"),
							div(
								textOutput("variant_covid"),
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #000000; background-color: #ffffff; border: 2px solid #000000;")
						)
					), # fluidRow (status blocks)
					fluidRow(
						#style = "margin-left: 0px; margin-right: 0px;", 
						style = "margin-top: 5px;", 
						column(12,
							div(
								textOutput("plot_title_covid"), 
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (COVID plot title)
					fluidRow(
						column(12,
							# Plot of COVID change over time
							plotlyOutput("plot_covid", height="358px", width="100%")
						)
					), # fluidRow (COVID plot)
					fluidRow(
						#style = "margin-left: 0px; margin-right: 0px;", 
						style = "margin-top: 5px;", 
						column(12,
							div(
								textOutput("plotsq_title_covid"), 
								style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"
							)
						)
					), # fluidRow (COVID variant plot title)
					fluidRow(
						style = "margin-top: 5px;",
						column(12,
							# Plot of COVID variant proportions over time
							plotlyOutput("plotsq_covid", height="325px", width="100%")
						)
					) # fluidRow (variant plot)
				),
				column(2,
					style = "margin-top: 8px;",
					fluidRow(
						column(4,
							style = "background-color: #303D4E;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
							div("County", style = "padding: 2px;")
						),
						column(2,
							style = "background-color: #303D4E;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
							div("Level", style = "padding: 2px;")
						),
						column(3,
							style = "background-color: #303D4E;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
							div("Trend", style = "padding: 2px;")
						),
						column(3,
							style = "background-color: #303D4E;color: #FFFFFF;text-align: left;font-size: 15px;font-weight: 800;padding: 3px;",
							div("Age (d)", style = "padding: 2px;")
						)
					), # fluidRow (table header)
					fluidRow(
						column(12,
							style = "overflow: auto; padding: 3px;",
							tableHTML(
								df_regions %>% filter(region_geolevel == "county") %>% select(region_name, any_of(DISEASES[1])) %>% rename(County = region_name), 
								#collapse = "separate_shiny", 
								#spacing = "5px 2px", 
								rownames = FALSE, 
								border = 0
								#widths = c(96, 60, 60, 65, 60)
							) %>% 
							add_css_thead(css = list("background-color", "#000000")) %>% 
							add_css_thead(css = list("color", "#000000")) %>% 
							add_css_thead(css = list("font-size", "1px")) %>% 
							add_css_column(css = list("text-align", "left"), columns=names(df_regions[2:5])) %>% 
							add_css_table(css = list("width", "100%")) %>% 
							add_css_table(css = list("background-color", "#ffffff")) %>% 
							add_css_row(css = list("background-color", "#f2f2f2"), rows = odd(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1)) %>%
							add_css_row(css = list("background-color", "#e6f0ff"), rows = even(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1))
						) # column
					), # fluidRow (table)
					fluidRow(
						style = "margin-top: 12px; background-color: #000000; color: #ffffff; padding: 5px;",
						column(12,
							# explanation of plots
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px 20px 2px 2px;text-align: left;",
								"Plot values are reported as the copies of viral particles per person after correction for daily flow, averaged across all selected sites."
							)
						)
					), # fluidRow (abundance plot explanation)
					fluidRow(
						style = "margin-top: 12px; background-color: #000000; color: #ffffff; padding: 5px;",
						column(12,
							# explanation of plots
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px 20px 2px 2px;text-align: left;",
								"Variant proportions are shown as the percent of total identified variants, averaged across all selected sites."
							)
						)
					), # fluidRow (variant plot explanation)
					fluidRow(
						style = "margin-top: 12px; background-color: #000000; color: #ffffff; padding: 5px;",
						column(12,
							# explanation of trend lines
							style = "font-size: 14px;padding: 3px;font-weight: 400;",
							div(
#								style = "padding: 2px;text-align: left;",
								span("The "),
								span(style="color: #EAAA00;", "solid gold line"),
								span(paste0("on each plot represents the ", VIEW_RANGE_PRIMARY, " month average abundance. The ")),
								span(style="color: #00B140;", "green dashed line"),
								span("represents the average abundance over the most recent month.")
							)
						)
					) # fluidRow (trend line descriptions)
				)
			),

			hidden(
				absolutePanel(
					id = "risk_level_key_popup",
					class = "mdinfo",
					top = 365, left = 610, width = 580, height = 350,
					div("Risk Level Color Key", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "alertinfo",
						id = "level_1",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
						span(paste0("LOW. The latest amount of this disease agent is less than 50% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be low."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_2",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
						span(paste0("MODERATE. The latest amount of this disease agent is between 50% and 100% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be moderate."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_3",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
						span(paste0("HIGH. The latest amount of this disease agent is between 100% and 150% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_4",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
						span(paste0("VERY HIGH. The latest amount of this disease agent is greater than 150% of the ", VIEW_RANGE_PRIMARY, " month average. Community transmission is estimated to be very high."), style="font-size: 13px;")
					),
					div(
						class = "alertinfo",
						id = "level_5",
						span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
						span("UNKNOWN. The data for this disease agent is either missing or too old to make an accurate determination.", style="font-size: 13px;")
					),
					div(
						style="padding-top: 15px; padding-right: 5px; float: right;",
						actionBttn(inputId="risk_level_key_popup_close", label="Close", style="pill", size="xs", color="success")
					) # button div
				)
			) # hidden (risk level key popup)

		) # tabPanel (COVID)

	)) # navbarPage and container div
))

