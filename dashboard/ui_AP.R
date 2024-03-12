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
#			div(
#				class="outer",

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
					actionBttn(inputId="center_map_rs", label="Recenter Map", style="pill", size="xs", color="default")
				), # absolutePanel
				

				absolutePanel(
					class = "controls",
					top = 65, left = 235, height = 95, width = 400,
					fixed=TRUE, draggable=FALSE,
					div("MAP & PLOT CONTROLS", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
					div(
						class = "map_embed",
						style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:125px;",
						"Map View",
						selectInput(
							"geo_level_rs",
							label = NULL,
							choices = GEOLEVELS, 
							selected = GEOLEVELS_DEFAULT
						)
					),
					div(
						class = "map_embed",
						style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:125px;",
						"Map Color",
						selectInput(
							"map_color_rs",
							label = NULL,
							choices = c("By Lab" = "lab"),
#							choices = c("By Lab" = "lab", "FluA" = "FluA", "SARS-CoV-2" = "SARS-CoV-2", "RSV" = "RSV"),
							selected = "lab"
						)
					),
					div(
						class = "map_embed",
						style = "display: inline-block;font-size: 14px;font-weight: 800;text-align: center;width:125px;",
						"Plot View",
						selectInput(
							"view_range_rs",
							label = NULL,
							choices = c("1 month" = 1, "3 months" = 3, "6 months" = 6, "1 year" = 12, "2 years" = 24),
							selected = VIEW_RANGE_PRIMARY
						)
# 					),
# 					div(
# 						class = "map_embed",
# 						style = "font-size: 15px;font-weight: 800;text-align: left;padding: 0px;color: #FF2600;",
# 						checkboxInput(
# 							inputId = "trendline_1mo_rs",
# 							label = "Show 1 Month Mean Line (green)",
# 							value = TRUE
# 						)
# 					),
# 					div(
# 						class = "map_embed",
# 						style = "font-size: 15px;font-weight: 800;text-align: left;padding: 0px;color: #0096FF;",
# 						checkboxInput(
# 							inputId = "trendline_1yr_rs",
# 							label = "Show 1 Year Mean Line (gold)",
# 							value = TRUE
# 						)
					)
				), # absolutePanel

				# SARS seqr data
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 590, left = 5, width = 637, height = 290, 
					div(textOutput("plotsq_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
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
					div(id="rs_CSL_label", "Alert Level: ", style="display: inline-block;font-size: 13px;font-weight: 800;text-align: right;padding: 2px;width:145px;"),
					div(id="rs_CSL_1", textOutput("rs_CSL_1"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin-right: 1px;margin-left: 3px;padding: 2px;width:83px;background-color: #ee4444;"),
					div(id="rs_CSL_2", textOutput("rs_CSL_2"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #44ee44;"),
					div(id="rs_CSL_3", textOutput("rs_CSL_3"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #4444ee;"),
					div(id="rs_CSL_4", textOutput("rs_CSL_4"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #eeee44;")
				), # absolutePanel
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 138, left = 648, width = 505, height = 30, 
					div(id="rs_delta_label", "Change from Yr. Avg.:", style="display: inline-block;font-size: 13px;font-weight: 800;text-align: right;padding: 2px;width:145px;"),
					div(id="rs_delta_1", textOutput("rs_delta_1"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin-right: 1px;margin-left: 3px;padding: 2px;width:83px;background-color: #fafafa;"),
					div(id="rs_delta_2", textOutput("rs_delta_2"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #fafafa;"),
					div(id="rs_delta_3", textOutput("rs_delta_3"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #fafafa;"),
					div(id="rs_delta_4", textOutput("rs_delta_4"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #fafafa;")
				), # absolutePanel
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 170, left = 648, width = 505, height = 30, 
					div(id="rs_fresh_label", "Data Age (days):", style="display: inline-block;font-size: 13px;font-weight: 800;text-align: right;padding: 2px;width:145px;"),
					div(id="rs_fresh_1", textOutput("rs_fresh_1"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin-right: 1px;margin-left: 3px;padding: 2px;width:83px;background-color: #fafafa;"),
					div(id="rs_fresh_2", textOutput("rs_fresh_2"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #fafafa;"),
					div(id="rs_fresh_3", textOutput("rs_fresh_3"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #fafafa;"),
					div(id="rs_fresh_4", textOutput("rs_fresh_4"), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: center;margin: 1px;padding: 2px;width:83px;background-color: #fafafa;")
				), # absolutePanel

				absolutePanel(
					id = "alert_scale", 
					class = "panel panel-default",
					style = "background-color: #000000;",
					fixed = TRUE, draggable = FALSE, 
					top = 202, left = 648, width = 125, height = 150, 
					div("Alert Color Key", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;color: #FFFFFF;"),
					div("LOW", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:119px;background-color: #3288BD;"),
					div("MODERATE", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:119px;background-color: #E6F598;"),
					div("HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:119px;background-color: #FDAE61;"),
					div("VERY HIGH", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:119px;background-color: #D53E4F;"),
					div("UNKNOWN", style="font-size: 12px;font-weight: 800;text-align: center;padding: 2px;margin: 3px 2px;width:119px;background-color: #EEEEEE;")
				), # absolutePanel
				absolutePanel(
					class = "panel panel-default",
					style = "background-color: #000000;",
					fixed = TRUE, draggable = FALSE, 
					top = 202, left = 775, width = 378, height = 150, 
					div(textOutput("site_rs_info"), style="font-size: 14px; font-weight: 400;text-align: center;padding: 15px 2px; color: #ffffff;"),
					div(textOutput("site_rs_flow"), style="font-size: 14px; font-weight: 400;text-align: center;padding: 0px 2px; color: #ffffff;")
#					div(downloadButton("downloadData", "Download"), style="font-size: 14px; font-weight: 400;text-align: center;padding-top: 15px;")
#					div(downloadLink(outputId = "downloadData", label = NULL, class = "dlsmall"), "Download (csv)", style="display: inline-block; font-size: 9px; font-weight: 400;text-align: left;width: 25px;height: 25px;padding-top: 5px;")
				), # absolutePanel

				# Top right
				absolutePanel(
					class = "panel panel-default",
					style = "background-color: #000000;",
					fixed = TRUE, draggable = FALSE, 
					top = 63, left = 1154, width = 345, height = 58, 
					div("Percent Change from Yearly Mean (All Counties)", style="font-size: 13px;font-weight: 800;text-align: center;padding: 2px;color: #FFFFFF;"),
					div("County", style="display: inline-block;font-size: 15px;font-weight: 800;text-align: left;margin-right: 1px;margin-left: 3px;padding: 2px;width:96px;background-color: #000000; color: #ffffff;"),
					div(paste0(DISEASE_RS[1]), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: left;margin: 1px;padding: 2px;width:56px;background-color: #000000; color: #ffffff;"),
					div(paste0(DISEASE_RS[2]), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: left;margin: 1px;padding: 2px;width:54px;background-color: #000000; color: #ffffff;"),
					div(paste0(DISEASE_RS[3]), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: left;margin: 1px;padding: 2px;width:65px;background-color: #000000; color: #ffffff;"),
					div(paste0(DISEASE_RS[4]), style="display: inline-block;font-size: 15px;font-weight: 800;text-align: left;margin: 1px;padding: 2px;width:46px;background-color: #000000; color: #ffffff;")
				),
				absolutePanel(
					class = "panel panel-default",
					style = "overflow: auto;",
					fixed = TRUE, draggable = FALSE, 
					top = 122, left = 1154, width = 345, height = 230, 
					tableHTML(
						df_regions %>% filter(region_geolevel == "county") %>% select(region_name, FLUA, FLUB, COVID, RSV) %>% rename(County = region_name), 
						#collapse = "separate_shiny", 
						#spacing = "5px 2px", 
						rownames = FALSE, 
						border = 0,
						widths = c(96, 60, 60, 65, 60)
					) %>% 
					add_css_thead(css = list("background-color", "#000000")) %>% 
					add_css_thead(css = list("color", "#000000")) %>% 
					add_css_thead(css = list("font-size", "1px")) %>% 
					add_css_column(css = list("text-align", "center"), columns=names(df_regions[2:5])) %>% 
					add_css_table(css = list("background-color", "#ffffff")) %>% 
					add_css_row(css = list("background-color", "#f2f2f2"), rows = odd(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1)) %>%
					add_css_row(css = list("background-color", "#e6f0ff"), rows = even(1:length((df_regions %>% filter(region_geolevel == "county"))$region_name)+1))
				), # absolutePanel


				# Middle left
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 358, left = 648, width = 425, height = 225, 
					div(textOutput("plot1_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
					plotlyOutput("plot1_rs", height="210px", width="100%")
				), # absolutePanel

				# Middle right
				absolutePanel(
					class = "panel panel-default",
					fixed = TRUE, draggable = FALSE, 
					top = 358, left = 1075, width = 425, height = 225, 
					div(textOutput("plot2_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
					plotlyOutput("plot2_rs", height="210px", width="100%")
				), # absolutePanel

				# Bottom left
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 590, left = 648, width = 425, height = 225, 
					div(textOutput("plot3_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
					plotlyOutput("plot3_rs", height="210px", width="100%")
				), # absolutePanel

				# Bottom right
				absolutePanel(
					class = "panel panel-default",
					fixed=TRUE, draggable = FALSE,
					top = 590, left = 1075, width = 425, height = 225, 
					div(textOutput("plot4_rs_title"), style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;color: #ffffff; background-color: #303D4E;"),
					plotlyOutput("plot4_rs", height="210px", width="100%")
				), # absolutePanel

				# Tag line (below the plot grid)
				absolutePanel(
					class = "panel panel-default",
					style = "background-color: #303D4E;",
					fixed=TRUE, draggable = FALSE,
					top = 815, left = 648, width = 505, height = 66, 
					div(
						style = "font-size: 14px;padding: 2px;font-weight: 400;font-style: normal;text-align: left;color: #ffffff;",
						"Plot values are shown as copies of viral particles per person, after correction for daily flow. Variant proportions are shown as percent of the total identified SARS variants."
					)
				), # absolutePanel

				absolutePanel(
					class = "panel panel-default",
					style = "background-color: #303D4E;",
					fixed=TRUE, draggable = FALSE,
					top = 815, left = 1154, width = 345, height = 66, 
					div(
						style = "font-size: 12px;padding: 2px;font-weight: 400;font-style: normal;text-align: left;color: #ffffff;",
						span("The "),
						span(style="color: #EAAA00;", "solid gold line"),
						span("on each plot represents the annual mean level of each target, used to generate the percent change. The "),
						span(style="color: #00B140;", "green dashed line"),
						span("represents the mean over the most recent month.")
					)
				), # absolutePanel
				
				hidden(
					absolutePanel(
						id = "alert_scale_info",
						class = "mdinfo",
						top = 150, left = 775, width = 450, height = 350,
						div("Alert Color Explanations", style="font-size: 13px;padding: 4px;font-weight: 800;text-align: center;margin-bottom: 5px;color: #ffffff; background-color: #303D4E;"),
						div(
							class = "alertinfo",
							id = "level_1",
							span("XX", style="height: 40px; width: 40px; margin: 5px; color: #3288BD; background-color: #3288BD; border: 1px solid black; border-radius: 3px;"),
							span("The latest amount of this disease agent is less than 50% of the yearly average. Community transmission is estimated to be low.", style="font-size: 13px;")
						),
						div(
							class = "alertinfo",
							id = "level_2",
							span("XX", style="height: 40px; width: 40px; margin: 5px; color: #E6F598; background-color: #E6F598; border: 1px solid black; border-radius: 3px;"),
							span("The latest amount of this disease agent is between 50% and 100% of the yearly average. Community transmission is estimated to be moderate.", style="font-size: 13px;")
						),
						div(
							class = "alertinfo",
							id = "level_3",
							span("XX", style="height: 40px; width: 40px; margin: 5px; color: #FDAE61; background-color: #FDAE61; border: 1px solid black; border-radius: 3px;"),
							span("The latest amount of this disease agent is between 100% and 150% of the yearly average. Community transmission is estimated to be high.", style="font-size: 13px;")
						),
						div(
							class = "alertinfo",
							id = "level_4",
							span("XX", style="height: 40px; width: 40px; margin: 5px; color: #D53E4F; background-color: #D53E4F; border: 1px solid black; border-radius: 3px;"),
							span("The latest amount of this disease agent is greater than 150% of the yearly average. Community transmission is estimated to be very high.", style="font-size: 13px;")
						),
						div(
							class = "alertinfo",
							id = "level_5",
							span("XX", style="height: 40px; width: 40px; margin: 5px; color: #EEEEEE; background-color: #EEEEEE; border: 1px solid black; border-radius: 3px;"),
							span("The data for this disease agent is either missing or too old to make an accurate determination.", style="font-size: 14px;")
						),
						div(
							style="padding-top: 15px; padding-right: 5px; float: right;",
							actionBttn(inputId="alert_scale_info_close", label="Close", style="pill", size="xs", color="success")
						) # button div
					)
				)
					

				# Alert message boxes (one for each plot)
# 				conditionalPanel(
# 					id = "alert_info_4",
# 					condition = "plot3_rs_err.output == 'y'",
# 					class = "mdinfo",
# 					top = 640, left = 1100, height="auto", width = 150,
# 					div(
# 						class = "alertinfo", 
# 						style="font-size: 14px;text-align: center;", 
# 						"No data available."
# 					)
# 						div(
# 							style="padding-top: 15px; padding-right: 5px; float: right;",
# 							actionBttn(inputId="alert_level_info_close", label="Close", style="pill", size="xs", color="success")
# 						) # button div
# 				)
				
#			) # tabPanel div
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
					actionBttn(inputId="center_map_mg", label="Recenter Map", style="pill", size="xs", color="default")
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

