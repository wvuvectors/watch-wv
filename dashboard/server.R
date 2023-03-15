#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


shinyServer(function(input, output, session) {
	
#	getTrendData <- function(target) {
#		
#	}

	#
	#
	# Create the map
	#
	generateMap <- function(center_lat, center_lng, zoom_level) {
		print("generateMap called!")

		mymap = leaflet() %>% 
		addTiles() %>% 
		setView(lng = center_lng, lat = center_lat, zoom = zoom_level) %>% 
		addPolylines(data=county_sf, fill=FALSE, weight=3, color="#999999", layerId="countiesLayer") %>%
		addPolygons(data=state_sf, fill=FALSE, weight=2, color="#000000", layerId="stateLayer")
		
		return(mymap)
	}


	#
	# Plot of pathogen
	#
	plotIt <- function(targ, locus) {
		print("plotIt called!")
		
		df_plot <- df_result %>% filter(target == targ & target_genetic_locus == locus) %>% 
		           group_by(date_to_plot) %>% 
		           arrange(date_to_plot) %>%
							 summarize(val := mean(target_copies_per_l, na.rm = TRUE))

		gplot <- ggplot(df_plot) + labs(y = "", x = "") + 
											scale_y_continuous(labels = comma) + 
											scale_x_date() + 
											#scale_color_manual(name = "Target", values = TARGET_COLORS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											#scale_fill_manual(name = "Target", values = TARGET_FILLS, labels = c("n1" = "SARS-CoV-2 N1", "n1n2" = "SARS-CoV-2 N1N2", "n2" = "SARS-CoV-2 N2")) + 
											plot_theme() + 
											geom_point(aes(x = date_to_plot, y = val, color = val), shape = 1, size = 2, alpha=0.9) + 
											geom_col(aes(x = date_to_plot, y = val, fill = val), alpha=0.4, na.rm = TRUE)

		ggplotly(gplot)# %>% layout(legend = list(orientation = "h"))
	}
	

	#
	# Render the WW table
	#
	output$tableWW <- renderDataTable(
		#df_result %>% filter(target == "SARS-CoV-2") %>% 
		#select(date_to_plot, location_id, target_genetic_locus, target_copies_flownorm_per_person) %>% 
		#rename(Date = date_to_plot, Location = location_id, Locus = target_genetic_locus, Copies = target_copies_flownorm_per_person) %>% 
		df_county %>% filter(county_lab_code != "") %>% select(county_name, county_population, county_lab_code) %>% rename(county = county_name, population = county_population, status = county_lab_code),
		options = list(paging = FALSE,   ## paginate the output
									 pageLength = 9,   ## number of rows to output for each page
									 scrollX = TRUE,   ## enable scrolling on X axis
									 scrollY = TRUE,   ## enable scrolling on Y axis
									 autoWidth = TRUE, ## use smart column width handling
									 server = FALSE,   ## use server- or client-side processing
									 dom = 't',
									 #buttons = c('csv', 'excel'),
									 columnDefs = list(list(targets = '_all', className = 'dt-center'))
		),
		#extensions = 'Buttons',
		selection = 'single',						## enable selection of a single row
		#filter = 'bottom',            	## include column filters at the bottom
		rownames = FALSE                ## don't show row numbers/names
	)	


	output$plotWW <- renderPlotly({
		#print("watch_plot top")

		#writeMetadata()
		#print("writeMetadata() just ran!")

		plotIt("SARS-CoV-2", "N2") %>% config(displayModeBar = FALSE) %>% style(hoverinfo = "skip")
		#print("watch_plot bottom")
	})	

	#
	# Render the map
	#
	output$watch_map <- renderLeaflet({
		center <- MAP_CENTERS %>% filter(layer == "WWTP")
		generateMap(center_lat = center$lat, center_lng = center$lng, zoom_level = center$zoom)
	})




})

