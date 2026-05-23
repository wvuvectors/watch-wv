#! /usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", dependencies = TRUE)
    
BiocManager::install("Biobase")

list.of.packages <- c(
	"tidyverse", 
	"dplyr", 
	"data.table", 
	"DT", 
	"zoo", 
	"rlang", 
	"glue", 
	"readxl", 
	"scales", 
	"lubridate", 
	"rstatix", 
	"ggplot2", 
	"ggthemes", 
	"viridis", 
	"RColorBrewer", 
	"plotly", 
	"shiny", 
	"shinyjs", 
	"shinytest", 
	"shinythemes", 
	"shinyWidgets", 
	"leaflet", 
	"tableHTML", 
	"fontawesome", 
	"rsconnect", 
	"sf", 
	"tigris"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)
