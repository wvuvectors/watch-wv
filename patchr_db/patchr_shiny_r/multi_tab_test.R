library(shiny)
library(DT)

# ----------------------------
# SOURCE HELPER SCRIPTS
# ----------------------------

source("R/api_helpers.R")
source("R/samples_api.R")
source("R/concentration_api.R")
source("R/assay_api.R")

# ----------------------------
# UI
# ----------------------------
ui <- navbarPage(
  "PATCHR Dashboard",
  
  # ---------------- SAMPLES TAB ----------------
  tabPanel(
    "Samples",
    sidebarLayout(
      sidebarPanel(
        textInput("sample_id", "Sample ID"),
        textInput("sample_status", "Sample Status"),
        textInput("location_id", "Location ID"),
        textInput("sample_qc", "Sample QC"),
        dateInput("recovered_start", "Recovered Start"),
        dateInput("recovered_end", "Recovered End"),
        actionButton("search_sample", "Search"),
        actionButton("clear_sample", "Clear")
      ),
      mainPanel(DTOutput("samples_table"))
    )
  ),
  
  # ---------------- CONCENTRATION TAB ----------------
  tabPanel(
    "Concentration",
    sidebarLayout(
      sidebarPanel(
        textInput("concentration_id", "Concentration ID"),
        textInput("sample_id_c", "Sample ID"),
        textInput("batch_id", "Batch ID"),
        textInput("batch_location_c", "Batch Location"),
        actionButton("search_conc", "Search"),
        actionButton("clear_conc", "Clear")
      ),
      mainPanel(DTOutput("concentration_table"))
    )
  ), 
  
  # ---------------- EXTRACTIONS TAB ----------------
  tabPanel(
    "Extractions",
    sidebarLayout(
      sidebarPanel(
        textInput("extraction_id", "Extraction ID"),
        textInput("concentration_id", "Concentration ID"),
        textInput("extraction_batch_id", "Batch ID"),
        textInput("extraction_location_in_batch", "Batch Location"),
        textInput("extraction_location_int_storage", "Storage Location"),
        actionButton("search_extr", "Search"),
        actionButton("clear_extr", "Clear")
      ),
      mainPanel(DTOutput("extractions_table"))
    )
  ),
  
  # ---------------- ASSAYS TAB ----------------
  tabPanel(
    "Assays",
    sidebarLayout(
      sidebarPanel(
        textInput("assay_id", "Assay ID"),
        textInput("extraction_id", "Extraction ID"),
        textInput("sample_id", "Sample ID"),
        textInput("assay_batch_id", "Batch ID"),
        textInput("assay_location_in_batch", "Batch Location"),
        numericInput("assay_input_ul", "Input ul", value = NA),
        textInput("assay_class", "Assay Class"),
        textInput("assay_type", "Assay Type"),
        numericInput("assay_target_copies_per_ul_reaction", "Copies/ul", value = NA),
        actionButton("search_assay", "Search"),
        actionButton("clear_assay", "Clear")
      ),
      mainPanel(DTOutput("assay_table"))
    )
  )
)

# ----------------------------
# SERVER
# ----------------------------
server <- function(input, output, session) {
  
  # =========================
  # Reactive storage
  # =========================
  samples_data <- reactiveVal(data.frame())
  concentration_data <- reactiveVal(data.frame())
  extractions_data <- reactiveVal(data.frame())
  assay_data <- reactiveVal(data.frame())
  
  # =========================
  # LOAD INITIAL DATA
  # =========================
  observe({
    data <- get_samples(limit = as.integer(100000))
    samples_data(data)
  })
  
  observe({
    data <- get_concentration(limit = as.integer(100000))
    concentration_data(data)
  })
  
  observe({
    data <- get_extractions(limit = as.integer(100000))
    extractions_data(data)
  })
  
  observe({
    data <- get_assays(limit = as.integer(100000))
    assay_data(data)
  })
  
  # =========================
  # SAMPLES SEARCH
  # =========================
  observeEvent(input$search_sample, {
    data <- query_samples_api(
      sample_id = input$sample_id,
      sample_status = input$sample_status,
      location_id = input$location_id,
      sample_qc = input$sample_qc,
      recovered_start = input$recovered_start,
      recovered_end = input$recovered_end
    )
    
    samples_data(data)
  })
  
  observeEvent(input$clear_sample, {
    updateTextInput(session, "sample_id", value = "")
    updateTextInput(session, "sample_status", value = "")
    updateTextInput(session, "location_id", value = "")
    updateTextInput(session, "sample_qc", value = "")
    updateDateInput(session, "recovered_start", value = "")
    updateDateInput(session, "recovered_end", value = "")
    
    data <- get_samples(limit = as.integer(100000))
    
    samples_data(data)
  })
  
  # =========================
  # CONCENTRATION SEARCH
  # =========================
  observeEvent(input$search_conc, {
    data <- query_concentration_api(
      concentration_id = input$concentration_id,
      concentration_batch_id = input$batch_id,
      sample_id = input$sample_id_c,
      concentration_location_in_batch = input$batch_location_c
    )
    
    concentration_data(data)
  })
  
  observeEvent(input$clear_conc, {
    updateTextInput(session, "concentration_id", value = "")
    updateTextInput(session, "batch_id", value = "")
    updateTextInput(session, "sample_id_c", value = "")
    updateTextInput(session, "batch_location_c", value = "")
    
    data <- get_concentration(limit = as.integer(100000))
    
    concentration_data(data)
  })
  
  # =========================
  # EXTRACTIONS SEARCH
  # =========================
  observeEvent(input$search_extr, {
    data <- query_extractions_api(
      extraction_id = input$extraction_id,
      concentration_id = input$concentration_id_e,
      extraction_batch_id = input$extraction_batch_id,
      extraction_location_in_batch = input$extraction_location_in_batch,
      extraction_location_in_storage = input$extraction_location_in_storage
    )
    
    extractions_data(data)
  })
  
  observeEvent(input$clear_extr, {
    updateTextInput(session, "extraction_id", value = "")
    updateTextInput(session, "concentration_id_e", value = "")
    updateTextInput(session, "extraction_batch_id", value = "")
    updateTextInput(session, "extraction_location_in_batch", value = "")
    updateTextInput(session, "extraction_location_in_storage", value = "")
    
    data <- get_extractions(limit = as.integer(100000))
    
    extractions_data(data)
  })
  
  # =========================
  # ASSAYS SEARCH
  # =========================
  observeEvent(input$search_assay, {
    data <- query_assays_api(
      assay_id = input$assay_id,
      extraction_id = input$extraction_id,
      sample_id = input$sample_id,
      assay_batch_id = input$assay_batch_id,
      assay_location_in_batch = input$assay_location_in_batch,
      assay_input_ul = input$assay_input_ul,
      assay_class = input$assay_class,
      assay_type = input$assay_type,
      assay_target_copies_per_ul_reaction = input$assay_target_copies_per_ul_reaction
    )
    
    assay_data(data)
  })
  
  observeEvent(input$clear_assay, {
    updateTextInput(session, "assay_id", value = "")
    updateTextInput(session, "extraction_id", value = "")
    updateTextInput(session, "sample_id", value = "")
    updateTextInput(session, "assay_batch_id", value = "")
    updateTextInput(session, "assay_location_in_batch", value = "")
    updateNumericInput(session, "assay_input_ul", value = NULL)
    updateTextInput(session, "assay_class", value = "")
    updateTextInput(session, "assay_type", value = "")
    updateNumericInput(session, "assay_target_copies_per_ul_reaction", value = NULL)
    
    data <- get_assays(limit = as.integer(100000))
    
    assay_data(data)
  })
  
  # =========================
  # RENDER TABLES
  # =========================
  output$samples_table <- renderDT({
    datatable(samples_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$concentration_table <- renderDT({
    datatable(concentration_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$extractions_table <- renderDT({
    datatable(extractions_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$assay_table <- renderDT({
    datatable(assay_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
}

# ----------------------------
# RUN
# ----------------------------
shinyApp(ui, server)