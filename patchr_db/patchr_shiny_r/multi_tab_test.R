library(shiny)
library(DT)

# ----------------------------
# SOURCE HELPER SCRIPTS
# ----------------------------

source("R/api_helpers.R")
source("R/samples_api.R")
source("R/concentration_api.R")

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
        textInput("batch_id", "Batch ID"),
        textInput("sample_id_c", "Sample ID"),
        actionButton("search_conc", "Search"),
        actionButton("clear_conc", "Clear")
      ),
      mainPanel(DTOutput("concentration_table"))
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
  
  
  # =========================
  # LOAD INITIAL DATA
  # =========================
  observe({
    data <- get_samples(limit = 1500)
    samples_data(data)
  })
  
  observe({
    data <- get_concentration(limit = 10000)
    concentration_data(data)
  })
  
  # =========================
  # SAMPLES SEARCH
  # =========================
  observeEvent(input$search_sample, {
    data <- query_samples_api(
      sample_id = input$sample_id
    )
    
    samples_data(data)
  })
  
  observeEvent(input$clear_sample, {
    updateTextInput(session, "sample_id", value = ""
    )
    
    data <- get_samples(limit = 100)
    
    samples_data(data)
  })
  
  # =========================
  # CONCENTRATION SEARCH
  # =========================
  observeEvent(input$search_conc, {
    data <- query_concentration_api(
      concentration_id = input$concentration_id,
      concentration_batch_id = input$batch_id,
      sample_id = input$sample_id_c
    )
    
    concentration_data(data)
  })
  
  observeEvent(input$clear_conc, {
    updateTextInput(session, "concentration_id", value = ""
    )
    
    updateTextInput(session, "batch_id", value = ""
    )
    
    updateTextInput(session, "sample_id_c", value = ""
    )
    
    data <- get_concentration(limit = 100)
    
    concentration_data(data)
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
}

# ----------------------------
# RUN
# ----------------------------
shinyApp(ui, server)