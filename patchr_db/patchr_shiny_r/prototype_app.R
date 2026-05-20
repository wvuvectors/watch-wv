library(shiny)
library(httr)
library(jsonlite)
library(DT)

# ----------------------------
# User Interface
# ----------------------------

ui <- fluidPage(
  titlePanel("PATCHR Sample Query Prototype"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("sample_id", "Sample ID", placeholder = "e.g. A0028"),
      actionButton("search", "Search"),
      actionButton("clear", "Clear")
    ),
    
    mainPanel(
      DTOutput("samples_table")
    )
  )
)

# ----------------------------
# Server Logic
# ----------------------------

server <- function(input, output, session) {
  
  # Reactive data container
  samples_data <- reactiveVal(data.frame())
  
  # Load first 100 rows on startup 
  observe({
    tryCatch({
      
      response <- GET("http://127.0.0.1:8000/samples/?limit=100")
      
      if (status_code(response) == 200) {
        data <- fromJSON(content(response, "text", encoding = "UTF-8"))
        samples_data(do.call(rbind, lapply(data, as.data.frame)))
      } else {
        showNotification("Failed to load initial data", type = "error")
      }
      
    }, error = function(e) {
      print(e)
      showNotification("Error connecting to API", type = "error")
    })
  })
  
  # Search button
  observeEvent(input$search, {
    
    sample_id <- trimws(input$sample_id)
    
    if (sample_id == "") {
      showNotification("Please enter a sample ID", type = "warning")
      return()
    }
    
    url <- paste0(
      "http://127.0.0.1:8000/samples/query?sample_id=",
      URLencode(sample_id)
    )
    
    response <- GET(url)
    
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, "text", encoding = "UTF-8"))
      
      if (length(data) == 0) {
        showNotification("No matching sample found", type = "message")
        samples_data(data.frame())
      } else {
        samples_data(as.data.frame(data))
      } 
      
    } else {
      showNotification("API request failed", type = "error")
    }
  })
  
  # Clear button
  observeEvent(input$clear, {
    updateTextInput(session, "sample_id", value = "")
    
    response <- GET("http://127.0.0.1:8000/samples/?limit=100")
    
    if (status_code(response) == 200) {
      data <- fromJSON(content(response, "text", encoding = "UTF-8"))
      samples_data(as.data.frame(data))
    }
  })
  
  # Render table
  output$samples_table <- renderDT({
    datatable(
      samples_data(),
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
}

# ----------------------------
# Run app
# ----------------------------
shinyApp(ui = ui, server = server)