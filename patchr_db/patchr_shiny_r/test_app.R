# test_app.R
library(shiny)

ui <- fluidPage("Hello PATCHR")
server <- function(input, output, session) {}

shinyApp(ui, server)