library(shiny)
library(shinyWidgets)


# Load the ui.R and server.R scripts
source("ui.R")
source("server.R")

# Create and run the Shiny app
shinyApp(ui = ui, server = server)