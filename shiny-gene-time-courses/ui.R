library(shiny)
library(ggplot2)
 
# dataset <- diamonds

# dataset <- readRDS("plot-time-courses-all-genes.rds")
 
shinyUI(pageWithSidebar(
 
  headerPanel("Plot genes"),
  
  sidebarPanel(
 
    textInput('geneName', 'Gene Ensembl ID')

  ),
 
  mainPanel(
    plotOutput('plot')
  )
))