library(shiny)
library(ggplot2)
 
# dataset <- diamonds

# dataset <- readRDS("plot-time-courses-all-genes.rds")
 
shinyUI(pageWithSidebar(
        headerPanel(
                "Plot RNA-Seq profiles for a selected gene"
        ),
        sidebarPanel(
                textInput('geneName', 'Gene HGNC symbol / Ensembl ID'),
                checkboxInput("checkbox", label = "Normalize by the gene length", value = FALSE),
                checkboxInput('savePlot', "Check to save the plot as PDF"),
                actionButton("goButton", "Plot!"),
                div("After pressing this button wait until the data is loaded", style = "color:blue")
        ),
 
        mainPanel(
                plotOutput('plot')
        )
))