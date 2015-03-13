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
                checkboxInput('savePlot', "Save the plot as PDF"),
                actionButton("goButton", "Plot!"),
                div("Please click on Download link below only after the plot appears on the right.", style = "color:black"),
                conditionalPanel(
                        condition = "input.savePlot == true & input.goButton != 0",
                        downloadLink('pdflink')
                )
        ),
 
        mainPanel(
                plotOutput('plot')
        )
))