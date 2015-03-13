library(shiny)
library(ggplot2)
library(data.table)

shinyServer(
        function(input, output) {
                # this will load the data only once when the webpage is open
                d <- readRDS("data/plot-data.rds")
                # the function below will be executed each time a user presses goButton
                output$plot <- renderPlot({
                        if (input$goButton == 0)
                                return()
                        data <- isolate({
                                if(input$checkbox) {
                                        d <- d[d$norm_by_glength]
                                } else {
                                        d <- d[!d$norm_by_glength]
                                }
                                validate(
                                        need((toupper(input$geneName) %in% d$id) | (toupper(input$geneName) %in% d$hgnc_symbol),
                                             "Your gene is not in our list. We only accept HGNC symbols and Ensembl IDs. Please check that your gene name is either of the above...")
                                )
                                d[(d$id %in% toupper(input$geneName)) | (d$hgnc_symbol %in% toupper(input$geneName))]
                        })
                        limits <- aes(ymax = ymax, ymin = ymin)
                        p <- ggplot(data,
                                aes(time, value, group = Condition, color = Condition)) +
                                geom_line(size = 1) +
                                geom_point(size = 3) +
                                geom_errorbar(limits, size = 0.5, width = 5) +
                                facet_grid(. ~ hgnc_symbol) +
                                labs(x = "Time, min", y = "Read counts") +
                                theme_bw()
                        isolate({
                                if(input$savePlot) {
                                        pdf("plot.pdf", width = 7, height = 4)
                                        print(p)
                                        dev.off()
                                }
                        })
                        print(p)
                }, height=400, width = 700)
                output$pdflink <- downloadHandler(
                        filename <- function() {
                                paste0(input$geneName, "-rna-seq-profile.pdf")
                        },
                        content <- function(file) {
                                file.copy("plot.pdf", file)
                        }
                )
        }
)