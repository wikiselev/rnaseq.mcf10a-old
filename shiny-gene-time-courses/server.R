library(shiny)
library(ggplot2)

shinyServer(function(input, output) {
        output$plot <- renderPlot({
                if (input$goButton == 0)
                        return()
                data <- isolate({
                        d <- readRDS("data/plot-data.rds")
                        if(input$checkbox) {
                                d <- d[d$norm_by_glength]
                        } else {
                                d <- d[!d$norm_by_glength]
                        }
                        validate(
                                need((input$geneName %in% d$id) | (input$geneName %in% d$hgnc_symbol),
                                     "Your gene is not in our list. We only accept HGNC symbols and Ensembl IDs. Please check that your gene name is either of the above...")
                        )
                        d[(d$id %in% input$geneName) | (d$hgnc_symbol %in% input$geneName)]
                        
                })
                limits <- aes(ymax = ymax, ymin = ymin)
                p <- ggplot(data,
                        aes(time, value, group = Condition, color = Condition)) +
                        geom_line(size = 1) +
                        geom_point(size = 3) +
                        geom_errorbar(limits, size = 0.5, width = 5) +
                        labs(x = "Time, min", y = "Read counts") +
                        theme_bw()
                isolate({
                        if(input$savePlot) {
                                ggsave(paste0(input$geneName, "-rna-seq-profile.png"), p, type="cairo-png")
                        }
                })
                print(p)
        }, height=300)
})