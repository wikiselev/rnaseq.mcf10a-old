library(shiny)
library(ggplot2)
 
shinyServer(function(input, output) {

  tmp <- readRDS("data/plot-time-courses-all-genes.rds")
  dataset <- reactive({
    tmp[tmp$id == input$geneName,]
  })
 
  output$plot <- renderPlot({

    limits <- aes(ymax = ymax, ymin = ymin)

    p <- ggplot(dataset(),
      aes(time, value, group = cond, color = cond)) +
        geom_line() + facet_wrap(~ id, scale = "free_y") +
        geom_errorbar(limits, width = 0.25) +
        theme_bw()
    
    # p <- ggplot(dataset(), aes_string(x=input$x, y=input$y)) + geom_point()
    
    # if (input$color != 'None')
    #   p <- p + aes_string(color=input$color)
    
    # facets <- paste(input$facet_row, '~', input$facet_col)
    # if (facets != '. ~ .')
    #   p <- p + facet_grid(facets)
    
    # if (input$jitter)
    #   p <- p + geom_jitter()
    # if (input$smooth)
    #   p <- p + geom_smooth()
    
    print(p)
    
  }, height=400)
  
})