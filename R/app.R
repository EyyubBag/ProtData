# TODO
#   * export Figures
#   * make the Table more readable
#   * maybe fix reactive binding of selected columns and PCA + Datatable
#   * Correct way to trim Sample Names
#   * fix ratio of plots(take advantage of the full page)
#   * Feedback
#   * more efficient imports 

# Future Ideas 
#   * Enrichments 
#   * move away from PerseusR 
#      * implement own Parser -> might not need if we go forward with not using Perseus anymore 
#   * general analysis platform for DIANN and Maxquant output 



# installed_packages <- "librarian" %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#      install.packages("librarian")
# }
# 
# librarian::shelf(tidyverse,
#                  shiny,
#                  shinyWidgets,
#                  shinydashboard,
#                  DT,
#                  pcaMethods,
#                  ggpubr,
#                  ggrepel, 
#                  pheatmap, 
#                  EnhancedVolcano,
#                  #BiocManager,
#                  #devtools,
#                  Biobase,
#                  cox-labs/PerseusR,
#                  quiet=TRUE)


# installed_packages <- "renv" %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages("renv")
# }
# 
# renv::restore()

library(magrittr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(pcaMethods)
library(DT)
library(ggpubr)
library(EnhancedVolcano)
library(PerseusR)
library(pheatmap)

trimSamples <- function(x){
  y <- unlist(strsplit(x,"[.]"))
  return(y[[length(y)-1]])
}

#' @import DT
#' @import shinydashboard
#' @import shiny
#' @import shinyWidgets
#' @import stringr
#' @import ggplot2
#' @import dplyr
#' @import magrittr
#' @import Biobase
#' @importFrom ggrepel geom_text_repel
#' @importFrom pcaMethods pca scores
#' @importFrom ggpubr ggboxplot
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom PerseusR read.perseus.default
#' @importFrom pheatmap pheatmap
#' @importFrom utils stack 
#' @export 
ProtData <- function(...){

  
  ui <- dashboardPage(
    header = dashboardHeader(title = "ProtData", titleWidth = 450),
    sidebar = dashboardSidebar(
      sidebarMenu(
        fileInput("file1", "Choose Perseus .txt File", accept = c(".csv",".txt")),
        checkboxInput("trim","Trim SampleIDs",value = FALSE),
        uiOutput("picker"),
        actionButton("filter", "Load Data"),
        uiOutput("pcaTab"),
        uiOutput("volcanoTab")
      )
    ),
    body = dashboardBody(
      tabsetPanel(id="tabs",
        tabPanel("Datatable",DT::DTOutput("contents")),
        tabPanel("Summary",DT::DTOutput("summary")),
        tabPanel("Boxplots",plotOutput("boxplots")),
        tabPanel("PCA",plotOutput("pcaPlot")),
        tabPanel("Heatmap",plotOutput("heatmap")),
        tabPanel("Volcano",plotOutput("volcano"))
      )
    )
  )

server <- function(input, output, session) {
  
  # Load in Perseus Output File 
  dataset <- reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext  %in%  c("csv","txt"), "Please upload a Perseus file"))
    
    # Perseus seems to produce files with mixed decimal seperator values
    # therefore we replace all commas with dots 
    data_raw <- readLines(file$datapath)
    data_corrected <- str_replace_all(data_raw,",",".")
    
    
    perseus <- PerseusR::read.perseus.default(textConnection(data_corrected))
    
    if (input$trim){
      cols <- colnames(perseus$main)
      if (grepl(".",cols[1],fixed=TRUE)){
      cols_split <- lapply(cols,trimSamples)
      cols_split <- unlist(cols_split)
      
      colnames(perseus$main) <- cols_split
      rownames(perseus$annotRows) <- cols_split
      }
    }
    
    
    perseus
  })
  
  # Select the wanted Columns for the PCA plot and displayed Datatable 
  output$picker <- renderUI({
    data <- dataset()
    pickerInput("pick","Choose Columns",choices = colnames(data$main),selected=colnames(data$main),multiple=TRUE)
  })
  
  
  # render the appropriate UI elements for the PCA and Boxplot tabs
  # only shows them when the PCA or Boxplot tab are the current active tabs
  output$pcaTab <- renderUI({
    data <- dataset()
    groupings <- colnames(data$annotRows)
    
    tagList(
      conditionalPanel(condition='input.tabs == "PCA" || input.tabs == "Boxplots"',
                       pickerInput("pickGrouping", "Choose Grouping",choices=groupings,selected = groupings[1],multiple = FALSE)
                       )
    )
  })
  
  # render the appropriate UI elements for the Volcano tab
  # only shows them when the VolcanoTab is the current active tab
  output$volcanoTab <- renderUI({
    data <- dataset()
    
    pvalueColumns <- columns <- lapply(colnames(data$annotCols),function(x) {grepl("p.value",x)})
    
    pvalueColumn <- NULL
    if (sum(pvalueColumns == TRUE) > 0){
      pvalueColumn <- colnames(data$annotCols)[pvalueColumns == TRUE][1]
    }
    
    
    FCColumns <- lapply(colnames(data$annotCols),function(x) {grepl("Difference",x)})
    
    lfcColumn <- NULL
    if (sum(FCColumns == TRUE) > 0){
      lfcColumn <- colnames(data$annotCols)[FCColumns == TRUE][1]
    }
    
    tagList(
      conditionalPanel(condition = 'input.tabs == "Volcano"',
                       pickerInput("pick1","Choose pvalue",choices=colnames(data$annotCols),selected = pvalueColumn,multiple = FALSE),
                       pickerInput("pick2","Choose log2FC",choices=colnames(data$annotCols),selected=lfcColumn,multiple = FALSE),
                       numericInput("fcCutoff","+/- FC cutoff",
                                    value = 2,
                                    min=0,
                                    max=10,
                                    step=0.5),
                       numericInput("pValueCutoff","pvalue cutoff",
                                    value = 0.05,
                                    min=0,
                                    max=0.1,
                                    step=0.01),
                       sliderInput("overlaps",
                                   "Number of Overlaps:",
                                   min = 1,
                                   max = 50,
                                   value = 30),
                       actionButton("volcanoButton","Draw Volcano")
                       )
    )
  })
  
  


  # prepare a custom Dataframe volcano_df to draw the VolcanoPlot 
  volcano_df <- eventReactive(input$volcanoButton,{
    data <- dataset()
    df <- data$annotCols  %>% select(c(input$pick1,input$pick2,"Genes")) 
    colnames(df) = c("pValue","FC","Genes")
    
    #generally Perseus outputs the pValues in -log10(p), which we need to transform back to normal pvalues
    #for the Volcano
    df$FC <- as.numeric(df$FC)
    df$pValue <- 10**-(as.numeric(df$pValue))
    
    return(df)
    })
  
  # Render the VolcanoPlot from the volcano_df
  output$volcano <- renderPlot({
    
    df <- volcano_df()
    
    fcCutoff <- input$fcCutoff
    pvalueCutoff <- input$pValueCutoff
    
    EnhancedVolcano(df,
                    lab = df$Genes,
                    x = 'FC',
                    y = 'pValue',
                    pCutoff = pvalueCutoff,
                    FCcutoff = log2(fcCutoff),
                    # xlim = c(-10, 10),
                    # ylim = c(0, -log10(10e-12)),
                    pointSize = 1.5,
                    labSize = 2.5,
                    title = 'Volcano',
                    subtitle = 'Differential expression',
                    caption = sprintf('FC cutoff, %.2f; p-value cutoff, %.2f',fcCutoff,pvalueCutoff),
                    legendPosition = "bottom",
                    legendLabSize = 14,
                    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                    colAlpha = 0.9,
                    max.overlaps = input$overlaps,
                    drawConnectors = TRUE,
                    #hline = c(10e-8),
                    widthConnectors = 0.5)
    
  })
  
  
  # filter the Dataset according to the selected Columns
  # return the index of the kept columns
  dataset_filtered <- eventReactive(input$filter,{
    data <- dataset()
    filter_idx <- colnames(data$main)  %in% input$pick
    return(filter_idx)
    
  })
  
  
  output$summary <- DT::renderDT({
    data <- dataset()
    filter <- dataset_filtered()
    
    main_filtered <- data$main[,filter]
    annot_filtered <- data$annotRow[filter,]
    
    temp_df <- main_filtered %>% 
      reframe(across(where(is.numeric),.names="{.col}-{.fn}",
                     .fns = list(Median = median,
                                 Mean=mean,
                                 SD=sd,
                                 SE = ~sd(.)/sqrt(n()),
                                 Min = min,
                                 Max = max,
                                 q25 = ~quantile(., 0.25), 
                                 q75 = ~quantile(., 0.75),
                                 Valid = ~sum(!is.na(.)),
                                 "Valid %" = ~sum(!is.na(.))/length(.) * 100
                                 ))) %>% 
      pivot_longer(everything(), names_sep = "-", names_to = c( "Samples", ".value")) %>% 
      left_join(annot_filtered %>% mutate(Samples = row.names(annot_filtered)),by="Samples")
    
    
    table <- DT::datatable(temp_df,
                           escape=FALSE,
                           options = list(
                             pageLength = 20, 
                             autoWidth = TRUE,
                             columnDefs = list(list( targets = 2, width = '200px')),
                             scrollX = TRUE
                           ))
    
    DT::formatRound(table,
                    columns=colnames(temp_df[unlist(lapply(temp_df,is.double))]),
                    digits=2)
    
  })
  
  # Render the Table from the filtered Datatable
  output$contents <- DT::renderDT({
    data <- dataset()
    filter <- dataset_filtered()
    
    temp_df <- data$main[,filter]
    if ("Genes"  %in%  colnames(data$annotCols)){
      temp_df <- cbind(data$annotCols$Genes,data$main[,filter])
    }
    
    colnames(temp_df)[1] <- "Genes"
    
    table <- DT::datatable(temp_df,
                  escape=FALSE,
                  options = list(
                    pageLength = 20, 
                    autoWidth = TRUE,
                    columnDefs = list(list( targets = 2, width = '100px')),
                  scrollX = TRUE
                   ))
    
    DT::formatRound(table,
                    columns=colnames(temp_df[unlist(lapply(temp_df,is.double))]),
                    digits=3)
    })

  
  # PCA plot of the selected columns, colored according to the selected grouping   
  output$pcaPlot <- renderPlot({
    data <- dataset()
    filter <- dataset_filtered()
    
    data$main[,filter] %>% 
      t() %>% 
      pcaMethods::pca() %>% 
      scores() %>% 
      as.data.frame() %>% 
      ggplot(aes(PC1,PC2,label=colnames(data$main[,filter]),col=data$annotRows[input$pickGrouping][filter,])) + 
      geom_point() + 
      geom_text_repel(show.legend = FALSE) +
      labs(col=input$pickGrouping)
    
    
  })
  
  # Boxplots of the selected columns
  # colored according to the selected Grouping
  output$boxplots <- renderPlot({
    data <- dataset()
    filter_list <- dataset_filtered()
    
    data_filtered <- data$main[,filter_list]
    
    annot_df <- as.data.frame(data$annotRows)
    annot_df$ind <- row.names(annot_df)
    annot_df_filtered <- annot_df[filter_list,]

    boxplot_df <- right_join(stack(data_filtered),annot_df_filtered,by="ind")
    
    ggboxplot(boxplot_df, x="ind",y="values", col=input$pickGrouping,
              xlab="Samples",ylab="Expression") + 
      scale_x_discrete(guide = guide_axis(angle=45))
    
  })
  
  # Heatmap of the selected columns
  # all annotations are automatically drawn as colorbars 
  output$heatmap <- renderPlot({
    data <- dataset()
    data_filtered <- dataset_filtered()
    
    x <- as.matrix(data$main[,data_filtered])
    annot_filtered <- data$annotRows[data_filtered,,drop=FALSE]
    
    pheatmap(x,show_rownames = FALSE,cluster_rows = TRUE,treeheight_row = 0,annotation_col = annot_filtered)
  })
  
  #this makes sure the program exits properly when pressing the "x" button 
  session$onSessionEnded(function() { stopApp() })
  }
  # Run the application 
  #shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
  shinyApp(ui = ui, server = server)
}
ProtData()