# TODO
#   * export Figures
#   * make the Table more readable
#   * maybe fix reactive binding of selected columns and PCA + Datatable
#   * Correct way to trim Sample Names
#   * fix ratio of plots(take advantage of the full page)
#   * swap to shinydashboard for sidebar 
#   * create project structure with separate files 
#   * fix issue with import PerseusR on a new machine 

# Future Ideas 
#   * Enrichments 
#   * move away from PerseusR 
#      * implement own Parser
#   * general analysis platform for DIANN and Maxquant output 


installed_packages <- "librarian" %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
     install.packages("librarian")
}

librarian::shelf(tidyverse,
                 shiny, 
                 DT, 
                 shinyWidgets, 
                 pcaMethods, 
                 ggrepel, 
                 pheatmap, 
                 EnhancedVolcano,
                 #BiocManager,
                 #devtools,
                 Biobase,
                 cox-labs/PerseusR,
                 quiet=TRUE)

trimSamples <- function(x){
  y <- unlist(strsplit(x,"[.]"))
  return(y[[length(y)-1]])
}

ProtData <- function(...){
  ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose Perseus .txt File", accept = c(".csv",".txt")),
      checkboxInput("trim","Trim SampleIDs",value = FALSE),
      uiOutput("picker"),
      actionButton("filter", "Load Data"),
      width=3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Datatable",DT::dataTableOutput("contents")),
        tabPanel("PCA",sidebarLayout(sidebarPanel(uiOutput("pickerGrouping")),
                                     mainPanel(plotOutput("pcaPlot")))),
        tabPanel("Heatmap",plotOutput("heatmap")),
        tabPanel("Volcano",sidebarLayout(sidebarPanel(uiOutput("picker1"),
                                                      uiOutput("picker2"),
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
                                                      actionButton("volcanoButton","Draw Volcano"),
                                                      width=3),
                                         mainPanel(
                                           plotOutput("volcano")
                                         ))))
     
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
      
      cols_split <- lapply(cols,trimSamples)
      cols_split <- unlist(cols_split)
      
      colnames(perseus$main) <- cols_split
      rownames(perseus$annotRows) <- cols_split
    }
    
    
    perseus
  })
  
  
  dataset_trim <- eventReactive(input$trim,{
    data <- dataset()
    
    cols <- colnames(data$main)
    
    cols_split <- lapply(cols,trimSamples)
    cols_split <- unlist(cols_split)
    
    colnames(dataset()$main) <- cols_split
    rownames(dataset()$annotRows) <- cols_split
    
  })
  
  # Select the wanted Columns for the PCA plot and displayed Datatable 
  output$picker <- renderUI({
    data <- dataset()
    pickerInput("pick","Choose Columns",choices = colnames(data$main),selected=colnames(data$main),multiple=TRUE)
  })
  
  # Select the Grouping for coloring the PCA plot
  output$pickerGrouping <- renderUI({
    data <- dataset()
    groupings <- colnames(data$annotRows)
    pickerInput("pickGrouping", "Choose Grouping",choices=groupings,selected = groupings[1],multiple = FALSE)
  })
  
  # Manually select the correct columns with the  pvalues
  output$picker1 <- renderUI({
    data <- dataset()
    
    columns <- lapply(colnames(data$annotCols),function(x) {grepl("p.value",x)})
    
    pvalueColumn <- NULL
    if (sum(columns == TRUE) > 0){
      pvalueColumn <- colnames(data$annotCols)[columns == TRUE][1]
    }
    
    pickerInput("pick1","Choose pvalue",choices=colnames(data$annotCols),selected = pvalueColumn,multiple = FALSE)
  })
  
  # Manually select the correct columns with the log2FC values
  output$picker2 <- renderUI({
    data <- dataset()
    
    
    columns <- lapply(colnames(data$annotCols),function(x) {grepl("Difference",x)})
    
    lfcColumn <- NULL
    if (sum(columns == TRUE) > 0){
      lfcColumn <- colnames(data$annotCols)[columns == TRUE][1]
    }
    
    pickerInput("pick2","Choose log2FC",choices=colnames(data$annotCols),selected=lfcColumn,multiple = FALSE)
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
  dataset_filtered <- eventReactive(input$filter,{
    data <- dataset()
    filter_idx <- colnames(data$main)  %in% input$pick
    
    return(filter_idx)
    
  })
  
  # Render the Table from the filtered Datatable
  output$contents <- DT::renderDataTable({
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
                    columnDefs = list(list( targets = 2, width = '200px')),
                  scrollX = TRUE
                   ))
    
    DT::formatRound(table,
                    columns=colnames(temp_df[-1]),
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
  
  output$heatmap <- renderPlot({
    data <- dataset()
    data_filtered <- dataset_filtered()
    
    x <- as.matrix(data$main[,data_filtered])
    annot_filtered <- data$annotRows[data_filtered,,drop=FALSE]
    
    pheatmap(x,show_rownames = FALSE,cluster_rows = TRUE,treeheight_row = 0,annotation_col = annot_filtered)
  })
  
  session$onSessionEnded(function() { stopApp() })
  }
  # Run the application 
  #shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
  shinyApp(ui = ui, server = server)
}

ProtData()