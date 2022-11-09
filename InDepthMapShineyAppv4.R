#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# load packages and main data objects
library(shiny)
library(shinycssloaders)
library(shinydashboard)

# library for loading the big matrices
library(gtable)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(igraph)
library(networkD3)
library(ggraph)
library(htmlwidgets)
library(gprofiler2)
library(plotly)
library(Rtsne)
library(dbscan)
library(viridis)

# create 3 matrices, one for positive, negative and combined correlation

# open matrices
#CoDepCorrMatrix <- readRDS("data/TriangleMatrix.rds")
#CoDepCorrMatrix <- readRDS("data/CoDepMapCorMatrixFastCorr.rds")
CoDepCorrMatrix <- readRDS("data/TriMat2.rds")

# import dataframe with overal cancer gene dependency scores
DependencyScoreDF <- readRDS("data/DependencyScoreDF.rds")

# import functions
source("InDepthMapHelperv5.R")

inputSearch <- colnames(CoDepCorrMatrix)
#?selectizeInput
# Define UI for application that draws a histogram
ui <- navbarPage("InDepthMap",
                
    tabPanel("Introduction",
             column(4,
               h2("Welcome to InDepthMap!"),
               h6("InDepthMap is an R shiny application for exploring cancer cell line gene co-dependencies,
                  based on data provided by DepMap."),
               h6("Specifically, InDepthMap enables users to cluster genes based on their inter-related co-dependencies,
                  revealing protein interactors, metabolic regulators and pathways important to cancer cell line proliferation."),
               
               h2("Data"),
               h6("DepMap crispr knock out data was extracted using the depmap R-package from bioconductor. Subsequently, the co-dependencies were calculated by pearson correlation."),
               h6("Currently, the resulting pearson correlation matrix is not available from this site, but will be available for download soon."),
               h6("Cancer cell gene dependencies were computed differently as compared to DepMap. For simplicity, the gene dependency score was defined as the percentage of cell lines where a knock out resulted in a Chronos score below -0.5."),
               
               
               h2("Known bugs"),
               h6("Some genes will return an error, when searched. If multiple genes are searched, one of such genes will block the entire network from being displayed."),
               
               h2("References"),
               h6("R libraries: shiny, shinycssloaders, shinydashboard, gtable, tidyr, dplyr, stringr, data.table, igraph, networkd3, ggraph, htmlwidgets, gprofiler2, plotly, Rtsne, dbscan and viridis"),
               h6(a("DepMap", href = "https://depmap.org/portal/")),
               h6(a("ShinyDepMap", href = "https://labsyspharm.shinyapps.io/depmap/"))
                )
             ),
  
    tabPanel("Analysis",

      # Sidebar with a slider input for number of bins 
      sidebarLayout(
          
          # create widgets for determining the number of top codependencis
          # Negative of positively correlated dependencies
          # Lastly a method of loading in strings for gene names
          
          
          sidebarPanel(
            
            # side panel width, 12 is max
            width = 2,
            
            selectizeInput(inputId = 'Genes',
                           label = "Genes",
                           choices = NULL,
                           multiple = TRUE,
                           width = "250px"),
            
            # reactive action button for initiating search
            actionButton("SearchGenes", "Search"),
            
            # input for the number of top codependencies to choose from
            numericInput("nTopCo", 
                         h5("Top Codepencenies n"), 
                         value = 10,
                         min = 2,
                         max = 25),
            
            # input for the correlation threshold to display
            numericInput("CorrThreshold", 
                         h5("Correlation threshold for display"), 
                         value = 0.25,
                         min = 0.01,
                         max = 1.0),
            
            # select box for correlation type
            selectInput("CorrType",
                        h5("Correlation type"),
                        choices = list("Positive",
                                       "Negative",
                                       "Both"),
                        selected = "Both"),
            
            # decide tSNE-DBSCAN input
            checkboxInput("ClustMethod",
                          h5("Use  t-SNE and DBSCAN clustering?"),
                          value = FALSE),
            
            # input for the eps variable of DBSCAN
            sliderInput("eps", h5("DBSCAN eps:"),
                        min = 0, max = 300, value = 50)
          
            
            # side panel end
            ), 
          
  
          
          
  
          # Show a plot of the generated distribution
          mainPanel(
            tabsetPanel(
                tabPanel("Network",
                      
                         width = 12,
                      
                         forceNetworkOutput("network"),
                          
                         imageOutput("depscale",
                                     width = "400px",
                                     height = "200px")
                       ),
                
                tabPanel("t-SNE + DBSCAN clustering",
                         width = 12,
                         
                         # make a tiny tSNE clustering plot
                         br(),
                         plotOutput("tSNE")
                         
                ),
                
                tabPanel("Enrichment", 
                        width = 12,
                        
                        # print a table displaying the available clusters from the network diagram
                        br(),
                        dataTableOutput("ClusterMemberTable"),
                        
                        br(),
                        h5("Below, select clusters for enrichment."),
                        h5("When the All option is chosen, the search defaults to All."),
                        selectizeInput(inputId = 'AvClusters',
                                       label = "Choose clusters to enrich",
                                       choices = NULL,
                                       multiple = TRUE,
                                       width = "250px"),
                        
                        # reactive action button for initiating search
                        actionButton("GO", "Enrich"),
                        
                        br(),
                        # generate manhattan plot with spinner
                        plotlyOutput("ManHatP") %>% withSpinner(color="#0dc5c1")
                       ),
                
                tabPanel("EnrichmentTable",
                          width = 12,
                          dataTableOutput("EnrichTable") %>% withSpinner(color="#0dc5c1")
                )
            )
        )
    )  
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ######################### citations for mainpage ##########################################

  
    ############################# reactive search input ########################################
  
    # create a reactive function for loading CoDepMatrix entries.
   
    ##### reactive search entity #####
    updateSelectizeInput(session,
                         'Genes',
                         choices = inputSearch,
                         server = TRUE,
                         selected = c("BRCA1", "USP7", "USP39"))
      
  
    # gather Genes when button is pressed, output the listed genes
     GeneList <- eventReactive(input$SearchGenes, {
       input$Genes
     })
     
    
     ############################ plot network diagram #######################################
     
     # create network data when a new genelist is generated
     data <- reactive({
       NetworkMapper(GeneList = GeneList(), 
                     Matrix = CoDepCorrMatrix,
                     DependencyScoreDF = DependencyScoreDF,
                     CorrType = input$CorrType, 
                     nCoDeps = input$nTopCo, 
                     CorrThreshold = input$CorrThreshold, 
                     ClusteringMethod = input$ClustMethod,
                     eps = input$eps,
                     NodeBorders = TRUE)
     })
     
     # now use data to output network diagram
     output$network <- renderForceNetwork({
       data()[[1]]
     })
     
     # render an external scale for dependency borders
     output$depscale <- renderImage({
       n = 100
       
       outfile <- tempfile(fileext = ".png")
       png(outfile, width=500, height=150)
       img <- image(
         x = 1:n,
         y = 1,
         z = as.matrix(1:n),
         col = viridis(n, option = "magma"),
         xlab = "Cancer cell line dependency(%)", 
         ylab = "" ,
         #ylim = 1,
         #xaxt = "n",
         yaxt = "n",
         bty = "n"
       )
       dev.off()
       list(src = outfile,
            alt = "This is alternate text")
       #return(img)
     },
     deleteFile = TRUE)
     
     ######################### plot tiny t-SNE DBSCAN output ######################
     
     # should also be dependent on slider input
     
     
     # tSNE <- Reactive(
     #   c(input$eps, input$ClustMethod, 1), {
     #   
     # })
     
     output$tSNE <- renderPlot(
       data()[[3]]
     )
     
     
     
           
     ######################### Enrichment output ##################################
     
     # before we run any enrichment, lets give the user a summary of the clusters
     # whenever data() is returned, the below output will be available
     output$ClusterMemberTable <- renderDataTable({
       CreateClusterDF(data()[[2]])
     })
    
     # create a reactive search entity for available clusters for enrichment
     observe({
       ClusterList <- c("All", paste0("C", unique(data()[[2]])))
       
       updateSelectizeInput(session,
                            'AvClusters',
                            choices = ClusterList,
                            server = TRUE,
                            selected = "All")
     })
     

     
     # Run Go Term enrichment when button is pressed, based on apparent clustering and participants of network plot
     # cluster list is element 2 of data
     GoList <- eventReactive(input$GO, {
       
       FullClusterList <- data()[[2]]
       
       # selected clusters
       if ("All" %in% input$AvClusters){
         
         return(FullClusterList)
         
       } else {
         
         # strip all text from input
         reformattedClusters <- input$AvClusters %>%
                                  
                                # extract the numbers of each list element and convert to numeric
                                lapply(function(x) stringr::str_extract(x, "\\(?[0-9,.]+\\)?")) %>%
                                as.numeric()
         
         outputList <- FullClusterList[FullClusterList %in% reformattedClusters]
         
         return(outputList)
         
       }
     })
     
     
     # dependent on GoList. When a new list is created, generate GOST enrichment object
     EnrichmentData <- reactive({
       EnrichResult <- GoFunction(GoList())
     })
     
      
    # use enrichment data to output Manhattan plots 
    output$ManHatP <- renderPlotly({
        ManHatP <- gostplot(EnrichmentData()[[1]]) %>% plotly::layout(autosize = TRUE,
                                                                 height = 200 * EnrichmentData()[[2]])
        
        return(ManHatP)
      })
    
    ######################################### table output ##############################################
    
    output$EnrichTable <- renderDataTable({
        
        EnrichTable <- publish_gosttable(EnrichmentData()[[1]], ggplot = FALSE, dataframe = TRUE)
        
        EnrichTable <- base::as.data.frame(EnrichTable)

        return(EnrichTable)
        })

}

# Run the application 
shinyApp(ui = ui, server = server)
