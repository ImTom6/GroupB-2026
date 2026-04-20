# Last edited: 20/04/26. Tom

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(Seurat)
library(patchwork)
library(stringr)

### Preparation for original figure 3D
originalpipeline_tsne <- read.csv("./files/OriginalPipeline_tSNE.csv")


### Preparation for original figure 3E
# Load seurat object for figureplot
seu <- readRDS("./files/OriginalPipeline.rds")

# Load correct gene names and corresponding ids
originalDF_3E <- read.csv("./files/Original3E_GeneNames", header = TRUE, row.names = 1)

# Set gene names to show user in ui menu
genes_fig3E <- c("FOXG1", "OTX2", "RSPO2", "DCN", "ASPM", "LIN28A", "MYT1L", "NEUROD6")


### Preparation for original figure 3F
# Load in csvs into dataframes and save gene names to be shown to user
FOXG1_3F <- read.csv("./files/FOXG1_3F.csv", header = TRUE, row.names = 1)
NEUROD6_3F <- read.csv("./files/NEUROD6_3F.csv", header = TRUE, row.names = 1)
OTX2_3F <- read.csv("./files/OTX2_3F.csv", header = TRUE, row.names = 1)

# Set gene names to show user in ui menu
genes_fig3F <- c("FOXG1", "NEUROD6", "OTX2")


### Preparation for alternative figure 3D
# Load alternate scanpy pipeline output to plot alt tsne graph
altpipeline_tsne <- read.csv("./files/export_tsne.csv")

# Descriptions
altDescription1 <- "This is where the description of alternate pipeline graph 1 will go. This is where the description of alternate pipeline graph 1 will go. This is where the description of alternate pipeline graph 1 will go. "

### ui tagList gives a nice page layout
ui <- tagList(
  navbarPage(
    
    # Website title
    "Group B",
    
    # Original pipeline tab
    tabPanel("Original Pipeline",
             fillPage(
               tabsetPanel(
                 
                 # Original pipeline graph one goes here
                 tabPanel("Figure 3D",
                          sidebarLayout(
                            sidebarPanel(
                              
                              # Allow the user to colour the datapoints in different ways
                              selectInput("originalpipeline_colour", "Colour by:",
                                          choices = c("Cluster" = "cluster_code",
                                                      "Experiment" = "shape_group",
                                                      "Stage" = "day")),
                            ),
                            
                              #selectizeInput(
                                #"cellsID",
                                #"Cell ID",
                                #choices = originalpipeline_tsnec$cell,
                                #multiple = TRUE,
                              #),
                            mainPanel(
                              plotlyOutput("originalfig_3D"),
                            )
                          )
                 ),
                 
                 # Original pipeline graph two goes here
                 tabPanel("Figure 3E",
                          sidebarLayout(
                            sidebarPanel(
                              selectizeInput( 
                                "original3E_featureSelect", 
                                "Select Gene:", 
                                choices = c(genes_fig3E), 
                                multiple = FALSE
                              ), 
                            ),
                          
                            mainPanel(
                              plotlyOutput("originalfig_3E"),
                          )
                        )
                 ),
                 
                 # Original pipeline graph three goes here
                 tabPanel("Figure 3F",
                          sidebarLayout(
                            sidebarPanel(
                              
                              # Allow the user to colour the datapoints in different ways
                              selectizeInput( 
                                "original3F_featureSelect", 
                                "Select Gene:", 
                                choices = c("FOXG1","NEUROD6","OTX2"), 
                                multiple = FALSE
                              ), 
                            ),
                            
                            mainPanel(
                              plotlyOutput("originalfig_3F"),
                            )
                          ))
               )
             )
    ),
    
    # Alternative pipeline tBab
    tabPanel("Alternative Pipeline",
             fillPage(
               tabsetPanel(
                 
                 # Alternative pipeline scanpy t-SNE clusters tab
                 tabPanel("Figure 3D",
                          sidebarLayout(
                            sidebarPanel(
                              
                              # Allow the user to colour the datapoints in different ways
                              selectInput("altpipeline_colour", "Colour by:",
                                          choices = c("Cluster" = "cluster",
                                                      "Experiment" = "experiment",
                                                      "Stage" = "stage")),
                            ),
                            mainPanel(
                              plotlyOutput("altfig_3D")
                            )
                          )
                 ),
                 
                 # Alternative pipeline graph 2
                 tabPanel("Graph 2",)
               ),
             )
             
             
    )
  ),
)

# Issues with consistently having this command work, change to correctly denote location of output.rds if needed
# This has to be loaded in both ui and server to correctly work


### SERVER ###

server <- function(input, output) {
  
  # Original pipeline figure 3D graph
  output$originalfig_3D<- renderPlotly({
    Original_tSNEPlotlyObj <- plot_ly(originalpipeline_tsne,
                                      x = originalpipeline_tsne$tSNE1,
                                      y = originalpipeline_tsne$tSNE2,
                                      color = get(input$originalpipeline_colour, originalpipeline_tsne),
                                      colors = "Set3",
                                      type = "scatter",
                                      mode = "markers",
                                      text = paste("Cell: ", originalpipeline_tsne$sample, "\n","Cluster: ", originalpipeline_tsne$cluster_code, sep = "")
                                      )
    
    # Adding cluster label annotations at centroids
    # Setting location of cluster label based on the mean of datapoints for that cluster
    originalcentroids <- originalpipeline_tsne |> group_by(cluster_code)|> summarise(x_mean = mean(tSNE1), y_mean = mean(tSNE2))
    
    # Adding the annotations onto the graph
    add_annotations(
      p= Original_tSNEPlotlyObj, 
      text = paste(originalcentroids$cluster_code),
      data = originalcentroids, 
      x = originalcentroids$x_mean, y = originalcentroids$y_mean, 
      showarrow = FALSE
    )
  })
  
  
  
  
  
  
  # Original pipeline figure 3E graph
  #Convert user gene selection into correct gene id
  features_3E <- reactive({
    switch(input$original3E_featureSelect,
           "FOXG1" = originalDF_3E$gene_ids_fig3E[1],
           "OTX2" = originalDF_3E$gene_ids_fig3E[2],
           "RSPO2" = originalDF_3E$gene_ids_fig3E[3],
           "DCN" = originalDF_3E$gene_ids_fig3E[4],
           "ASPM" = originalDF_3E$gene_ids_fig3E[5],
           "LIN28A" = originalDF_3E$gene_ids_fig3E[6],
           "MYT1L" = originalDF_3E$gene_ids_fig3E[7],
           "NEUROD6" = originalDF_3E$gene_ids_fig3E[8])
  })
  
  
  
  # Generate a feature plot for the selected gene
  output$originalfig_3E <- renderPlotly({
    OriginalPlotTwo <- FeaturePlot(
      seu,
      features = features_3E(),
      reduction = "tsne",
      dims = c(1,2),
      cols = c("lightgrey", "red"),
    )
    
    # Make selected gene name the title of the plot
    OriginalPlotTwo <- OriginalPlotTwo + ggplot2::ggtitle(input$original3E_featureSelect)

  })
  

  
  
  
  
  # Original pipeline figure 3F graph
  
  # Converts user input into correct dataframe
  fig3F_df <- reactive({
    switch(input$original3F_featureSelect,
           "FOXG1" = FOXG1_3F,
           "NEUROD6" = NEUROD6_3F,
           "OTX2" = OTX2_3F)
  })
  
  # Plot original pipeline figure 3F
  output$originalfig_3F <- renderPlotly({
    
    # Ensure one dataframe is selected
    req(fig3F_df())
    ggplot(data = fig3F_df(), aes(x = region_group, y = expr)) +
      geom_violin(
        fill = "grey70",    # violin fill color
        color = "black",    # outline color
        trim = FALSE,       # show full distribution tails
        width = 0.8,
        adjust = 1.5,       # smoothing factor
        scale = "width"     # normalize widths across groups
      ) +
      geom_jitter(
        width = 0.12,       # horizontal jitter
        height = 0,
        size = 0.3,
        alpha = 0.35,
        color = "black"
      ) +
      ggtitle(input$original3F_featureSelect) +      # use gene symbol as title
      xlab(NULL) +
      ylab("LOG2 FPKM") +
      
      # Fix y-axis range for consistent comparison across genes
      coord_cartesian(ylim = c(-1, 11), expand = FALSE) +
      # Custom y-axis ticks and labels (sparse labeling for clarity)
      scale_y_continuous(
        breaks = c(0.5, 2.5, 5, 7.5, 9.5),
        labels = c("", "0", "", "10", "")
      ) +
      # Apply clean theme and adjust styling
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
        panel.grid.minor = element_line(color = "grey94", linewidth = 0.2),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
      )
  })
  
  
  
  
  
  
  
  # Plot alternative pipline figure 3D
  output$altfig_3D <- renderPlotly({
    # Using plot_ly to render a similar graph to that of the original pipeline graph
    Alt_tSNEPlotlyObj <- plot_ly(altpipeline_tsne,
                                 x = altpipeline_tsne$tSNE1,
                                 y = altpipeline_tsne$tSNE2,
      
                                 # Changing colours based on user input selection, and using Set3 as it has enough unique colours to plot each cluster
                                 color = get(input$altpipeline_colour, altpipeline_tsne),
                                 colors = "Set3",
      
                                 # Plotting correct graph type
                                 type = "scatter",
                                 mode = 'markers',
      
                                 # Setting text that shows when the mouse hovers over a data point
                                 text = paste("Cell: ", altpipeline_tsne$cell, "\n","Cluster: ", altpipeline_tsne$cluster, sep = "")
                                 ) 
    
    # Adding cluster label annotations at centroids
    # Setting location of cluster label based on the mean of datapoints for that cluster
    centroids <- altpipeline_tsne |> group_by(cluster)|> summarise(x_mean = mean(tSNE1), y_mean = mean(tSNE2))
    
    # Adding the annotations onto the graph
    add_annotations(
      p= Alt_tSNEPlotlyObj, 
      text = centroids$cluster, 
      data = centroids, 
      x = centroids$x_mean, y = centroids$y_mean, 
      showarrow = FALSE
    )
  })
  
  
  
  
}

shinyApp(ui = ui, server = server)
