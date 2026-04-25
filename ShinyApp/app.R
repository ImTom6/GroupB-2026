# Code to create R shiny website for group B. 
# Code split into three key sections: Preparation, UI, and Server. 
# Preparation prepares necessary objects for different figures and surrounding information. This mostly involves reading in many files
# UI creates the overall layout of the website and sets out a lot of the user interaction features
# Server puts this all together to produce live updating graphs based on user input
# Last edited: 24/04/26. Tom

library(shiny)
#library(bslib)
library(ggplot2)
library(plotly)
library(Seurat)
#library(patchwork)
#library(stringr)
#library(dplyr)
library(viridis)
library(DT)


### PREPARATION ###

# Colour palettes
colour_palette <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo')

### Preparation for original figure 2A
original2A_DF <- read.csv("./files/Original2A.csv", header = TRUE, row.names = 1)

### Preparation for original figure 2B (Zones)
original2B_Zones_DF <- read.csv("./files/Original2B_Zones.csv", header = TRUE, row.names = 1)

### Preparation for original figure 3D, 3E
seurat_organoid <- readRDS("./files/Camp2015_organoid_final.rds")

### Preparation for original figure 3E
genes_3e <- c("FOXG1", "OTX2", "RSPO2", "DCN",
              "ASPM", "LIN28A", "MYT1L", "NEUROD6")

### Preparation for original figure 3F
seurat_3f <- readRDS("./files/Camp2015_fig3f.rds")

genes_3f <- c("FOXG1", "NEUROD6", "OTX2")
genes_3f <- genes_3f[genes_3f %in% rownames(seurat_3f)]

Idents(seurat_3f) <- factor(seurat_3f$region_group,
                            levels = c("r1","r2","r3","r4","fetal"))

### Preparation for alternative figure 3D
# Load alternate scanpy pipeline output to plot alt tsne graph
altpipeline_tsne <- read.csv("./files/export_tsne.csv")

### Preparation for alternate figure 3E
# This requires the use of anndata.
scpy_data <- readRDS("./files/scpy_data.rds")

# Markers identified from scanpy pipeline
markers_present <- c("FOXG1", "NFIA", "NFIB", "NEUROD6",
                     "BCL11A", "DCX", "OTX2", "FABP7",
                     "BCAT1", "GAD1", "DLX6", "WNT2B",
                     "RSPO2", "RSPO3", "WLS", "COL3A1",
                     "LUM", "DCN", "SPARC")

### Preparation for alternate figure 3F
# Importing data and renaming to work just like the original pipeline for simplicity 
df_alt3F <- read.csv("./files/export_violin_data.csv", row.names = 1)
FOXG1_Alt_3F <- select(df_alt3F, c('FOXG1','group'))
NEUROD6_Alt_3F <- select(df_alt3F, c('NEUROD6','group'))
OTX2_Alt_3F<- select(df_alt3F, c('OTX2','group'))

FOXG1_Alt_3F <- rename(FOXG1_Alt_3F, expr = FOXG1, region_group = group)
NEUROD6_Alt_3F <- rename(NEUROD6_Alt_3F, expr = NEUROD6, region_group = group)
OTX2_Alt_3F <- rename(OTX2_Alt_3F, expr = OTX2, region_group = group)

# All figure legends are saved into a txt, this can be loaded into a dataframe with each row representing a different legend
figure_legend <- read.table("./files/figure_legends.txt", sep="\n",header=F)

# Text to be placed at bottom of each figure tab to credit original paper
citation <- read.table("./files/citations.txt", sep="\n", header=F)


### UI ###

### ui tagList gives a nice page layout
ui <- tagList(
  navbarPage(
    
    # Website title
    "Group B",
    
    # Original pipeline tab
    tabPanel("Original Pipeline",
             fluidPage(
               tabsetPanel(
                 
                 # Original pipeline figure 2A graph
                 tabPanel("Figure 2A",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "og_2A_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("originalfig_2A")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[1,])),
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "og_2A_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_2A_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,]))
                          ),),
                 
                 tabPanel("Figure 2B (Zones)",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "og_2B_Zones_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("originalfig_Zones_2B")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[1,])),
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "og_2B_zone_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_2B_Zones_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,]))
                          ),),
                 
                 tabPanel("Figure 3D",
                          fluidRow(
                            
                              # Select colour palette
                              column(2, radioButtons(
                                "og_3D_palette", "Select Colour Palette",
                                choices = colour_palette, selected = colour_palette[4],
                              )),
                              
                              # Plot figure
                              column(5, plotlyOutput("originalfig_3D")),
                              
                              # Display figure legend
                              column(2, paste(figure_legend[1,])),
                            ),
                            fluidRow(
                              
                              # Colour based on different groupings
                              column(2, selectInput("originalpipeline_colour", "Colour by:",
                                                    choices = c("Cluster" = "cluster_code",
                                                                "Experiment" = "shape_group",
                                                                "Stage" = "day"))),
                              
                              # Slider for data point opacity
                              column(2, sliderInput(
                                "og_3D_alpha", "Data point opacity",
                                min = 0, max = 1,
                                value = 1
                              ))
                            ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_3D_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,])),
                          ),),
                 
                 # Original pipeline figure 3E
                 tabPanel("Figure 3E",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "og_3E_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("originalfig_3E")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[2,])),
                          ),
                          
                          fluidRow(
                            
                            # Select gene to plot figure for
                            column(3, selectizeInput( 
                              "original3E_featureSelect", 
                              "Select Gene:", 
                              choices = c(genes_3e), 
                              multiple = FALSE)),
                            
                            column(2, sliderInput(
                              "og_3E_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_3E_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,])),
                          ),),
                 
                 # Original pipeline figure 3F
                 tabPanel("Figure 3F",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "og_3F_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5,plotlyOutput("originalfig_3F")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[3,])),
                          ),
                          fluidRow(
                            
                            # Select gene to plot figure for
                            column(6, selectizeInput( 
                            "original3F_featureSelect", 
                            "Select Gene:", 
                            choices = c("FOXG1","NEUROD6","OTX2"), 
                            multiple = FALSE))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_3F_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,])),
                          ),),

               ))
             ),

    # Alternative pipeline Tab
    tabPanel("Alternative Pipeline",
             fluidPage(
               tabsetPanel(
                 
                 # Alternative pipeline figure 3D
                 tabPanel("Figure 2A",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_2A_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("altfig_2A")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[1,])),
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_2A_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("alt_2A_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,]))
                          ),),
                 
                 tabPanel("Figure 2B",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_2B_Zone_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("alt_fig_2B_zone")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[1,])),
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_2A_zone_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("alt_2A_zone_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,]))
                          ),),
                 
                 tabPanel("Figure 3D",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_3D_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4]
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("altfig_3D")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[4,])),
                          ),
                          
                          fluidRow(
                            
                            # Colour data points by different groupings
                            column(2, selectInput("altpipeline_colour", "Colour by:",
                                                  choices = c("Cluster" = "cluster",
                                                              "Experiment" = "experiment",
                                                              "Stage" = "stage"))),
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_3D_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(5, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Alt_3D_Table"))
                            
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,])),
                          ),),
                              
                 # Alternative pipeline figure 3E
                 tabPanel("Figure 3E",
                          
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_3E_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4]
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("altfig_3E")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[5,])),
                          ),
                          fluidRow(
                            
                            # Slider to change minimum dot size
                            column(2, sliderInput(
                              "alt_3E_min", "Min Dot Size",
                              min = 0, max = 1,
                              value = 0
                          )),
                          
                            # Slider to change maximum dot size
                            column(2, sliderInput(
                              "alt_3E_max", "Dot Scale",
                              min = 0, max = 10,
                              value = 4
                          ))),

                          fluidRow(
                            column(5, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Alt_3E_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,])),
                          ),),
                              
                 # Alternate pipeline figure 3F
                 tabPanel("Figure 3F",
                          
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_3F_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("alternatefig_3F")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[6,])),
                          ),
                          
                          
                          fluidRow(
                            # Select gene to generate figure of
                            column(3, selectizeInput( 
                              "alternate3F_featureSelect", 
                              "Select Gene:", 
                              choices = c("FOXG1","NEUROD6","OTX2"), 
                              multiple = FALSE))
                          ),
                          
                          fluidRow(
                            column(5, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Alt_3F_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,])),
                          ),),
                      ))),
    
    # Final tab to give credit to group members and original paper authors
    tabPanel("Credit",
             fluidPage(
               paste("Credit group members, authors, and github .This research used the ALICE High Performance Computing facility at the University of Leicester.")
             )
    )
))








### SERVER ###

server <- function(input, output) {
  
  ### Original pipeline 2A graph
  # Show interactive datatable for data plotted in original fig 2A
  output$Og_2A_Table <- renderDataTable({datatable(original2A_DF)})
  
  # Plot original 2A graph
  output$originalfig_2A <- renderPlotly({
    original_2A_obj <- plot_ly(original2A_DF, 
            x = original2A_DF$Component_1,
            y = original2A_DF$Component_2,
            
            # Colour data points based on their label
            color = original2A_DF$paper_like_label,
            
            # User input for colour palette
            colors = viridis_pal(option = input$og_2A_palette)(7),
            type = "scatter",
            mode = "markers",
            
            # User input for data point opacity
            alpha = input$og_2A_alpha
            ) |>
    layout(xaxis = list(title = "Component 1"),
           yaxis = list(title = "Component 2"))

  })
  
  ### Original pipeline 2B Zones graph
  # Show interactive datatable for data plotted in original fig 2A
  output$Og_2B_Zones_Table <- renderDataTable({datatable(original2B_Zones_DF)})
  
  # Plot original 2A graph
  output$originalfig_Zones_2B  <- renderPlotly({
    original_2B_Zones_obj <- plot_ly(original2B_Zones_DF, 
                               x = original2B_Zones_DF$Component_1,
                               y = original2B_Zones_DF$Component_2,
                               
                               # Colour data points based on their zone
                               color = original2B_Zones_DF$Zone,
                               
                               # Get user input for colour palette 
                               colors = viridis_pal(option = input$og_2B_Zones_palette)(7),
                               type = "scatter",
                               mode = "markers",
                               
                               # User input for data point opacity
                               alpha = input$og_2B_zone_alpha
    ) |>
      # Add x and y axis titles
      layout(xaxis = list(title = "Component 1"),
             yaxis = list(title = "Component 2"))
    
  })
  
  ### Original pipeline figure 3D graph
  # Show interactive datatable for data plotted in original fig 3D
  output$Og_3D_Table <- renderDataTable({datatable(originalpipeline_tsne)}) 
  
  # Plot original 3D graph
  output$originalfig_3D <- renderPlotly({
    fig3d <- DimPlot(seurat_organoid,
                     reduction = "tsne",
                     label = TRUE,
                     label.size = 3,
                     pt.size = 1.5,
                     repel = TRUE,
                     cols = viridis_pal(option = input$og_3D_palette)(11),
                     alpha = input$og_3D_alpha,
                     shape.by = "shape_group") +
      scale_shape_manual(values = c(
        "33d" = 16, "35d" = 17, "37d" = 25,
        "65d" = 15, "41d" = 18,
        "r1 53d" = 0, "r2 53d" = 5,
        "r3 58d" = 1, "r4 58d" = 2)) +
      theme_classic() +
      labs(title = "Figure 3D. Organoid cell clusters",
           x = "tSNE 1", y = "tSNE 2",
           shape = "Stage / Region")
  })
  
  

  ### Original pipeline figure 3E graph

  
  output$originalfig_3E <- renderPlotly({
  fig3e <- FeaturePlot(
    seurat_organoid,
    features = input$original3E_featureSelect,
    reduction = "tsne",
    alpha = input$og_3E_alpha,
    cols = viridis_pal(option = input$og_3E_palette)(5),
    pt.size = 1,
    order = TRUE,
    ncol = 4,
    slot = "data"
  ) &
    theme_classic() &
    theme(plot.title = element_text(face = "italic", hjust = 0.5))
  })
  
  
  
  ### Original pipeline figure 3F graph
  # Plot original pipeline figure 3F
  output$originalfig_3F <- renderPlotly({
    fig3f <- VlnPlot(seurat_3f,
                     features = input$original3F_featureSelect,
                     pt.size = 0.1,
                     cols = viridis_pal(option = input$og_3F_palette)(5),
                     layer = "data") &
      theme_classic() &
      labs(y = "log2(FPKM + 1)") &
      theme(legend.position = "none")
  
  
  })
  
  # Plot alternative pipline figure 3D
  output$Alt_3D_Table <- renderDataTable({datatable(altpipeline_tsne)}) 
  output$altfig_3D <- renderPlotly({
    # Using plot_ly to render a similar graph to that of the original pipeline graph
    Alt_tSNEPlotlyObj <- plot_ly(altpipeline_tsne,
                                 x = altpipeline_tsne$tSNE1,
                                 y = altpipeline_tsne$tSNE2,
      
                                 # Changing colours based on user input selection, and using Set3 as it has enough unique colours to plot each cluster
                                 color = get(input$altpipeline_colour, altpipeline_tsne),
                                 colors =  viridis_pal(option = input$alt_3D_palette)(11),
      
                                 # Plotting correct graph type
                                 type = "scatter",
                                 mode = 'markers',
                                 
                                 # User input changing opacity of data points
                                 alpha = input$alt_3D_alpha,
      
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
  
  # Plot alternative pipeline figure 3E
  #output$Alt_3E_Table <- renderDataTable({datatable(adata)}) 
  
  # Requires use of reticulate to run scanpy code
  output$altfig_3E <- renderPlotly({
    DotPlot(
      scpy_data,
      features = markers_present,
      group.by = "cluster",
      scale.by = "radius",
      scale = TRUE,
      dot.min = input$alt_3E_min,
      dot.scale = input$alt_3E_max
      #cols = viridis_pal(option = input$alt_3E_palette)(2)
      
    ) + RotatedAxis()

  })
  
  
  ### Plot alternative pipeline figure 3F
  
  # Converts user input into correct dataframe
  Alt_fig3F_df <- reactive({
    switch(input$alternate3F_featureSelect,
           "FOXG1" = FOXG1_Alt_3F,
           "NEUROD6" = NEUROD6_Alt_3F,
           "OTX2" = OTX2_Alt_3F)
  })
  
  # Plot alternate pipeline figure 3F
  output$Alt_3F_Table <- renderDataTable({datatable(Alt_fig3F_df())}) 
  output$alternatefig_3F <- renderPlotly({
    
    # Ensure one dataframe is selected
    req(Alt_fig3F_df())
    ggplot(data = Alt_fig3F_df(), aes(x = region_group, y = expr, fill=region_group)) +
      geom_violin(
        linewidth = 0.8, 
        trim = FALSE,
        width = 0.8,
        adjust = 1.5, 
        scale = "width") +
      
      # Adding jitter
      geom_jitter(
        width = 0.12,
        height = 0,
        size = 0.3,
        alpha = 0.35,
        color = "black"
      ) +
      
      # User input for colour palette
      scale_fill_viridis(discrete = TRUE, option = input$alt_3F_palette) +
      
      # Changes the gene selected to be the title
      ggtitle(input$alternate3F_featureSelect) +      
      xlab("Region Group") +
      ylab("LOG2 FPKM") +

      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 9)
      ) 
  })
}

shinyApp(ui = ui, server = server)
