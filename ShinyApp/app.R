# Code to create R shiny website for group B. 
# Code split into three key sections: Preparation, UI, and Server. 
# Preparation prepares necessary objects for different figures and surrounding information. This mostly involves reading in many files
# UI creates the overall layout of the website and sets out a lot of the user interaction features
# Server puts this all together to produce live updating graphs based on user input
# Last edited: 24/04/26. Tom

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(Seurat)
library(patchwork)
library(stringr)
library(dplyr)
library(viridis)
library(DT)

# Needed for python to R conversion (For scanpy figure)
library(anndata)
library(reticulate)
library(Matrix)
py_require("scanpy")
py_require("matplotlib")
sc <- import("scanpy")
plt <- import("matplotlib.pyplot")






### PREPARATION ###

# Colour palettes
colour_palette <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo')
# Some viridis colours dont workin matplotlib, therefore different colour scheme created for python graph
py_colour_palette <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'turbo')

### Preparation for original figure 2A
original2A_DF <- read.csv("./files/Original2A.csv", header = TRUE, row.names = 1)

### Preparation for original figure 2B (Zones)
original2B_Zones_DF <- read.csv("./files/Original2B_Zones.csv", header = TRUE, row.names = 1)

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

### Preparation for alternate figure 3E
# This requires the use of anndata.
adata = read_h5ad("./ScanPy/output.h5ad")

markers_present <- c("FOXG1", "NFIA", "NFIB", "NEUROD6",
                     "BCL11A", "DCX", "OTX2", "FABP7",
                     "BCAT1", "GAD1", "DLX6", "WNT2B",
                     "RSPO2", "RSPO3", "WLS", "COL3A1",
                     "LUM", "DCN", "SPARC")

### Preparation for alternate figure 3F
# Importing data and renaming to work just like the original pipeline for simplicity 
df_alt3F <- read.csv("export_violin_data.csv", row.names = 1)
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
                              choices = c(genes_fig3E), 
                              multiple = FALSE))
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
                              choices = py_colour_palette, selected = colour_palette[4]
                            )),
                            
                            # Plot figure
                            column(5, plotOutput("altfig_3E")),
                            
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
                              "alt_3E_max", "Max Dot Size",
                              min = 0, max = 1,
                              value = 1
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
               paste("Credit group members, authors, and github")
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
    Original_tSNEPlotlyObj <- plot_ly(originalpipeline_tsne,
                                      x = originalpipeline_tsne$tSNE1,
                                      y = originalpipeline_tsne$tSNE2,
                                      
                                      # Gets user input for colour palette, and changes colours based on grouping selected
                                      color = get(input$originalpipeline_colour, originalpipeline_tsne),
                                      colors = viridis_pal(option = input$og_3D_palette)(7),
                                      
                                      type = "scatter",
                                      mode = "markers",
                                      
                                      # User input of data point opacity 
                                      alpha = input$og_3D_alpha,
                                      
                                      # Changes mouse hover text
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
  
  

  ### Original pipeline figure 3E graph
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
      
      # User input of colour palette
      cols = viridis_pal(option = input$og_3E_palette)(11),
    )
    
    # Make selected gene name the title of the plot
    OriginalPlotTwo <- OriginalPlotTwo + ggplot2::ggtitle(input$original3E_featureSelect)
  })
  
  
  
  ### Original pipeline figure 3F graph
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
  output$altfig_3E <- renderPlot({
    sc$pl$dotplot(
      adata,
      markers_present,
      groupby="cluster",
      
      # Interactive colour selection
      cmap = input$alt_3E_palette,
      
      # Interactive minimum dot size selection
      dot_min = input$alt_3E_min,
      
      # Interactive maximum dot size selection
      dot_max = input$alt_3E_max,
      
      title="Marker gene expression per cluster"
    )

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
