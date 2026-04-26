# Code to create R shiny website for group B. 
# Code split into three key sections: Preparation, UI, and Server. 
# Preparation prepares necessary objects for different figures and surrounding information. This mostly involves reading in many files
# UI creates the overall layout of the website and sets out a lot of the user interaction features
# Server puts this all together to produce live updating graphs based on user input
# Last edited: 26/04/26. Tom

library(shiny)
library(shinyjs)
#library(bslib)

# Graph ploitting
library(ggplot2)
library(plotly)
library(Seurat)
library(viridis)
library(DT)
#library(patchwork)
#library(stringr)

# Used to rename DF columns
library(dplyr)



### PREPARATION ###

# Colour palettes
colour_palette <- c('magma', 'inferno', 'plasma', 'viridis', 'cividis', 'rocket', 'mako', 'turbo')

### Preparation for original figure 2A
original2A_DF <- read.csv("./files/Original2A.csv", header = TRUE, row.names = 1)

### Preparation for original figure 2B (Zones)
original2B_Zones_DF <- read.csv("./files/Original2B_Zones.csv", header = TRUE, row.names = 1)

### Preparation for original figure 2B (Genes)
original2B_DF <- read.csv("./files/Figure2_cells.csv", row.names = 1)

### Preparation for original figure 3D, 3E
seurat_organoid <- readRDS("./files/Camp2015_organoid_final.rds")

### Preparation for original figure 3E
genes_3e <- c("FOXG1", "OTX2", "RSPO2", "DCN",
              "ASPM", "LIN28A", "MYT1L", "NEUROD6")

### Preparation for original figure 3F
seurat_3f <- readRDS("./files/Camp2015_fig3f.rds")

genes_3f <- c("FOXG1", "NEUROD6", "OTX2")
#genes_3f <- genes_3f[genes_3f %in% rownames(seurat_3f)]

Idents(seurat_3f) <- factor(seurat_3f$region_group,
                            levels = c("r1","r2","r3","r4","fetal"))


### Preparation for original figure 4A
original4A_DF <- read.csv("./files/Figure4_monocle_coordinates.csv", header = TRUE, row.names = 1)

### Preparation for original figure 4B
original4B_DF <- read.csv("./files/Figure4_cells.csv", row.names = 1)


### Preparation for alternative figure 2A and 2B (Zones)
alt2a <- read.csv("./files/Alt2A.csv", header = TRUE, row.names = 1)

### Preparation for alternative figure 2A and 2B (Zones)
alt2b_zones <- read.csv("./files/Alt2B_zones.csv", header = TRUE, row.names = 1)

### Preparation for alternative figure 2B (Genes)
alt2B_DF <- read.csv("./files/Figure2_Slingshot_coordinates.csv", row.names = 1)

### Preparation for alternative figure 3D
# Load alternate scanpy pipeline output to plot alt tsne graph
altpipeline_tsne <- read.csv("./files/export_tsne.csv")

### Preparation for alternate figure 3E
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

FOXG1_Alt_3F <- dplyr::rename(FOXG1_Alt_3F, expr = FOXG1, region_group = group)
NEUROD6_Alt_3F <- dplyr::rename(NEUROD6_Alt_3F, expr = NEUROD6, region_group = group)
OTX2_Alt_3F <- dplyr::rename(OTX2_Alt_3F, expr = OTX2, region_group = group)

### Preparation for alternative figure 4A
alt4a <- read.csv("./files/Alt4A.csv", header = TRUE, row.names = 1)

### Preparation for alternative figure 4B 
alt4B_DF <- read.csv("./files/Figure4_Slingshot_coordinates.csv", row.names = 1)

# All figure legends are saved into a txt, this can be loaded into a dataframe with each row representing a different legend
figure_legend <- read.table("./files/figure_legends.txt", sep="\n",header=F)

# Text to be placed at bottom of each figure tab to credit original paper
citation <- read.table("./files/citations.txt", sep="\n", header=F)

# Reading in addition csvs
Original3D_DF <- read.csv("./files/Original3D_DF.csv", header = TRUE)
Original3D_DF <- dplyr::rename(Original3D_DF, ID = X)
Original3E_DF <- read.csv("./files/Original3E_DF.csv", header = TRUE, row.names = 1)
Original3F_DF_FOXG1 <- read.csv("./files/Original3F_DF_FOGG1.csv", header = TRUE, row.names = 1)
Original3F_DF_NEUROD6 <- read.csv("./files/Original3F_DF_NEUROD6.csv", header = TRUE, row.names = 1)
Original3F_DF_OTX2 <- read.csv("./files/Original3F_DF_OTX2.csv", header = TRUE, row.names = 1)
alt3E_DF <- read.csv("./files/alt3E_DF.csv", header = TRUE, row.names = 1)



### UI ###

### ui tagList gives a nice page layout
ui <- tagList(
  useShinyjs(),
  navbarPage(
    
    # Website title
    "Group B",
    
    # Final tab to give credit to group members and original paper authors
    tabPanel("Home",
             fluidPage(
               fluidRow(
               column(6, paste("Welcome to the webpage for 2026s Group-B Group Project.")),
               ),
               fluidRow(
                 column(6, paste("This project aims to replicate the findings of Camp et al. (2015). Camp, J.G. et al. (2015). Human cerebral organoids recapitulate gene expression programs of fetal neocortex development. PNAS, 112(51), pp.15672–15677.")),
               ),
               fluidRow(
               column(6, paste("Tabs are available for either replicated original pipeline figures, or alternate pipeline figures.")),
               ),
               fluidRow(
                 column(6, paste("Figures Replicated: 2A, 2B, 3D, 3E, 3F, 4A, 4B.")),
               ),
               fluidRow(
               column(6, paste("This research used the ALICE High Performance Computing facility at the University of Leicester."))
             ))),
    
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
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                            column(2, paste(figure_legend[2,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                 
                 tabPanel("Figure 2B (Genes)",
                          fluidRow(

                            # Select colour palette
                            column(2, radioButtons(
                              "og_2B_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),

                            # Plot figure
                            column(5, plotlyOutput("originalfig_2B")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[3,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          fluidRow(
                            column(2, selectizeInput( 
                              "original2B_featureSelect", 
                              "Select Gene:", 
                              choices = c("SOX2", "EOMES", "MYT1L"),
                              multiple = FALSE)),
                            
                            column(2, sliderInput(
                              "og_2B_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                            
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_2B_Table"))
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
                              column(2, paste(figure_legend[4,])),
                              
                              column(2, actionButton("Reset_Button", "Reset All Inputs"))
                            ),
                            fluidRow(
                              # Slider for data point opacity
                              column(2, sliderInput(
                                "og_3D_alpha", "Data point opacity",
                                min = 0, max = 1,
                                value = 1
                              )),
                              column(2, selectizeInput(
                                "Original_3D_highlight",
                                "Highlight Cell ID",
                                Original3D_DF$ID,
                                multiple = TRUE
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
                            column(2, paste(figure_legend[5,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                            column(2, paste(figure_legend[6,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                 
                 # Original pipeline figure 4A graph
                 tabPanel("Figure 4A",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "og_4A_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("originalfig_4A")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[7,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "og_4A_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_4A_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,]))
                          ),),
                 
                 
                 # Original figure 4B
                 tabPanel("Figure 4B",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "og_4B_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("originalfig_4B")),
                            
                            # Display figure legend
                            column(2, paste(figure_legend[8,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          fluidRow(
                            column(2, selectizeInput( 
                              "original4B_featureSelect", 
                              "Select Gene:", 
                              choices = c("PAX6", "EOMES", "MYT1L"),
                              multiple = FALSE)),
                            
                            column(2, sliderInput(
                              "og_4B_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                            
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("Og_4B_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[1,]))
                          ),),

               ))
             ),

    # Alternative pipeline Tab
    tabPanel("Alternative Pipeline",
             fluidPage(
               tabsetPanel(
                 
                 # Alternative pipeline figure 2A
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
                            column(2, paste(figure_legend[9,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                            column(12, paste(citation[2,]))
                          ),),
                 
                 # Alternative pipeline figure 2B
                 tabPanel("Figure 2B (Zones)",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_2B_Zone_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("alt_fig_2B_zone")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[10,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_2B_zone_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("alt_2B_zone_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,]))
                          ),),
                 
                 # Alternative pipeline figure 2B
                 tabPanel("Figure 2B (Genes)",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_2B_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),

                            # Plot figure
                            column(5, plotlyOutput("altfig_2B")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[11,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          
                          fluidRow(
                            
                            # Select gene
                            column(2, selectizeInput( 
                              "alternate2B_featureSelect", 
                              "Select Gene:", 
                              choices = c("SOX2", "EOMES", "MYT1L"),
                              multiple = FALSE)),
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_2B_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                            
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("alt_2B_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,]))
                          ),),

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
                            column(2, paste(figure_legend[12,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                            column(2, paste(figure_legend[13,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                            column(2, paste(figure_legend[14,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
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
                 
                 # Alternative pipeline figure 4A
                 tabPanel("Figure 4A",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_4A_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("altfig_4A")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[15,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          
                          fluidRow(
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_4A_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("alt_4A_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,]))
                          ),),
                 
                 # Alternative pipeline figure 4B
                 tabPanel("Figure 4B",
                          fluidRow(
                            
                            # Select colour palette
                            column(2, radioButtons(
                              "alt_4B_palette", "Select Colour Palette",
                              choices = colour_palette, selected = colour_palette[4],
                            )),
                            
                            # Plot figure
                            column(5, plotlyOutput("altfig_4B")),
                            
                            # # Display figure legend
                            column(2, paste(figure_legend[16,])),
                            
                            column(2, actionButton("Reset_Button", "Reset All Inputs"))
                          ),
                          
                          fluidRow(
                            
                            # Select gene
                            column(2, selectizeInput( 
                              "alternate4B_featureSelect", 
                              "Select Gene:", 
                              choices = c("PAX6", "EOMES", "MYT1L"),
                              multiple = FALSE)),
                            
                            # Slider for data point opacity
                            column(2, sliderInput(
                              "alt_4B_alpha", "Data point opacity",
                              min = 0, max = 1,
                              value = 1
                            ))
                            
                          ),
                          
                          fluidRow(
                            column(6, paste("Data Table"))
                          ),
                          
                          # Show data frame used to plot figure
                          fluidRow(
                            column(12, dataTableOutput("alt_4B_Table"))
                          ),
                          
                          # Bottom of page original paper credit
                          fluidRow(
                            column(12, paste(citation[2,]))
                          ),),
               ))),
))








### SERVER ###

server <- function(input, output) {
  
  observeEvent(input$Reset_Button, {
    reset()
  })
  
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
  
  ### Original pipeline 2B (Genes)
  
  og_fig2b_df <- reactive({
    switch(input$original2B_featureSelect,
           "SOX2" = original2B_DF$SOX2,
           "EOMES" = original2B_DF$EOMES,
           "MYT1L" = original2B_DF$MYT1L)
  })
  
  output$Og_2B_Table <- renderDataTable({datatable(original2B_DF)})
  
  output$originalfig_2B <- renderPlotly({
    req(og_fig2b_df())
    plot_ly(original2B_DF, 
            x = original2B_DF$Component_1,
            y = original2B_DF$Component_2,
            # Colour data points based on their zone
            color = og_fig2b_df(),
            
            # Get user input for colour palette 
            colors = viridis_pal(option = input$og_2B_palette)(10),
            type = "scatter",
            mode = "markers",
            alpha = input$og_2B_alpha
    )
    
  })
  
  ### Original pipeline figure 3D graph
  # Show interactive datatable for data plotted in original fig 3D
  output$Og_3D_Table <- renderDataTable({datatable(Original3D_DF)}) 
  
  test <- c("!","£")
  # Plot original 3D graph
  output$originalfig_3D <- renderPlotly({
    fig3d <- DimPlot(seurat_organoid,
                     reduction = "tsne",
                     label = TRUE,
                     label.size = 3,
                     pt.size = 1.5,
                     repel = TRUE,
                     cells.highlight = input$Original_3D_highlight,
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
  output$Og_3E_Table <- renderDataTable({datatable(Original3E_DF)}) 
  
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
  og_3F_options <- reactive({
    switch(input$original3F_featureSelect,
           "FOXG1" = Original3F_DF_FOXG1,
           "NEUROD6" = Original3F_DF_NEUROD6,
           "OTX2" = Original3F_DF_OTX2)
  })
  
  output$Og_3F_Table <- renderDataTable({datatable(og_3F_options())}) 
  
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
  
  ### Original pipeline 4A graph
  # Show interactive datatable for data plotted in original fig 2A
  output$Og_4A_Table <- renderDataTable({datatable(original4A_DF)})
  
  # Plot original 4A graph
  output$originalfig_4A <- renderPlotly({
    original_4A_obj <- plot_ly(original4A_DF, 
                               x = original4A_DF$Component_1,
                               y = original4A_DF$Component_2,
                               
                               # Colour data points based on their label
                               color = original4A_DF$paper_like_label,
                               
                               # User input for colour palette
                               colors = viridis_pal(option = input$og_4A_palette)(7),
                               type = "scatter",
                               mode = "markers",
                               
                               # User input for data point opacity
                               alpha = input$og_4A_alpha
    ) |>
      layout(xaxis = list(title = "Component 1"),
             yaxis = list(title = "Component 2"))
    
  })
  
  ### Original pipeline 4B 
  og_fig4b_df <- reactive({
    switch(input$original4B_featureSelect,
           "PAX6" = original4B_DF$PAX6,
           "EOMES" = original4B_DF$EOMES,
           "MYT1L" = original4B_DF$MYT1L)
  })
  
  output$Og_4B_Table <- renderDataTable({datatable(original4B_DF)})
  
  output$originalfig_4B <- renderPlotly({
    req(og_fig2b_df())
    plot_ly(original4B_DF, 
            x = original4B_DF$Component_1,
            y = original4B_DF$Component_2,
            # Colour data points based on their zone
            color = og_fig4b_df(),
            
            # Get user input for colour palette 
            colors = viridis_pal(option = input$og_4B_palette)(10),
            type = "scatter",
            mode = "markers",
            alpha = input$og_4B_alpha
    )
    
  })
  
  ### ALTERNATE GRAPHS ###
  
  ### Alternate pipeline 2A
  # Display data table for alternate 2A data
  output$alt_2A_Table <- renderDataTable({datatable(alt2a)}) 
  # Plot alternate 2A figure
  output$altfig_2A <- renderPlotly({
    alt_2a_obj <- plot_ly(alt2a, 
                           x = alt2a$X,
                           y = alt2a$Y,
                           
                           # Colour data points based on their label
                           color = alt2a$paper_state,
                           
                           # User input for colour palette
                           colors = viridis_pal(option = input$alt_2A_palette)(7),
                           type = "scatter",
                           mode = "markers",
                           
                           # User input for data point opacity
                           alpha = input$alt_2A_alpha
    ) |>
    # Add x and y axis titles
    layout(xaxis = list(title = "Component 1"),
           yaxis = list(title = "Component 2"))
    
  })
  
  ###Alternate pipeline 2B
  # Display data table for alternate 2A data
  output$alt_2B_zone_Table <- renderDataTable({datatable(alt2b_zones)}) 
  # Plot alternate 2A figure
  output$alt_fig_2B_zone <- renderPlotly({
    alt_2b_zone_obj <- plot_ly(alt2b_zones , 
                          x = alt2b_zones$X,
                          y = alt2b_zones$Y,
                          
                          # Colour data points based on their label
                          color = alt2b_zones$Zone,
                          
                          # User input for colour palette
                          colors = viridis_pal(option = input$alt_2B_Zone_palette)(7),
                          type = "scatter",
                          mode = "markers",
                          
                          # User input for data point opacity
                          alpha = input$alt_2B_zone_alpha
    ) |>
      # Add x and y axis titles
      layout(xaxis = list(title = "Component 1"),
             yaxis = list(title = "Component 2"))
    
  })
  
  ### Alternate pipeline 2B (Genes)
  
  alt_fig2b_df <- reactive({
    switch(input$alternate2B_featureSelect,
           "SOX2" = alt2B_DF$SOX2,
           "EOMES" = alt2B_DF$EOMES,
           "MYT1L" = alt2B_DF$MYT1L)
  })
  
  output$alt_2B_Table <- renderDataTable({datatable(alt2B_DF)})
  
  output$altfig_2B <- renderPlotly({
    req(alt_fig2b_df())
    plot_ly(alt2B_DF, 
            x = alt2B_DF$X,
            y = alt2B_DF$Y,
            # Colour data points based on their zone
            color = alt_fig2b_df(),
            
            # Get user input for colour palette 
            colors = viridis_pal(option = input$alt_2B_palette)(10),
            type = "scatter",
            mode = "markers",
            alpha = input$alt_2B_alpha
    )
    
  })
  
  ### Alternate pipeline 3D
  
  # Display datatable for alternate 3D data
  output$Alt_3D_Table <- renderDataTable({datatable(altpipeline_tsne)}) 
  # Plot alternative pipline figure 3D
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
  output$Alt_3E_Table <- renderDataTable({datatable(alt3E_DF)}) 
  
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
      
    ) + RotatedAxis() +
      scale_color_viridis_c(option = input$alt_3E_palette)

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
  
  ### Alternate pipeline 4A
  # Display data table for alternate 4A data
  output$alt_4A_Table <- renderDataTable({datatable(alt4a)}) 
  # Plot alternate 4A figure
  output$altfig_4A <- renderPlotly({
    alt_4a_obj <- plot_ly(alt4a, 
                          x = alt4a$X,
                          y = alt4a$Y,
                          
                          # Colour data points based on their label
                          color = alt4a$paper_state,
                          
                          # User input for colour palette
                          colors = viridis_pal(option = input$alt_4A_palette)(7),
                          type = "scatter",
                          mode = "markers",
                          
                          # User input for data point opacity
                          alpha = input$alt_4A_alpha
    ) |>
      # Add x and y axis titles
      layout(xaxis = list(title = "Component 1"),
             yaxis = list(title = "Component 2"))
    
  })
  
  ### Alternate pipeline 4B
  
  alt_fig4b_df <- reactive({
    switch(input$alternate4B_featureSelect,
           "PAX6" = alt4B_DF$PAX6,
           "EOMES" = alt4B_DF$EOMES,
           "MYT1L" = alt4B_DF$MYT1L)
  })
  
  output$alt_4B_Table <- renderDataTable({datatable(alt4B_DF)})
  
  output$altfig_4B <- renderPlotly({
    req(alt_fig4b_df())
    plot_ly(alt4B_DF, 
            x = alt4B_DF$X,
            y = alt4B_DF$Y,
            # Colour data points based on their zone
            color = alt_fig4b_df(),
            
            # Get user input for colour palette 
            colors = viridis_pal(option = input$alt_4B_palette)(10),
            type = "scatter",
            mode = "markers",
            alpha = input$alt_4B_alpha
    )
    
  })
}

shinyApp(ui = ui, server = server)
