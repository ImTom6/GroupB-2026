# Converting dataframes to be shown in data tables to user
# Last edited 26/04/26: TOM

library(Seurat)

seurat_organoid <- readRDS("./files/Camp2015_organoid_final.rds")
seurat_3f <- readRDS("./files/Camp2015_fig3f.rds")

Idents(seurat_3f) <- factor(seurat_3f$region_group,
                            levels = c("r1","r2","r3","r4","fetal"))

# Extracting key data
coords <- Embeddings(seurat_organoid, "tsne")
meta <- seurat_organoid@meta.data

Original3F_DF <- FetchData(seurat_3f, vars = c("FOXG1", "NEUROD6", "OTX2", "region_group"))
                                            
# Needed for 3E
exp_data <- FetchData(seurat_organoid, vars = c("FOXG1", "OTX2", "RSPO2", "DCN", "ASPM", "LIN28A",
                                                "MYT1L", "NEUROD6", "seurat_clusters"))

# Building DF needed to be shown for original 3D figure
Original3D_DF <- data.frame(
  cluster =  meta$seurat_clusters,
  x = coords[,1],
  y = coords[,2]
)
write.csv(Original3D_DF, file = "./files/Original3D_DF.csv", row.names = TRUE)

# Building DF for original 3E
Original3E_DF <- data.frame(
  cluster = exp_data$seurat_clusters,
  x = coords[,1],
  y = coords[,2],
  exp_data
)
write.csv(Original3E_DF, file = "./files/Original3E_DF.csv", row.names = TRUE)

# Building Original 3F Dataframes 
Original3F_DF_FOGG1 <- data.frame(
  id = rownames(Original3F_DF),
  cluster = Original3F_DF$region_group,
  gene_expression = Original3F_DF$FOXG1
)
write.csv(Original3F_DF_FOGG1, file = "./files/Original3F_DF_FOGG1.csv", row.names = FALSE)

Original3F_DF_NEUROD6 <- data.frame(
  id = rownames(Original3F_DF),
  cluster = Original3F_DF$region_group,
  gene_expression = Original3F_DF$NEUROD6
)
write.csv(Original3F_DF_NEUROD6, file = "./files/Original3F_DF_NEUROD6.csv", row.names = FALSE)

Original3F_DF_OTX2 <- data.frame(
  id = rownames(Original3F_DF),
  cluster = Original3F_DF$region_group,
  gene_expression = Original3F_DF$OTX2
)
write.csv(Original3F_DF_OTX2, file = "./files/Original3F_DF_OTX2.csv", row.names = FALSE)



# For alternate plot 3E
scpy_data <- readRDS("./files/scpy_data.rds")

markers_present <- c("FOXG1", "NFIA", "NFIB", "NEUROD6",
                     "BCL11A", "DCX", "OTX2", "FABP7",
                     "BCAT1", "GAD1", "DLX6", "WNT2B",
                     "RSPO2", "RSPO3", "WLS", "COL3A1",
                     "LUM", "DCN", "SPARC")

alt3E <- DotPlot(
  scpy_data,
  features = markers_present,
  group.by = "cluster",
  scale.by = "radius",
  scale = TRUE,
  #dot.min = input$alt_3E_min,
  #dot.scale = input$alt_3E_max
  #cols = viridis_pal(option = input$alt_3E_palette)(2)
  
) + RotatedAxis()
alt3E_DF <-alt3E$data
write.csv(alt3E_DF, file = "./files/alt3E_DF.csv", row.names = TRUE)

# For monocle
library(SingleCellExperiment)
load("./Hari/monocle_figures.RData")

# Prep for original figures
write.csv(plot_df, file = "./files/Original2A.csv", row.names = TRUE)
write.csv(zone_df, file = "./files/Original2B_Zones.csv", row.names = TRUE)


load("./slingshot_figures.RData")

`# Prep for alt figures
write.csv(plot_df_fetal, file = "./files/Alt2A.csv", row.names = TRUE)
write.csv(plot_df_fetal, file = "./files/Alt2B_zones.csv", row.names = TRUE)
write.csv(plot_df_org, file = "./files/Alt4A.csv", row.names = TRUE)

