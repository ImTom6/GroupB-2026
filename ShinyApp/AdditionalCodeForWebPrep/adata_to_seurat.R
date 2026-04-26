# Code to take input anndata file (.h5ad) and output Seurat Object (.rds). 
# Last Edited: Tom 25/04/26

#library(anndataR)
library(Seurat)
library(anndata)
#library(reticulate)
#library(Matrix)

data <- read_h5ad("./ScanPy/output.h5ad")


# Converting object to seurat object, using same filtering as used in scanpy pipeline
data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 200, min.cells = 3)

# Ensuring data is normalised to avoid errors when plotting dotplot
data <- NormalizeData(data)

saveRDS(data,"./files/scpy_data.rds")

#DotPlot(
#  data,
#  features = markers_present,
#  group.by = "cluster"
#)
