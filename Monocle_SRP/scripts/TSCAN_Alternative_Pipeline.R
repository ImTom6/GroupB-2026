# Alternative trajectory analysis using TSCAN
# Generates Figure 2A/2B and Figure 4A/4B
# Saves .plot.rds files instead of saving the full sce object

# loading required libraries
library(SingleCellExperiment)
library(TSCAN)
library(scater)
library(S4Vectors)
library(SummarizedExperiment)
library(ggplot2)
library(gridExtra)
library(grid)
library(igraph)

# file path
file_path <- "/home/haripriya/Downloads/expression_matrix_full.csv"

# reading the file
expr_raw <- read.csv(
  file_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("Loaded matrix dimensions:", dim(expr_raw), "\n")

# checking required columns
if (!all(c("gene_id", "gene_name") %in% colnames(expr_raw))) {
  stop("Input must contain gene_id and gene_name columns.")
}

gene_id <- expr_raw$gene_id
gene_name <- expr_raw$gene_name

# keeping only required expression values
expr_df <- expr_raw[, !(colnames(expr_raw) %in% c("gene_id", "gene_name")), drop = FALSE]
expr_df[] <- lapply(expr_df, function(x) as.numeric(as.character(x)))
expr_df[is.na(expr_df)] <- 0

expr_mat <- as.matrix(expr_df)

# fixing gene names
clean_gene_name <- gene_name
missing_name_idx <- is.na(clean_gene_name) | clean_gene_name == ""
clean_gene_name[missing_name_idx] <- gene_id[missing_name_idx]

dup_idx <- duplicated(clean_gene_name)
clean_gene_name[dup_idx] <- paste0(clean_gene_name[dup_idx], "_", gene_id[dup_idx])

rownames(expr_mat) <- clean_gene_name

# filtering low-expression genes
expr_mat <- expr_mat[rowSums(expr_mat) > 0, , drop = FALSE]
expr_mat <- expr_mat[rowSums(expr_mat > 0) >= 3, , drop = FALSE]

cat("After filtering:", dim(expr_mat), "\n")

# log transform
expr_mat_log <- log2(expr_mat + 1)

cat("Total genes after filtering:", nrow(expr_mat_log), "\n")
cat("Total cells after filtering:", ncol(expr_mat_log), "\n")

# helper function to extract gene expression
get_gene_safe <- function(mat, g) {
  if (g %in% rownames(mat)) {
    return(as.numeric(mat[g, ]))
  } else {
    return(rep(0, ncol(mat)))
  }
}

# setting colours for cell states and cortical zones
paper_colors <- c(
  "AP1" = "#4C7F35",
  "AP2" = "#9FD37B",
  "BP1" = "#58C7C4",
  "BP2" = "#BFE7E3",
  "N1"  = "#8AB7E8",
  "N2"  = "#5C8AD8",
  "N3"  = "#2949A8"
)

zone_colors <- c(
  "VZ"   = "#F3C843",
  "iSVZ" = "#E7952F",
  "oSVZ" = "#C9473A",
  "CP"   = "#8E44AD"
)

# selecting fetal-like cells
pecam1 <- get_gene_safe(expr_mat_log, "PECAM1")
gad1   <- get_gene_safe(expr_mat_log, "GAD1")
dlx1   <- get_gene_safe(expr_mat_log, "DLX1")
dlx2   <- get_gene_safe(expr_mat_log, "DLX2")
dlx5   <- get_gene_safe(expr_mat_log, "DLX5")
dlx6   <- get_gene_safe(expr_mat_log, "DLX6")
erbb4  <- get_gene_safe(expr_mat_log, "ERBB4")

interneuron_score <- gad1 + dlx1 + dlx2 + dlx5 + dlx6 + erbb4
endothelial_score <- pecam1

keep_cells_fetal <- (interneuron_score < 1.5) & (endothelial_score < 1.5)
expr_mat_log_fetal <- expr_mat_log[, keep_cells_fetal, drop = FALSE]

cat("Fetal-like cells retained:", ncol(expr_mat_log_fetal), "\n")

# selecting organoid-like cells
foxg1   <- get_gene_safe(expr_mat_log, "FOXG1")
otx2    <- get_gene_safe(expr_mat_log, "OTX2")
nfia    <- get_gene_safe(expr_mat_log, "NFIA")
nfib    <- get_gene_safe(expr_mat_log, "NFIB")
neurod6 <- get_gene_safe(expr_mat_log, "NEUROD6")

dorsal_score <- foxg1 + nfia + nfib + neurod6 - otx2
threshold <- quantile(dorsal_score, 0.60, na.rm = TRUE)

keep_cells_org <- dorsal_score >= threshold
expr_mat_log_org <- expr_mat_log[, keep_cells_org, drop = FALSE]

cat("Organoid-like cells retained:", ncol(expr_mat_log_org), "\n")

# running TSCAN trajectory analysis
build_tscan <- function(mat_log, n_clusters = 6, marker_genes = c()) {
  
  # removing genes with zero variance
  gene_var <- apply(mat_log, 1, var, na.rm = TRUE)
  mat_log <- mat_log[gene_var > 0, , drop = FALSE]
  
  # removing empty cells
  cell_sum <- colSums(mat_log, na.rm = TRUE)
  mat_log <- mat_log[, cell_sum > 0, drop = FALSE]
  
  if (nrow(mat_log) < 10) stop("Too few variable genes remain.")
  if (ncol(mat_log) < 10) stop("Too few cells remain.")
  
  # creating SingleCellExperiment object temporarily
  sce <- SingleCellExperiment(
    assays = list(logcounts = mat_log)
  )
  
  # running PCA
  pca <- prcomp(
    t(mat_log),
    center = TRUE,
    scale. = FALSE
  )
  
  if (ncol(pca$x) < 2) {
    stop("PCA did not return enough components.")
  }
  
  reducedDims(sce) <- SimpleList(
    PCA = pca$x[, 1:2, drop = FALSE]
  )
  
  # clustering cells
  n_pc <- min(10, ncol(pca$x))
  pc_use <- pca$x[, 1:n_pc, drop = FALSE]
  
  set.seed(42)
  km <- kmeans(pc_use, centers = n_clusters)
  colData(sce)$cluster <- factor(km$cluster)
  
  # calculating cluster centres
  coords <- reducedDim(sce, "PCA")
  cluster_levels <- levels(colData(sce)$cluster)
  
  centroids <- do.call(rbind, lapply(cluster_levels, function(cl) {
    colMeans(coords[colData(sce)$cluster == cl, , drop = FALSE])
  }))
  
  rownames(centroids) <- cluster_levels
  
  # building minimum spanning tree between clusters
  mst <- TSCAN::createClusterMST(centroids, clusters = NULL)
  
  # choosing start cluster using progenitor markers
  ap_markers <- intersect(c("SOX2", "PAX6", "HES1", "VIM", "PROM1"), rownames(mat_log))
  ap_score <- rep(0, ncol(mat_log))
  
  if (length(ap_markers) > 0) {
    ap_score <- colMeans(mat_log[ap_markers, , drop = FALSE])
  }
  
  cluster_means <- tapply(ap_score, colData(sce)$cluster, mean)
  start_cluster <- names(which.max(cluster_means))
  
  # assigning pseudotime using distance from start cluster
  sp <- igraph::distances(mst, v = start_cluster)
  cluster_pt <- as.numeric(sp[1, as.character(colData(sce)$cluster)])
  
  local_pt <- sapply(seq_len(nrow(coords)), function(i) {
    cl <- as.character(colData(sce)$cluster[i])
    sqrt(sum((coords[i, ] - centroids[cl, ])^2))
  })
  
  pseudotime <- cluster_pt + local_pt
  
  # making plotting dataframe
  plot_df <- data.frame(
    cell_id = colnames(sce),
    X = coords[, 1],
    Y = coords[, 2],
    cluster = as.factor(colData(sce)$cluster),
    Pseudotime = pseudotime,
    stringsAsFactors = FALSE
  )
  
  # rotating coordinates for cleaner plotting
  rot <- prcomp(plot_df[, c("X", "Y")], center = TRUE, scale. = FALSE)
  rot_coords <- rot$x[, 1:2]
  
  plot_df$X <- rot_coords[, 1]
  plot_df$Y <- rot_coords[, 2]
  
  # applying same rotation to cluster centres
  cent_scaled <- scale(
    centroids,
    center = rot$center,
    scale = FALSE
  )
  
  cent_rot <- cent_scaled %*% rot$rotation[, 1:2, drop = FALSE]
  
  cent_plot <- data.frame(
    cluster = rownames(centroids),
    x = cent_rot[, 1],
    y = cent_rot[, 2],
    stringsAsFactors = FALSE
  )
  
  # making early cells appear on the left
  cor_test <- cor(plot_df$X, plot_df$Pseudotime, use = "complete.obs")
  
  if (!is.na(cor_test) && cor_test < 0) {
    plot_df$X <- -plot_df$X
    cent_plot$x <- -cent_plot$x
  }
  
  # adjusting vertical direction
  q_low <- quantile(plot_df$Pseudotime, 0.1, na.rm = TRUE)
  q_high <- quantile(plot_df$Pseudotime, 0.9, na.rm = TRUE)
  
  early_mean <- mean(plot_df$Y[plot_df$Pseudotime < q_low], na.rm = TRUE)
  late_mean <- mean(plot_df$Y[plot_df$Pseudotime > q_high], na.rm = TRUE)
  
  if (!is.na(early_mean) && !is.na(late_mean) && late_mean < early_mean) {
    plot_df$Y <- -plot_df$Y
    cent_plot$y <- -cent_plot$y
  }
  
  # making paper-style AP / BP / N groups
  q <- quantile(
    plot_df$Pseudotime,
    probs = c(0, 0.14, 0.28, 0.45, 0.62, 0.78, 0.90, 1),
    na.rm = TRUE
  )
  
  if (length(unique(q)) < 8) {
    ranks <- rank(plot_df$Pseudotime, ties.method = "first") / length(plot_df$Pseudotime)
    plot_df$paper_state <- cut(
      ranks,
      breaks = c(0, 0.14, 0.28, 0.45, 0.62, 0.78, 0.90, 1),
      include.lowest = TRUE,
      labels = c("AP1", "AP2", "BP1", "BP2", "N1", "N2", "N3")
    )
  } else {
    plot_df$paper_state <- cut(
      plot_df$Pseudotime,
      breaks = q,
      include.lowest = TRUE,
      labels = c("AP1", "AP2", "BP1", "BP2", "N1", "N2", "N3")
    )
  }
  
  plot_df$paper_state <- factor(
    plot_df$paper_state,
    levels = c("AP1", "AP2", "BP1", "BP2", "N1", "N2", "N3")
  )
  
  # making cortical zone groups
  qz <- quantile(
    plot_df$Pseudotime,
    probs = c(0, 0.25, 0.50, 0.75, 1),
    na.rm = TRUE
  )
  
  if (length(unique(qz)) < 5) {
    ranks_z <- rank(plot_df$Pseudotime, ties.method = "first") / length(plot_df$Pseudotime)
    plot_df$Zone <- cut(
      ranks_z,
      breaks = c(0, 0.25, 0.50, 0.75, 1),
      include.lowest = TRUE,
      labels = c("VZ", "iSVZ", "oSVZ", "CP")
    )
  } else {
    plot_df$Zone <- cut(
      plot_df$Pseudotime,
      breaks = qz,
      include.lowest = TRUE,
      labels = c("VZ", "iSVZ", "oSVZ", "CP")
    )
  }
  
  # adding marker expression directly into the plotting dataframe
  for (g in marker_genes) {
    if (g %in% rownames(mat_log)) {
      plot_df[[g]] <- as.numeric(mat_log[g, plot_df$cell_id])
    } else {
      plot_df[[g]] <- 0
    }
  }
  
  # making MST line segments from rotated cluster centres
  edges <- as.data.frame(igraph::as_edgelist(mst))
  colnames(edges) <- c("from", "to")
  
  segs <- merge(edges, cent_plot, by.x = "from", by.y = "cluster", all.x = TRUE)
  segs <- merge(
    segs,
    cent_plot,
    by.x = "to",
    by.y = "cluster",
    all.x = TRUE,
    suffixes = c("_from", "_to")
  )
  
  segs <- segs[complete.cases(segs[, c("x_from", "y_from", "x_to", "y_to")]), ]
  
  # returning only plot-ready data
  list(
    plot_df = plot_df,
    segs = segs
  )
}

# running TSCAN for Figure 2 and Figure 4
res_fetal <- build_tscan(
  expr_mat_log_fetal,
  n_clusters = 6,
  marker_genes = c("SOX2", "EOMES", "MYT1L")
)

plot_df_fetal <- res_fetal$plot_df
segs_fetal <- res_fetal$segs

res_org <- build_tscan(
  expr_mat_log_org,
  n_clusters = 6,
  marker_genes = c("PAX6", "EOMES", "MYT1L")
)

plot_df_org <- res_org$plot_df
segs_org <- res_org$segs

# making trajectory plots from dataframe
make_traj_plot <- function(plot_df, segs, title_text) {
  
  start_point <- plot_df[which.min(plot_df$Pseudotime), ]
  end_point <- plot_df[which.max(plot_df$Pseudotime), ]
  
  ggplot() +
    geom_segment(
      data = segs,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      color = "grey70",
      linewidth = 0.7,
      alpha = 0.9
    ) +
    geom_point(
      data = plot_df,
      aes(X, Y, color = paper_state),
      size = 2.3,
      alpha = 0.95
    ) +
    scale_color_manual(values = paper_colors, drop = FALSE) +
    theme_classic(base_size = 12) +
    labs(
      x = "Component 1",
      y = "Component 2",
      color = NULL,
      title = title_text
    ) +
    annotate(
      "text",
      x = start_point$X,
      y = start_point$Y,
      label = "start",
      hjust = -0.1,
      vjust = -0.6,
      size = 4
    ) +
    annotate(
      "text",
      x = end_point$X,
      y = end_point$Y,
      label = "end",
      hjust = -0.1,
      vjust = -0.6,
      size = 4
    ) +
    theme(
      legend.position = c(0.18, 0.95),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

# making marker expression plots from dataframe
make_marker_panel_from_df <- function(plot_df, gene_symbol) {
  
  if (!(gene_symbol %in% colnames(plot_df))) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste(gene_symbol, "not found"), size = 4) +
        theme_void() +
        ggtitle(gene_symbol)
    )
  }
  
  ggplot(plot_df, aes(X, Y, color = .data[[gene_symbol]])) +
    geom_point(size = 1.6, alpha = 0.95) +
    scale_color_gradient(low = "#1E2B7A", high = "#E29A47") +
    theme_void(base_size = 11) +
    ggtitle(gene_symbol) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "bottom"
    )
}

# making cortical zone plot
make_zone_panel <- function(plot_df, title_text = "Zone") {
  ggplot(plot_df, aes(X, Y, color = Zone)) +
    geom_point(size = 1.6, alpha = 0.95) +
    scale_color_manual(values = zone_colors, drop = FALSE) +
    theme_void(base_size = 11) +
    ggtitle(title_text) +
    theme(legend.position = "bottom")
}

# generating Figure 2A and Figure 2B
fig2A <- make_traj_plot(
  plot_df_fetal,
  segs_fetal,
  "Figure 2A-style (TSCAN)"
)

fig2B <- arrangeGrob(
  make_zone_panel(plot_df_fetal, "Zone"),
  make_marker_panel_from_df(plot_df_fetal, "SOX2"),
  make_marker_panel_from_df(plot_df_fetal, "EOMES"),
  make_marker_panel_from_df(plot_df_fetal, "MYT1L"),
  ncol = 4
)

fig2_combined <- arrangeGrob(
  fig2A,
  fig2B,
  nrow = 2,
  heights = c(2.2, 1)
)

grid.newpage()
grid.draw(fig2_combined)

ggsave(
  "Figure2A_TSCAN.png",
  plot = fig2A,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

png("Figure2B_TSCAN.png", width = 2600, height = 700, res = 250)
grid.draw(fig2B)
dev.off()

png("Figure2AB_TSCAN_combined.png", width = 2600, height = 1800, res = 250)
grid.draw(fig2_combined)
dev.off()

# generating Figure 4A and Figure 4B
fig4A <- make_traj_plot(
  plot_df_org,
  segs_org,
  "Figure 4A-style (TSCAN)"
)

fig4B <- arrangeGrob(
  make_zone_panel(plot_df_org, "Zone"),
  make_marker_panel_from_df(plot_df_org, "PAX6"),
  make_marker_panel_from_df(plot_df_org, "EOMES"),
  make_marker_panel_from_df(plot_df_org, "MYT1L"),
  ncol = 4
)

fig4_combined <- arrangeGrob(
  fig4A,
  fig4B,
  nrow = 2,
  heights = c(2.2, 1)
)

grid.newpage()
grid.draw(fig4_combined)

ggsave(
  "Figure4A_TSCAN.png",
  plot = fig4A,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

png("Figure4B_TSCAN.png", width = 2600, height = 700, res = 250)
grid.draw(fig4B)
dev.off()

png("Figure4AB_TSCAN_combined.png", width = 2600, height = 1800, res = 250)
grid.draw(fig4_combined)
dev.off()

# saving coordinate tables as CSV files
write.csv(
  plot_df_fetal,
  "Figure2_TSCAN_coordinates.csv",
  row.names = FALSE
)

write.csv(
  plot_df_org,
  "Figure4_TSCAN_coordinates.csv",
  row.names = FALSE
)

write.csv(
  segs_fetal,
  "Figure2_TSCAN_segments.csv",
  row.names = FALSE
)

write.csv(
  segs_org,
  "Figure4_TSCAN_segments.csv",
  row.names = FALSE
)

# saving .plot.rds plot objects only
saveRDS(fig2A, "Figure2A_TSCAN.plot.rds")
saveRDS(fig2B, "Figure2B_TSCAN.plot.rds")
saveRDS(fig2_combined, "Figure2AB_TSCAN.plot.rds")

saveRDS(fig4A, "Figure4A_TSCAN.plot.rds")
saveRDS(fig4B, "Figure4B_TSCAN.plot.rds")
saveRDS(fig4_combined, "Figure4AB_TSCAN.plot.rds")

tscan_plot_objects <- list(
  Figure2A_plot = fig2A,
  Figure2B_plot = fig2B,
  Figure2AB_plot = fig2_combined,
  Figure4A_plot = fig4A,
  Figure4B_plot = fig4B,
  Figure4AB_plot = fig4_combined,
  plot_df_fetal = plot_df_fetal,
  plot_df_org = plot_df_org,
  segs_fetal = segs_fetal,
  segs_org = segs_org,
  paper_colors = paper_colors,
  zone_colors = zone_colors
)

saveRDS(
  tscan_plot_objects,
  "TSCAN_all_plot_objects.plot.rds"
)

# checking saved output files
cat("\nDONE\n")

cat("Saved PNG files:\n")
print(list.files(pattern = "Figure[24].*TSCAN.*png"))

cat("\nSaved CSV files:\n")
print(list.files(pattern = "Figure[24]_TSCAN.*csv"))

cat("\nSaved .plot.rds files:\n")
print(list.files(pattern = "TSCAN.*\\.plot\\.rds|Figure[24].*TSCAN.*\\.plot\\.rds"))

cat("\nThe TSCAN plots were saved as PNG and .plot.rds files.\n")
cat("The marker panels are dataframe-based, so no sce object is needed in the website.\n")

cat("\nWorking directory:\n")
print(getwd())