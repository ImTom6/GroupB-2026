library(monocle)
library(ggplot2)
library(gridExtra)
library(igraph)
library(grid)


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
  stop("The input file must contain 'gene_id' and 'gene_name' columns.")
}

gene_id <- expr_raw$gene_id
gene_name <- expr_raw$gene_name


# keeping only expression values
expr_df <- expr_raw[, !(colnames(expr_raw) %in% c("gene_id", "gene_name")), drop = FALSE]
expr_df[] <- lapply(expr_df, function(x) as.numeric(as.character(x)))
expr_df[is.na(expr_df)] <- 0

expr_mat <- as.matrix(expr_df)


# fixing gene names
clean_gene_name <- gene_name
clean_gene_name[is.na(clean_gene_name) | clean_gene_name == ""] <-
  gene_id[is.na(clean_gene_name) | clean_gene_name == ""]

dup_idx <- duplicated(clean_gene_name)
clean_gene_name[dup_idx] <- paste0(clean_gene_name[dup_idx], "_", gene_id[dup_idx])

rownames(expr_mat) <- clean_gene_name


# removing genes with no expression
expr_mat <- expr_mat[rowSums(expr_mat) > 0, , drop = FALSE]


# keeping genes expressed in enough cells
min_cells <- max(10, round(0.05 * ncol(expr_mat)))
genes_keep <- rowSums(expr_mat > 0) >= min_cells
expr_mat <- expr_mat[genes_keep, , drop = FALSE]

cat("After filtering:", dim(expr_mat), "\n")


# log transform
expr_mat_log <- log2(expr_mat + 1)


# making metadata for cells
cell_metadata <- data.frame(
  cell_id = colnames(expr_mat_log),
  row.names = colnames(expr_mat_log),
  stringsAsFactors = FALSE
)


# making metadata for genes
gene_metadata <- data.frame(
  gene_short_name = rownames(expr_mat_log),
  row.names = rownames(expr_mat_log),
  stringsAsFactors = FALSE
)

pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_metadata)


# building the monocle object
cds <- newCellDataSet(
  expr_mat_log,
  phenoData = pd,
  featureData = fd,
  expressionFamily = uninormal()
)

cds <- estimateSizeFactors(cds)
cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= min_cells))
cat("Expressed genes retained:", length(expressed_genes), "\n")

if (length(expressed_genes) < 10) {
  stop("Too few expressed genes retained. Please check the input matrix.")
}


# choosing ordering genes using variance
expr_subset <- exprs(cds)[expressed_genes, , drop = FALSE]
gene_variance <- apply(expr_subset, 1, var)
gene_variance <- sort(gene_variance, decreasing = TRUE)

ordering_genes <- names(gene_variance)[1:min(3000, length(gene_variance))]
cds <- setOrderingFilter(cds, ordering_genes)

cat("Ordering genes selected:", length(ordering_genes), "\n")


# running trajectory
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree",
  norm_method = "none"
)


# force Monocle's orderCells() to avoid masking from other packages
cds <- monocle::orderCells(cds)

cat("Trajectory inference complete.\n")


# reduced dimension coordinates
rd <- reducedDimS(cds)

plot_df <- data.frame(
  cell_id = rownames(pData(cds)),
  Component_1 = rd[1, ],
  Component_2 = rd[2, ],
  State = as.factor(pData(cds)$State),
  Pseudotime = pData(cds)$Pseudotime,
  stringsAsFactors = FALSE
)


# making rough AP / BP / N groups from pseudotime
pt_breaks <- quantile(
  plot_df$Pseudotime,
  probs = c(0, 0.33, 0.66, 1),
  na.rm = TRUE
)

# avoid duplicate breaks issue
pt_breaks <- unique(pt_breaks)

if (length(pt_breaks) < 4) {
  plot_df$stage_group <- "BP_like"
} else {
  plot_df$stage_group <- cut(
    plot_df$Pseudotime,
    breaks = pt_breaks,
    include.lowest = TRUE,
    labels = c("AP_like", "BP_like", "N_like")
  )
}

plot_df$paper_like_label <- NA_character_

for (s in levels(plot_df$State)) {
  idx <- which(plot_df$State == s)
  dominant_group <- names(sort(table(plot_df$stage_group[idx]), decreasing = TRUE))[1]
  
  if (dominant_group == "AP_like") {
    plot_df$paper_like_label[idx] <- "AP"
  } else if (dominant_group == "BP_like") {
    plot_df$paper_like_label[idx] <- "BP"
  } else {
    plot_df$paper_like_label[idx] <- "N"
  }
}

ap_states <- unique(plot_df$State[plot_df$paper_like_label == "AP"])
bp_states <- unique(plot_df$State[plot_df$paper_like_label == "BP"])
n_states  <- unique(plot_df$State[plot_df$paper_like_label == "N"])

ap_map <- setNames(paste0("AP", seq_along(ap_states)), ap_states)
bp_map <- setNames(paste0("BP", seq_along(bp_states)), bp_states)
n_map  <- setNames(paste0("N", seq_along(n_states)), n_states)

for (s in names(ap_map)) plot_df$paper_like_label[plot_df$State == s] <- ap_map[[s]]
for (s in names(bp_map)) plot_df$paper_like_label[plot_df$State == s] <- bp_map[[s]]
for (s in names(n_map))  plot_df$paper_like_label[plot_df$State == s] <- n_map[[s]]

plot_df$paper_like_label <- factor(plot_df$paper_like_label)


# getting start and end points
start_idx <- which.min(plot_df$Pseudotime)
end_idx   <- which.max(plot_df$Pseudotime)

start_point <- plot_df[start_idx, ]
end_point   <- plot_df[end_idx, ]


# function to get trajectory line
get_mst_segments <- function(cds_obj) {
  mst_obj <- minSpanningTree(cds_obj)
  mst_edges <- as.data.frame(igraph::get.edgelist(mst_obj))
  colnames(mst_edges) <- c("from", "to")
  
  cell_coords <- reducedDimS(cds_obj)
  coord_df <- data.frame(
    cell_id = colnames(cds_obj),
    x = cell_coords[1, ],
    y = cell_coords[2, ],
    stringsAsFactors = FALSE
  )
  
  segs <- merge(mst_edges, coord_df, by.x = "from", by.y = "cell_id", all.x = TRUE)
  segs <- merge(
    segs, coord_df,
    by.x = "to", by.y = "cell_id",
    all.x = TRUE,
    suffixes = c("_from", "_to")
  )
  
  segs <- segs[complete.cases(segs[, c("x_from", "y_from", "x_to", "y_to")]), ]
  segs
}

mst_segments <- get_mst_segments(cds)


# colors for panel A
paper_colors <- c
  "AP1" = "#4A7C32",
  "AP2" = "#8CC76A",
  "AP3" = "#A9D98D",
  "BP1" = "#58C7C4",
  "BP2" = "#B7E3E0",
  "BP3" = "#7ED9D3",
  "N1"  = "#7FB4E8",
  "N2"  = "#5F8DD6",
  "N3"  = "#274DAD",
  "N4"  = "#9CC6F0"
)

present_labels <- levels(plot_df$paper_like_label)
present_colors <- paper_colors[present_labels]

# fallback if any labels missing from predefined color vector
missing_cols <- setdiff(present_labels, names(present_colors))
if (length(missing_cols) > 0) {
  extra_cols <- grDevices::rainbow(length(missing_cols))
  names(extra_cols) <- missing_cols
  present_colors <- c(present_colors, extra_cols)
}


# panel A
panel_A <- ggplot() +
  geom_segment(
    data = mst_segments,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    color = "grey70",
    linewidth = 0.4,
    alpha = 0.9
  ) +
  geom_point(
    data = plot_df,
    aes(Component_1, Component_2, color = paper_like_label),
    size = 2.4,
    alpha = 0.95
  ) +
  scale_color_manual(values = present_colors, drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    x = "Component 1",
    y = "Component 2",
    color = NULL
  ) +
  annotate(
    "text",
    x = start_point$Component_1,
    y = start_point$Component_2,
    label = "start",
    hjust = -0.1,
    vjust = -0.7,
    size = 4
  ) +
  annotate(
    "text",
    x = end_point$Component_1,
    y = end_point$Component_2,
    label = "end",
    hjust = -0.1,
    vjust = -0.7,
    size = 4
  ) +
  ggtitle("A") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    legend.position = c(0.18, 0.95),
    legend.background = element_blank(),
    legend.key = element_blank()
  )


# helper for marker plots
make_marker_panel <- function(cds_obj, gene_symbol, title_text, plot_data) {
  if (!(gene_symbol %in% rownames(exprs(cds_obj)))) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste(gene_symbol, "not found"), size = 4) +
        theme_void() +
        ggtitle(title_text)
    )
  }
  
  tmp <- plot_data
  tmp$marker <- as.numeric(exprs(cds_obj)[gene_symbol, rownames(pData(cds_obj))])
  
  ggplot(tmp, aes(Component_1, Component_2, color = marker)) +
    geom_point(size = 1.7, alpha = 0.95) +
    scale_color_gradient(low = "#1E2B7A", high = "#E29A47") +
    theme_void(base_size = 11) +
    ggtitle(title_text) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      legend.position = "bottom"
    )
}


# making rough zone groups from pseudotime
zone_df <- plot_df

zone_breaks <- quantile(
  zone_df$Pseudotime,
  probs = c(0, 0.25, 0.50, 0.75, 1),
  na.rm = TRUE
)

zone_breaks <- unique(zone_breaks)

if (length(zone_breaks) < 5) {
  zone_df$Zone <- "VZ"
} else {
  zone_df$Zone <- cut(
    zone_df$Pseudotime,
    breaks = zone_breaks,
    include.lowest = TRUE,
    labels = c("VZ", "iSVZ", "oSVZ", "CP")
  )
}

zone_colors <- c(
  "VZ" = "#F3C843",
  "iSVZ" = "#E7952F",
  "oSVZ" = "#C9473A",
  "CP" = "#8E44AD"
)


# panel B first plot
panel_B_zone <- ggplot(zone_df, aes(Component_1, Component_2, color = Zone)) +
  geom_point(size = 1.7, alpha = 0.95) +
  scale_color_manual(values = zone_colors, drop = FALSE) +
  theme_void(base_size = 11) +
  ggtitle("B") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    legend.position = "bottom"
  )


# marker panels
panel_B_sox2  <- make_marker_panel(cds, "SOX2", "SOX2", plot_df)
panel_B_eomes <- make_marker_panel(cds, "EOMES", "EOMES", plot_df)
panel_B_myt1l <- make_marker_panel(cds, "MYT1L", "MYT1L", plot_df)


# combining plots
bottom_row <- arrangeGrob(
  panel_B_zone, panel_B_sox2, panel_B_eomes, panel_B_myt1l,
  ncol = 4
)

combined_fig <- arrangeGrob(
  panel_A,
  bottom_row,
  nrow = 2,
  heights = c(2.3, 1)
)


# showing the figure
grid.newpage()
grid.draw(combined_fig)


# saving files
ggsave("Figure2A_panel.png", panel_A, width = 7.5, height = 5.5, dpi = 300)
ggsave("Figure2A_panel.pdf", panel_A, width = 7.5, height = 5.5)

png("Figure2AB_combined.png", width = 2800, height = 1900, res = 260)
grid.draw(combined_fig)
dev.off()

pdf("Figure2AB_combined.pdf", width = 11.5, height = 8)
grid.draw(combined_fig)
dev.off()

write.csv(plot_df, "Figure2_monocle_coordinates.csv", row.names = FALSE)
saveRDS(cds, "Figure2_monocle_object.rds")

cat("\n===== DONE =====\n")
cat("Saved files:\n")
print(list.files(pattern = "Figure2"))
cat("\nWorking directory:\n")
print(getwd())