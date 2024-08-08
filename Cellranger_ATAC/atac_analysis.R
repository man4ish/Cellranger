# Load necessary libraries
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork"))
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Define paths to Cell Ranger ATAC-seq output files
output_dir <- "/path/to/cellranger/atac/output"  # Update this path
peak_file <- file.path(output_dir, "atac_fragments.tsv.gz")
barcode_file <- file.path(output_dir, "barcodes.tsv.gz")
metadata_file <- file.path(output_dir, "metadata.tsv.gz")

# Load the data
# Read fragments file
fragments <- Read10X_h5(file.path(output_dir, "atac_peak_matrix.h5"))

# Read metadata (if available)
metadata <- read.csv(metadata_file, sep = "\t", header = TRUE, row.names = 1)

# Create a Seurat object
atac_seurat <- CreateSeuratObject(
  counts = fragments,
  assay = "ATAC",
  meta.data = metadata
)

# Perform quality control
# Add a new column to metadata with the total number of fragments per cell
atac_seurat$total_fragments <- Matrix::colSums(atac_seurat@assays$ATAC@counts)

# Filter cells based on quality control metrics (customize thresholds as needed)
atac_seurat <- subset(atac_seurat, subset = total_fragments > 1000 & total_fragments < 50000)

# Normalize the data
atac_seurat <- NormalizeData(atac_seurat, assay = "ATAC", normalization.method = "LogNormalize")

# Perform dimensionality reduction
atac_seurat <- RunPCA(atac_seurat, assay = "ATAC", npcs = 30)
atac_seurat <- RunUMAP(atac_seurat, dims = 1:30)

# Find clusters
atac_seurat <- FindNeighbors(atac_seurat, dims = 1:30)
atac_seurat <- FindClusters(atac_seurat, resolution = 0.5)

# Plot the results
# UMAP plot
umap_plot <- DimPlot(atac_seurat, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP of Single-Cell ATAC-seq Data")

# Feature plot (example for a specific peak)
# Replace 'peak_name' with a specific peak of interest
feature_plot <- FeaturePlot(atac_seurat, features = "peak_name") +
  ggtitle("Feature Plot for Peak")

# Combine plots
combined_plot <- umap_plot + feature_plot

# Save the plots
ggsave("umap_plot.png", plot = umap_plot)
ggsave("feature_plot.png", plot = feature_plot)
ggsave("combined_plot.png", plot = combined_plot)

# Save the Seurat object for future use
saveRDS(atac_seurat, file = "atac_seurat.rds")

# Print completion message
print("Analysis completed. Results saved to files.")

