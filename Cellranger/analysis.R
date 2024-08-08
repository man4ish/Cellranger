# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)  # For combining plots

# Define file paths
# Update these file paths to the location of your Cell Ranger output files
cellranger_directory <- "path/to/CellRanger/output/directory"
output_directory <- "path/to/output/directory"

# Load data into Seurat
# Assuming Cell Ranger output is in the 'filtered_feature_bc_matrix' directory
seurat_data <- Read10X(data.dir = file.path(cellranger_directory, "filtered_feature_bc_matrix"))

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = seurat_data, project = "SingleCellRNASeq")

# Quality control
# Add metadata for quality control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Plot QC metrics
qc_plots <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(output_directory, "QC_plots.pdf"), plot = qc_plots)

# Filter cells based on QC metrics
# Customize thresholds based on your specific dataset
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Identify highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_obj <- ScaleData(seurat_obj, features = all_of(VariableFeatures(seurat_obj)))

# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Visualize PCA
pca_plot <- DimPlot(seurat_obj, reduction = "pca") + ggtitle("PCA Plot")
ggsave(file.path(output_directory, "PCA_plot.pdf"), plot = pca_plot)

# Cluster the cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)  # Adjust dims based on elbow plot
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)  # Adjust resolution as needed

# Perform UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Visualize UMAP
umap_plot <- DimPlot(seurat_obj, reduction = "umap") + ggtitle("UMAP Plot")
ggsave(file.path(output_directory, "UMAP_plot.pdf"), plot = umap_plot)

# Save Seurat object
saveRDS(seurat_obj, file = file.path(output_directory, "seurat_obj.rds"))

# Optional: Identify cluster markers
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = file.path(output_directory, "cluster_markers.csv"))

# Optional: Heatmap of top markers
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap_plot <- DoHeatmap(seurat_obj, features = top_markers$gene) + ggtitle("Heatmap of Top Markers")
ggsave(file.path(output_directory, "Heatmap_top_markers.pdf"), plot = heatmap_plot)

# Optional: Differential expression analysis between clusters
# Example comparing cluster 1 and 2
diff_exp <- FindMarkers(seurat_obj, ident.1 = 1, ident.2 = 2)
write.csv(diff_exp, file = file.path(output_directory, "diff_exp_cluster_1_vs_2.csv"))

