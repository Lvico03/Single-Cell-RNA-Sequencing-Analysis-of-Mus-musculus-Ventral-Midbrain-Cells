
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(patchwork)

# Load the sparse count matrix (will create an object called 'sm')
load("SRA670243_SRS3078084.sparse.RData")

# Clean gene names
sm@Dimnames[[1]] <- sub("_[^_]*$", "", sm@Dimnames[[1]])
genenames <- rownames(sm)
rownames(sm) <- make.names(genenames, unique = TRUE)
ncol(sm)
# Create Seurat object
midbrain <- CreateSeuratObject(
  counts = sm,
  project = "ventral_midbrain",
  min.cells = 3,
  min.features = 200
)
midbrain

# Cell Quality Control 

# Check mitochondrial gene names in your data, it uses the "." as separator instead of the dash
grep("^mt\\.", rownames(midbrain), value = TRUE)

# Calculate % mitochondrial reads per cell
midbrain[["percent.mt"]] <- PercentageFeatureSet(midbrain, pattern = "^mt\\.")

# Calculate % ribosomal protein reads per cell
midbrain[["percent.rbp"]] <- PercentageFeatureSet(midbrain, pattern = "^Rp[ls]")

# Check the metadata table - one row per cell, should now show percent.mt and percent.rbp columns
head(midbrain@meta.data, 5)

# Visualize all 4 QC metrics
VlnPlot(midbrain, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4)

# Violin plots without the individual dots
VlnPlot(midbrain, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4,
        pt.size = 0)
############################## VIOLIN PLOT QC ########################################

# Pull the metadata into a data frame
qc <- midbrain@meta.data

# A helper function to build a nice violin plot for one metric
make_violin <- function(data, feature, ylab, fill) {
  ggplot(data, aes(x = orig.ident, y = .data[[feature]])) +
    geom_violin(fill = fill, colour = "grey30", linewidth = 0.4,
                scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.08, outlier.shape = NA,
                 fill = "white", colour = "grey30", linewidth = 0.4) +
    labs(title = feature, x = NULL, y = ylab) +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 11)
    )
}

# Build one plot per metric
p1 <- make_violin(qc, "nFeature_RNA", "Genes per cell",        "#E8998D")
p2 <- make_violin(qc, "nCount_RNA",   "Counts per cell",       "#E8998D")
p3 <- make_violin(qc, "percent.mt",   "% Mitochondrial reads", "#E8998D")
p4 <- make_violin(qc, "percent.rbp",  "% Ribosomal reads",     "#E8998D")

# Combine side by side
(p1 | p2 | p3 | p4) +
  plot_annotation(title = "Quality control metrics",
                  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))
ggsave("qc_violins.png", width = 11, height = 4, dpi = 300)
###############################################################################################################

# Scatter plots to visualize relationships between QC metrics
# These help spot outliers - e.g. cells with high counts but low genes (doublets) or cells with high mt% (dying cells)
plot1 <- FeatureScatter(midbrain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(midbrain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Tresholds as suggested in the Seurat vignette: 200 < nFeature_RNA < 2500 and percent.mt < 5
midbrain <- subset(midbrain, 
                   subset = nFeature_RNA > 200 & 
                     nFeature_RNA < 2500 & 
                     percent.mt < 5)

midbrain

# Normalize the data
midbrain <- NormalizeData(midbrain, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)

# Normalized and not normalized data example
midbrain@assays
midbrain@assays$RNA 

# Raw counts for a specific gene (first 10 cells)
head(midbrain[["RNA"]]$counts["Rpl10",], 10)

# Normalized counts for the same gene
head(midbrain[["RNA"]]$data["Rpl10",], 10)

# Take a look of the gene with the highest mean expression among oure cells
apply(midbrain[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)

# Plot distribution of Malat1 and Gapdh across cells (housekeeping genes)
VlnPlot(midbrain, features = c("Malat1", "Gapdh"))

# Seurat contains a pre-computed list of cell cycle genes (human format - UPPERCASE)
cc.genes.updated.2019

###############à CONVERSION UPPERCASE FROM HUMAN TO MOUSE ############################à
# The built-in list is human (uppercase). Convert to mouse format:
# capitalize only the first letter, lowercase the rest.
# str_to_title from stringr does this, e.g. MCM5 -> Mcm5
library(stringr)

s.genes.mouse <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes.mouse <- str_to_title(cc.genes.updated.2019$g2m.genes)
head(s.genes.mouse)
head(g2m.genes.mouse)

# Score each cell for cell cycle phase
# This assigns each cell an S score, a G2M score, and a predicted Phase
midbrain <- CellCycleScoring(midbrain, 
                             s.features = s.genes.mouse, 
                             g2m.features = g2m.genes.mouse, 
                             set.ident = TRUE)

# Look at the metadata - now should include S.Score, G2M.Score, and Phase
head(midbrain[[]])

# Count cells in each predicted phase
table(midbrain$Phase)

##########################################################################################

# Find the Most Variable Feature, chooses the 2000 genes with hte highest variance. 
midbrain <- FindVariableFeatures(midbrain, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(midbrain), 10)
top10

# Plot variable features with labels
plot1 <- VariableFeaturePlot(midbrain)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

# Scaling all the genes 
# >0 high expression,0 average expression, <0 under expression/ no expression
all.genes <- rownames(midbrain)
midbrain <- ScaleData(midbrain, features = all.genes)

# Confirm scale.data layer now exists (three layers: counts, data, scale.data)
midbrain@assays$RNA

# Retrieve scaled normalized counts for a specific gene (first 20 cells)
head(midbrain[["RNA"]]$scale.data["Lyz2",], 20)

# PCA! Performed only on the 2000 variable genes, uses the scaled expression values.
midbrain <- RunPCA(midbrain, features = VariableFeatures(object = midbrain))

# Examine the top genes (loadings) for the first 5 PCs
# Positive loadings = "driver" genes for that component
print(midbrain[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize loadings for the first 2 PCs
VizDimLoadings(midbrain, dims = 1:2, reduction = "pca")
ggsave("pca_loadings.png", width = 10, height = 6, dpi = 300)

# Plot cells projected onto the first 2 PCs, colored by cell cycle phase
DimPlot(midbrain, reduction = "pca")

# Elbow plot - shows variance explained by each PC
# We're looking for the "elbow" where the curve flattens out


# Extract the standard deviation of each PC
sdev <- midbrain[["pca"]]@stdev[1:30]
df <- data.frame(PC = 1:30, stdev = sdev)

ggplot(df, aes(x = PC, y = stdev)) +
  geom_point(color = "#E8998D", size = 2.5) +
  theme_classic(base_size = 13) +
  labs(x = "Principal Component", y = "Standard Deviation")

ggsave("elbow_plot.png", width = 8, height = 5, dpi = 300)

# Also compute the 75% variance rule from the vignette
pc.touse <- (midbrain[["pca"]]@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse >= 0.75))
pc.touse

# Variance explained by each PC (variance = stdev squared)
pc.var <- (midbrain[["pca"]]@stdev)^2

# Convert to proportion of total variance
pc.var.prop <- pc.var / sum(pc.var)

# Cumulative variance explained
cumsum(pc.var.prop)[1:20]   # look at the first 20

# Specifically at 15 PCs
cumsum(pc.var.prop)[15]

#Clustering
# Build kNN graph with 10 PCs
midbrain10 <- FindNeighbors(midbrain, dims = 1:10)
midbrain10_05 <- FindClusters(midbrain10, resolution = 0.5)

head(Idents(midbrain10_05), 5)
table(midbrain10_05$seurat_clusters)

# Build kNN graph with 15 PCs
midbrain15 <- FindNeighbors(midbrain, dims = 1:15)
midbrain15_05 <- FindClusters(midbrain15, resolution = 0.5)

head(Idents(midbrain15_05), 5)
table(midbrain15_05$seurat_clusters)

# tSNE for 10 PCs (t-SNE: t-Distributed Stochastic Neighbor Embedding)
midbrain10_05 <- RunTSNE(midbrain10_05, dims = 1:10)
DimPlot(midbrain10_05, reduction = "tsne", label = TRUE)

# tSNE for 15 PCs
midbrain15_05 <- RunTSNE(midbrain15_05, dims = 1:15)
DimPlot(midbrain15_05, reduction = "tsne", label = TRUE)

# UMAP for both solutions
midbrain10_05 <- RunUMAP(midbrain10_05, dims = 1:10)
DimPlot(midbrain10_05, reduction = "umap", label = FALSE)
ggsave("umap_10PC.png", width = 8, height = 6, dpi = 300)

midbrain15_05 <- RunUMAP(midbrain15_05, dims = 1:15)
DimPlot(midbrain15_05, reduction = "umap", label = FALSE)
ggsave("umap_15PC.png", width = 8, height = 6, dpi = 300)

# Cells per cluster - 10 PC solution
table(midbrain10_05$seurat_clusters)

# Cells per cluster - 15 PC solution
table(midbrain15_05$seurat_clusters)

# Check quality parameters across clusters
VlnPlot(midbrain15_05, features = "nCount_RNA")
VlnPlot(midbrain15_05, features = "nFeature_RNA")
VlnPlot(midbrain15_05, features = "percent.mt")
VlnPlot(midbrain15_05, features = "percent.rbp")

# Check cell cycle phase distribution across clusters
midbrain15_05@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters, y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

# Find marker genes for all clusters at once
# only.pos = TRUE: only return overexpressed genes (not underexpressed)
# min.pct = 0.25: gene must be expressed in at least 25% of cells in the cluster
# logfc.threshold = 0.25: minimum fold change threshold
midbrain.markers <- FindAllMarkers(midbrain15_05, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, 
                                   logfc.threshold = 0.25)

# Print all 90 rows to see top 5 markers for all 18 clusters
midbrain.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n = 90)

# First FeaturePlot - visualize known marker genes on the UMAP
# Let's pick key markers for the cell types we identified
FeaturePlot(midbrain15_05, reduction = "umap",
            features = c("Aqp4",      # Astrocytes
                         "Opalin",     # Mature Oligodendrocytes
                         "Pdgfra",     # OPCs
                         "Tmem119",    # Microglia
                         "Cldn5",      # Endothelial cells
                         "Acta2",      # Smooth muscle cells
                         "Dcn",        # Fibroblasts
                         "Stmn2",      # Neurons
                         "Stoml3"))    # Ependymal cells
ggsave("featureplot_markers.png", width = 12, height = 10, dpi = 300)

# Get top 10 markers per cluster for the heatmap
midbrain.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


# Generate heatmap (yellow = high expression, purple = low)
DoHeatmap(midbrain15_05, features = top10$gene) + NoLegend()
ggsave("marker_heatmap.png", width = 14, height = 12, dpi = 300)



# INVESTIGATE CLUSTER 7
cluster7.markers <- FindMarkers(midbrain15_05, ident.1 = 7, 
                                min.pct = 0.25, test.use = "wilcox")
cluster7.markers <- cluster7.markers[order(-cluster7.markers$avg_log2FC),]
head(cluster7.markers, n = 15)

# Where do cluster 7 markers light up?
FeaturePlot(midbrain15_05, reduction = "umap",
            features = c("Ctps", "Tmem141", "Onecut2"))

# INVESTIGATE CLUSTER 12
cluster12.markers <- FindMarkers(midbrain15_05, ident.1 = 12, 
                                 min.pct = 0.25, test.use = "wilcox")
cluster12.markers <- cluster12.markers[order(-cluster12.markers$avg_log2FC),]
head(cluster12.markers, n = 15)

# Highlight clusters 7 and 12 on UMAP to see their position
DimPlot(midbrain15_05, reduction = "umap", 
        cells.highlight = list(
          cluster7 = WhichCells(midbrain15_05, idents = 7),
          cluster12 = WhichCells(midbrain15_05, idents = 12)),
        cols.highlight = c("red", "blue"),
        sizes.highlight = 1) +
  ggtitle("Clusters 7 (red) and 12 (blue)")

# OLIGODENDROCYTES CLUSTERING 

# Compare clusters 1 vs 2 to see if they are truly different
cluster1vs2.markers <- FindMarkers(midbrain15_05, 
                                   ident.1 = 1, ident.2 = 2,
                                   min.pct = 0.25, test.use = "wilcox")
cluster1vs2.markers <- cluster1vs2.markers[order(-cluster1vs2.markers$avg_log2FC),]
head(cluster1vs2.markers, n = 10)

# Now flip it - what is overexpressed in cluster 2 vs cluster 1?
cluster2vs1.markers <- FindMarkers(midbrain15_05, 
                                   ident.1 = 2, ident.2 = 1,
                                   min.pct = 0.25, test.use = "wilcox")
cluster2vs1.markers <- cluster2vs1.markers[order(-cluster2vs1.markers$avg_log2FC),]
head(cluster2vs1.markers, n = 10)

# Compare cluster 13 vs clusters 1 and 2, oligodendrocites hierarchy
cluster13vs12.markers <- FindMarkers(midbrain15_05,
                                     ident.1 = 13, ident.2 = c(1, 2),
                                     min.pct = 0.25, test.use = "wilcox")
cluster13vs12.markers <- cluster13vs12.markers[order(-cluster13vs12.markers$avg_log2FC),]
head(cluster13vs12.markers, n = 10)

#NEURON CLUSTERING 

# Compare neuron clusters 5 vs 15
cluster5vs15.markers <- FindMarkers(midbrain15_05,
                                    ident.1 = 5, ident.2 = 15,
                                    min.pct = 0.25, test.use = "wilcox")
cluster5vs15.markers <- cluster5vs15.markers[order(-cluster5vs15.markers$avg_log2FC),]
head(cluster5vs15.markers, n = 10)

# And vice versa
cluster15vs5.markers <- FindMarkers(midbrain15_05,
                                    ident.1 = 15, ident.2 = 5,
                                    min.pct = 0.25, test.use = "wilcox")
cluster15vs5.markers <- cluster15vs5.markers[order(-cluster15vs5.markers$avg_log2FC),]
head(cluster15vs5.markers, n = 10)

# Also check key neuronal markers
FeaturePlot(midbrain15_05, reduction = "umap",
            features = c("Gabra1", "Snca", "Th"))

# Visualize neuron clusters on UMAP
DimPlot(midbrain15_05, reduction = "umap",
        cells.highlight = list(
          cluster5 = WhichCells(midbrain15_05, idents = 5),
          cluster15 = WhichCells(midbrain15_05, idents = 15)),
        cols.highlight = c("red", "blue"),
        sizes.highlight = 1) +
  ggtitle("Clusters 5 (red) and 15 (blue)")

# FINAL RELABALED UMAP 
# Rename clusters with cell type labels
midbrain15_05 <- RenameIdents(midbrain15_05,
                              "0" = "Astrocytes",
                              "1" = "Mature oligodendrocytes",
                              "2" = "Newly formed oligodendrocytes",
                              "3" = "OPCs",
                              "4" = "Microglia",
                              "5" = "GABAergic neurons",
                              "6" = "Endothelial cells",
                              "7" = "Oligodendrocyte lineage",
                              "8" = "Oligodendrocyte lineage",
                              "9" = "Committed OPCs",
                              "10" = "Pericytes",
                              "11" = "Smooth muscle cells",
                              "12" = "Unknown",
                              "13" = "Myelinating oligodendrocytes",
                              "14" = "Fibroblasts",
                              "15" = "Dopaminergic neurons",
                              "16" = "Macrophages",
                              "17" = "Ependymal cells")

# Final labeled UMAP
DimPlot(midbrain15_05, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Ventral Midbrain - Cell Types") +
  NoLegend()

ggsave("final_umap_labeled.png", width = 10, height = 8, dpi = 300)
# KEY MARKER DOT PLOT FINAL

# Dot plot showing key markers for each cell type
DotPlot(midbrain15_05, 
        features = c("Aqp4",      # Astrocytes
                     "Opalin",    # Mature oligodendrocytes
                     "Klk6",      # Newly formed oligodendrocytes
                     "Mbp",       # Myelinating oligodendrocytes
                     "Pdgfra",    # OPCs
                     "Gpr17",     # Committed OPCs
                     "Tmem119",   # Microglia
                     "Mrc1",      # Macrophages
                     "Cldn5",     # Endothelial cells
                     "Abcc9",     # Pericytes
                     "Acta2",     # Smooth muscle cells
                     "Dcn",       # Fibroblasts
                     "Gabra1",    # GABAergic neurons
                     "Snca",      # Dopaminergic neurons
                     "Stoml3"),   # Ependymal cells
        cols = c("lightgrey", "blue")) +
  RotatedAxis()

ggsave("dotplot_markers.png", width = 12, height = 7, dpi = 300)

DotPlot(midbrain15_05, 
        features = c("Aqp4",      # Astrocytes
                     "Opalin",    # Mature oligodendrocytes
                     "Klk6",      # Newly formed oligodendrocytes
                     "Mog",     # Myelinating oligodendrocytes
                     "Pdgfra",    # OPCs
                     "Gpr17",     # Committed OPCs
                     "Tmem119",   # Microglia
                     "Mrc1",      # Macrophages
                     "Cldn5",     # Endothelial cells
                     "Abcc9",     # Pericytes
                     "Acta2",     # Smooth muscle cells
                     "Dcn",       # Fibroblasts
                     "Gabra1",    # GABAergic neurons
                     "Snca",      # Dopaminergic neurons
                     "Stoml3"),   # Ependymal cells
        cols = c("lightgrey", "blue")) +
  RotatedAxis()






























