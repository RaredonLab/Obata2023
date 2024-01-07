# Set WD
setwd("/Users/msbr/T5/Pneumonectomy_For_Tuft_Project")

# Load packages
require(Seurat)

# Load data
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2022-08-17/pneum.project.data.sub.2022-08-18.Robj")

# Inspect data
table(Idents(data.sub))

# Pull pneumonectomy samples
SOI <- c('P0-f','P0-m','P3-f','P3-m','P7-f','P7-m','P14-f','P14-m')
pneum <- subset(data.sub,idents = SOI)
pneum <- NormalizeData(pneum)
rm(data.sub)
gc()

# Inspect QC (first time bringing all three together after initial filtration)
table(Idents(pneum))
pneum <- ScaleData(pneum)
pneum <- FindVariableFeatures(pneum)

pneum <- RunPCA(pneum,npcs = 100)
pdf(file='pneum.PCs.pdf',width=10,height=8)
ElbowPlot(pneum,ndims = 100)
PCHeatmap(pneum,cells=200,balanced=T,dims=1:9)
PCHeatmap(pneum,cells=200,balanced=T,dims=10:18)
PCHeatmap(pneum,cells=200,balanced=T,dims=19:27)
PCHeatmap(pneum,cells=200,balanced=T,dims=28:36)
PCHeatmap(pneum,cells=200,balanced=T,dims=37:45)
PCHeatmap(pneum,cells=200,balanced=T,dims=46:54)
PCHeatmap(pneum,cells=200,balanced=T,dims=55:63)
PCHeatmap(pneum,cells=200,balanced=T,dims=64:72)
PCHeatmap(pneum,cells=200,balanced=T,dims=73:81)
PCHeatmap(pneum,cells=200,balanced=T,dims=82:90)
PCHeatmap(pneum,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
pneum <- RunUMAP(pneum, reduction = "pca", dims = 1:40)
pneum <- FindNeighbors(pneum, reduction = "pca", dims = 1:40,force.recalc=T)
pneum <- FindClusters(pneum, resolution = 0.6)
p1 <- DimPlot(pneum, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(pneum, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p1 + p2

# Check QC
p2
FeaturePlot(pneum,'nFeature_RNA',label=T)
VlnPlot(pneum,'nFeature_RNA')
VlnPlot(pneum,'nCount_RNA',log=T)
VlnPlot(pneum,'nFeature_RNA',split.by = 'Sample')
DimPlot(pneum,cells = WhichCells(pneum,idents='12'))

# Remove cluster 0
pneum <- subset(pneum,idents = '0',invert=T)
table(Idents(pneum))

# Re-embed and cluster
pneum <- ScaleData(pneum)
pneum <- FindVariableFeatures(pneum)
pneum <- RunPCA(pneum,npcs = 100)
pdf(file='pneum.sub.PCs.pdf',width=10,height=8)
ElbowPlot(pneum,ndims = 100)
PCHeatmap(pneum,cells=200,balanced=T,dims=1:9)
PCHeatmap(pneum,cells=200,balanced=T,dims=10:18)
PCHeatmap(pneum,cells=200,balanced=T,dims=19:27)
PCHeatmap(pneum,cells=200,balanced=T,dims=28:36)
PCHeatmap(pneum,cells=200,balanced=T,dims=37:45)
PCHeatmap(pneum,cells=200,balanced=T,dims=46:54)
PCHeatmap(pneum,cells=200,balanced=T,dims=55:63)
PCHeatmap(pneum,cells=200,balanced=T,dims=64:72)
PCHeatmap(pneum,cells=200,balanced=T,dims=73:81)
PCHeatmap(pneum,cells=200,balanced=T,dims=82:90)
PCHeatmap(pneum,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
pneum <- RunUMAP(pneum, reduction = "pca", dims = 1:40)
pneum <- FindNeighbors(pneum, reduction = "pca", dims = 1:40,force.recalc=T)
pneum <- FindClusters(pneum, resolution = 0.2)
p1 <- DimPlot(pneum, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(pneum, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p1 + p2

# Find Markers Initial
pneum.mark <- FindAllMarkers(pneum,only.pos = T)
pneum.mark$ratio <- pneum.mark$pct.1/pneum.mark$pct.2
pneum.mark$power <- pneum.mark$ratio*pneum.mark$avg_log2FC
View(pneum.mark)

# Conclusion: looks clean enough to hand over to Allie for in-depth analysis
save(pneum,file = 'pneum.data.cleaned.2023-03-01.Robj')


