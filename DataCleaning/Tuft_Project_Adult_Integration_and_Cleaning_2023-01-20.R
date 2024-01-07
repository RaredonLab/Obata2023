# Set WD
setwd("/Users/msbr/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20")

# Set seed
set.seed(2)

# Packages
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(cowplot)
library(future)


# Load previously worked on data (see "Tuft_Project_Data_Cleaning_2022-08-17.R")
load("/Users/msbr/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2022-08-17/pneum.project.data.sub.2022-08-18.Robj")
gc()

# Subset to only adult samples ('Native Atlas'). Do not include trachea, per Taylor's advice
table(Idents(data.sub))
temp.sub <- subset(data.sub,idents = c('2-12','2-14','2-16','2-28','2-7','2-8','3-1','8-10','rat1A','rat811','fRat','mRat','P0-f','P0-m') ) # Do not include TE1 here, it is too different to integrate well
table(Idents(temp.sub))
rm(data.sub)
gc()

#### First Look at Data Together (without integrating) ####
rat.adult <- temp.sub
rm(temp.sub)
gc()
rat.adult <- NormalizeData(rat.adult)
rat.adult <- ScaleData(rat.adult)
rat.adult <- FindVariableFeatures(rat.adult)
rat.adult <- RunPCA(rat.adult, npcs = 100)
pdf(file='rat.adult.PCs.pdf',width=10,height=8)
ElbowPlot(rat.adult,ndims = 100)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=1:9)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=10:18)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=19:27)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=28:36)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=37:45)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=46:54)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=55:63)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=64:72)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=73:81)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=82:90)
PCHeatmap(rat.adult,cells=200,balanced=T,dims=91:99)
dev.off()
gc()

# Embed and cluster
rat.adult <- RunUMAP(rat.adult, reduction = "pca", dims = 1:40)
rat.adult <- FindNeighbors(rat.adult, reduction = "pca", dims = 1:40)
rat.adult <- FindClusters(rat.adult, resolution = 0.2)
DimPlot(rat.adult, reduction = "umap", label = TRUE, repel = TRUE)


# Visualization
png(file='rat.adult.UMAP.first.look.png',width=12,height=6,units = 'in',res=300)
p1 <- DimPlot(rat.adult, reduction = "umap", group.by = "Sample",shuffle = T)
p2 <- DimPlot(rat.adult, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()

##### Integration (by chemistry/tissue) #####
table(rat.adult$Sample)
table(rat.adult$orig.ident)
Idents(rat.adult) <- rat.adult$orig.ident

lung.drop <- subset(rat.adult,idents = c('2-12','2-14','2-16','2-28','2-7','2-8','3-1','8-10','rat1A','rat811'))
lung.v2 <- subset(rat.adult,idents = c('fRat','mRat'))
lung.v3 <- subset(rat.adult,idents = c('P0-f','P0-m'))

# List
lung.list <- list(lung.drop,lung.v2,lung.v3)

# Clean up
rm(rat.adult)
rm(lung.drop)
rm(lung.v2)
rm(lung.v3)
gc()

# Make smaller for integration
# for(i in 1:length(lung.list)){
#   lung.list[[i]]@assays$spliced <- NULL
#   lung.list[[i]]@assays$unspliced <- NULL
#   lung.list[[i]]@assays$GeneFull <- NULL
#   lung.list[[i]]@meta.data$nCount_GeneFull <- NULL
#   lung.list[[i]]@meta.data$nFeature_GeneFull <- NULL
#   lung.list[[i]]@meta.data$nCount_spliced <- NULL
#   lung.list[[i]]@meta.data$nFeature_spliced <- NULL
#   lung.list[[i]]@meta.data$nCount_unspliced <- NULL
#   lung.list[[i]]@meta.data$nFeature_unspliced <- NULL
# }
# gc()

# normalize and identify variable features for each dataset independently
lung.list <- lapply(X = lung.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst") 
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = lung.list)

# integration
gc()
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, anchor.features = features)
gc()
save(lung.anchors,file='lung.anchors.Robj')
gc()
# this command creates an 'integrated' data assay
lung.combined <- IntegrateData(anchorset = lung.anchors)
gc()

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(lung.combined) <- "integrated"
save(lung.combined,file='rat.adult.integrated.2023-01-20.Robj')

# Run the standard workflow for visualization and clustering
lung.combined <- ScaleData(lung.combined)
lung.combined <- RunPCA(lung.combined, npcs = 100)
pdf(file='lung.combined.PCs.pdf',width=10,height=8)
ElbowPlot(lung.combined,ndims = 100)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
lung.combined <- RunUMAP(lung.combined, reduction = "pca", dims = 1:45)
lung.combined <- FindNeighbors(lung.combined, reduction = "pca", dims = 1:45)
DefaultAssay(lung.combined) <- "integrated"
lung.combined <- FindClusters(lung.combined, resolution = 0.4)
DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)

# Check lymphatics
DefaultAssay(lung.combined) <- 'RNA'
FeaturePlot(lung.combined,'Prox1',label=T,order = T,pt.size = 1)
FeaturePlot(lung.combined,'Cdh5',label=T,order = T,pt.size = 1)

# Check BASCs and Tuft
FeaturePlot(lung.combined,'Epcam',label=T,order = F,pt.size = 1)
FeaturePlot(lung.combined,'Sox9',label=T,order = F,pt.size = 1)
FeaturePlot(lung.combined,'Lgr5',label=T,order = F,pt.size = 1)
FeaturePlot(lung.combined,'Dclk1',label=T,order = F,pt.size = 1)
FeaturePlot(lung.combined,'Sftpc',label=T,order = F,pt.size = 1)
FeaturePlot(lung.combined,'Scgb1a1',label=T,order = F,pt.size = 1)

# Visualization
png(file='rat.adult.combined.UMAP.first.look.png',width=18,height=6,units = 'in',res=300)
p1 <- DimPlot(lung.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p2 <- DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()


# Class annotations
DefaultAssay(lung.combined) <- 'RNA'
png(file='lung.combined.UMAP.class.png',width=10,height=10,units = 'in',res=300)
FeaturePlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,order = T)
dev.off()
png(file='lung.combined.VLN.class.png',width=20,height=6,units = 'in',res=300)
VlnPlot(lung.combined,pt.size=0,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 2)
dev.off()
png(file='lung.combined.VLN.QC.png',width=20,height=6,units = 'in',res=300)
VlnPlot(lung.combined,pt.size=0,log = T,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol = 2)
dev.off()
Idents(lung.combined) <- lung.combined$integrated_snn_res.0.4
lung.combined <- RenameIdents(lung.combined,
                              '0'='Endothelium',
                              '1'='Epithelium',
                              '2'='Epithelium',
                              '3'='Immune',
                              '4'='Immune',
                              '5'='Immune',
                              '6'='Immune',
                              '7'='LowInfo',
                              '8'='Immune',
                              '9'='Immune',
                              '10'='Mesenchyme',
                              '11'='Immune',
                              '12'='Epithelium',
                              '13'='Epithelium',
                              '14'='Mesenchyme',
                              '15'='Immune',
                              '16'='Immune',
                              '17'='Endothelium',
                              '18'='Immune',
                              '19'='Mesenchyme',
                              '20'='Immune',
                              '21'='Epithelium',
                              '22'='Endothelium',
                              '23'='Immune',
                              '24'='Immune',
                              '25'='Multiplet',
                              '26'='Multiplet')

# Check cluster 7 which looks like lowinfo
mark.7 <- FindMarkers(lung.combined,ident.1 = '7',min.pct = 0.5,only.pos = T)
pdf(file='LowInfo.cluster.7.on.first.integration.embedding.pdf',width = 14,height = 10)
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '7'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '7'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '7'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(lung.combined,features = c('percent.mt'))
VlnPlot(lung.combined,features = c('nFeature_RNA'))
VlnPlot(lung.combined,features = c('nCount_RNA'),log=T)
dev.off()

# Stash new class labels
lung.combined$CellClass <- Idents(lung.combined)

# Record landing point
Idents(lung.combined) <- lung.combined$CellClass
pdf(file='rat.adult.combined.UMAP.class.labeled.pdf',width=12,height=10)
DimPlot(lung.combined,group.by = 'CellClass')+ggtitle('Cell Class - Initial')
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Immune'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Epithelium'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Endothelium'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Mesenchyme'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'LowInfo'),features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA','percent.mt','perc.spliced')) 
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Multiplet'),features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
dev.off()

# save
save(lung.combined,file='rat.adult.combined.classed.for.Tomo.2023-01-21.Robj')
rm(lung.anchors)
rm(lung.list)
gc()
lung.combined.metadata <- lung.combined@meta.data
save(lung.combined.metadata,file='METADATA.rat.adult.combined.classed.for.Tomo.2023-01-21.Robj')

#### Class Separation ####
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/rat.adult.combined.classed.for.Tomo.2023-01-21.Robj")
lung.combined
table(lung.combined$orig.ident)
Idents(lung.combined) <- lung.combined$orig.ident
lung.combined <- RenameIdents(lung.combined,
                              '2-12'='DropSeq',   
                              '2-14'='DropSeq',   
                              '2-16'='DropSeq',   
                              '2-28'='DropSeq',    
                              '2-7'='DropSeq',    
                              '2-8'='DropSeq',    
                              '3-1'='DropSeq',   
                              '8-10'='DropSeq',   
                              'fRat'='10x_v2',   
                              'mRat'='10x_v2',   
                              'P0-f'='10x_v3',   
                              'P0-m'='10x_v3',  
                              'rat1A'='DropSeq', 
                              'rat811'='DropSeq')
lung.combined$Dataset <- Idents(lung.combined)
table(lung.combined$Dataset)

# Break out each cell class
Idents(lung.combined) <- lung.combined$CellClass
epi <- subset(lung.combined,idents = 'Epithelium')
endo <- subset(lung.combined,idents = 'Endothelium')
mes <- subset(lung.combined,idents = 'Mesenchyme')
imm <- subset(lung.combined,idents = 'Immune')
gc()

# Remove total object to save working space
#rm(lung.combined)
gc()

##### Epi Subclustering and Cleaning #####
DefaultAssay(epi) <- 'integrated' # Let's try using the old integration-space
epi <- ScaleData(epi)
epi<- RunPCA(epi, npcs = 100)
pdf(file='epi.PCs.pdf',width=10,height=8)
ElbowPlot(epi,ndims = 100)
PCHeatmap(epi,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
epi <- RunUMAP(epi, reduction = "pca", dims = 1:30)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:30)
DefaultAssay(epi) <- 'integrated'
epi <- FindClusters(epi, resolution = 0.4)

# Visualization
pdf(file='epi.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(epi, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi) <- 'RNA'
p4 <- FeaturePlot(epi,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(epi,features = c('Sox9','Sftpc','Ager','Dclk1',
                                   'Scgb1a1','Ccdc153','Lgr5','Krt5'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(epi,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(epi,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(epi,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

### Clean up ###
# Check cluster 8,10,11 which looks like multiplets
pdf(file='epi.cluster.8.10.11.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(epi,cells = WhichCells(epi,idents = c('8','10','11')))
FeaturePlot(epi,cells = WhichCells(epi),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi,idents = c('8','10','1')),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(epi,cells = WhichCells(epi,idents = c('8','10','11')),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(epi,features = c('percent.mt'))
VlnPlot(epi,features = c('nFeature_RNA'))
VlnPlot(epi,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
epi.clean <- subset(epi,idents = c('8','10','11'),invert = T)

# Re-Visualize
DefaultAssay(epi.clean) <- 'integrated' 
epi.clean <- ScaleData(epi.clean)
epi.clean<- RunPCA(epi.clean, npcs = 100)
pdf(file='epi.clean.PCs.pdf',width=10,height=8)
ElbowPlot(epi.clean,ndims = 100)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi.clean,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
DefaultAssay(epi.clean) <- 'integrated'
epi.clean <- RunUMAP(epi.clean, reduction = "pca", dims = 1:15)
epi.clean <- FindNeighbors(epi.clean, reduction = "pca", dims = 1:15)
epi.clean <- FindClusters(epi.clean, resolution = 0.4)

# Visualization
pdf(file='epi.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(epi.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi.clean) <- 'RNA'
p4 <- FeaturePlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(epi.clean,features = c('Sox9','Sftpc','Ager','Dclk1',
                                   'Scgb1a1','Ccdc153','Lgr5','Krt5'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(epi.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(epi.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()


# Identify Cell Types
DefaultAssay(epi.clean) <- 'RNA'
epi.clean.mark <- FindAllMarkers(epi.clean,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
epi.clean.mark$ratio <- epi.clean.mark$pct.1/epi.clean.mark$pct.2
save(epi.clean.mark,file='epi.clean.mark.Robj')
FeaturePlot(epi.clean,'Aqp5',label=T,label.size = 10)
epi.clean <- RenameIdents(epi.clean,
                            '0'='ATII',
                            '1'='ATII',
                            '2'='ATI',
                            '3'='ATI',
                            '4'='ATI',
                            '5'='Tuft',
                            '6'='Secretory',
                            '7'='Ciliated',
                            '8'='Tuft',
                            '9'='ATII-ATI',
                            '10'='BASC')

epi.clean$CellType <- Idents(epi.clean)

pdf(file='epi.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(epi.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi.clean) <- 'RNA'
p4 <- FeaturePlot(epi.clean,features = c('Ager','Sftpc','Ccdc153','Scgb1a1'),label = T,ncol = 3)
p5 <- VlnPlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
save(epi.clean,file='epi.clean.labeled.2023-01-29.Robj')
epi.clean.metadata <- epi.clean@meta.data
save(epi.clean.metadata,file = 'epi.clean.metadata.2023-01-29.Robj')
gc()
rm(epi)
rm(epi.clean.mark)
rm(epi.clean)
rm(epi.clean.metadata)
gc()

##### Endo Subclustering and Cleaning #####
DefaultAssay(endo) <- 'integrated' # Let's try using the old integration-space
endo <- ScaleData(endo)
endo<- RunPCA(endo, npcs = 100)
pdf(file='endo.PCs.pdf',width=10,height=8)
ElbowPlot(endo,ndims = 100)
PCHeatmap(endo,cells=200,balanced=T,dims=1:9)
PCHeatmap(endo,cells=200,balanced=T,dims=10:18)
PCHeatmap(endo,cells=200,balanced=T,dims=19:27)
PCHeatmap(endo,cells=200,balanced=T,dims=28:36)
PCHeatmap(endo,cells=200,balanced=T,dims=37:45)
PCHeatmap(endo,cells=200,balanced=T,dims=46:54)
PCHeatmap(endo,cells=200,balanced=T,dims=55:63)
PCHeatmap(endo,cells=200,balanced=T,dims=64:72)
PCHeatmap(endo,cells=200,balanced=T,dims=73:81)
PCHeatmap(endo,cells=200,balanced=T,dims=82:90)
PCHeatmap(endo,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
endo <- RunUMAP(endo, reduction = "pca", dims = 1:15)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:15)
endo <- FindClusters(endo, resolution = 0.2)

# Visualization
pdf(file='endo.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(endo, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo) <- 'RNA'
p4 <- FeaturePlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(endo,features = c('Prx','Vwf','Apln','Aplnr',
                                         'Ca4','Gja5','Efnb2','Nr2f2'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(endo,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(endo,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

### Clean up ###
# Check cluster 7 which looks like multiplets
pdf(file='endo.cluster.7.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(endo,cells = WhichCells(endo,idents = c('7')))
FeaturePlot(endo,cells = WhichCells(endo),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(endo,cells = WhichCells(endo),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
VlnPlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(endo,features = c('percent.mt'))
VlnPlot(endo,features = c('nFeature_RNA'))
VlnPlot(endo,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
endo.clean <- subset(endo,idents = c('7'),invert = T)

# Re-Visualize
DefaultAssay(endo.clean) <- 'integrated' 
endo.clean <- ScaleData(endo.clean)
endo.clean<- RunPCA(endo.clean, npcs = 100)
pdf(file='endo.clean.PCs.pdf',width=10,height=8)
ElbowPlot(endo.clean,ndims = 100)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=1:9)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=10:18)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=19:27)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=28:36)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=37:45)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=46:54)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=55:63)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=64:72)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=73:81)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=82:90)
PCHeatmap(endo.clean,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
DefaultAssay(endo.clean) <- 'integrated' 
endo.clean <- RunUMAP(endo.clean, reduction = "pca", dims = 1:13)
endo.clean <- FindNeighbors(endo.clean, reduction = "pca", dims = 1:13)
endo.clean <- FindClusters(endo.clean, resolution = 0.3)

# Visualization
pdf(file='endo.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(endo.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo.clean) <- 'RNA'
p4 <- FeaturePlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(endo.clean,features = c('Prx','Vcam1','Apln','Aplnr',
                                    'Ca4','Gja5','Efnb2','Nr2f2'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(endo.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(endo.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()


# Identify Cell Types
DefaultAssay(endo.clean) <- 'RNA'
endo.clean.mark <- FindAllMarkers(endo.clean,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
endo.clean.mark$ratio <- endo.clean.mark$pct.1/endo.clean.mark$pct.2
save(endo.clean.mark,file='endo.clean.mark.Robj')
FeaturePlot(endo.clean,'Nr2f2',label=T,label.size = 10)
endo.clean <- RenameIdents(endo.clean,
                           '0'='gCap',
                           '1'='gCap',
                           '2'='Arterial',
                           '3'='aCap',
                           '4'='Venous',
                           '5'='gCap',
                           '6'='Lymphatic')
endo.clean$CellType <- Idents(endo.clean)

pdf(file='endo.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(endo.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo.clean) <- 'RNA'
p4 <- FeaturePlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(endo.clean,features = c('Prx','Vcam1','Apln','Aplnr',
                                          'Ca4','Gja5','Efnb2','Nr2f2'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(endo.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(endo.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

# Save and Clean up
save(endo.clean,file='endo.clean.labeled.2023-01-29.Robj')
endo.clean.metadata <- endo.clean@meta.data
save(endo.clean.metadata,file = 'endo.clean.metadata.2023-01-29.Robj')
gc()
rm(endo)
rm(endo.clean.mark)
rm(endo.clean)
rm(endo.clean.metadata)
gc()






##### Mes Subclustering and Cleaning #####
DefaultAssay(mes) <- 'integrated' # Let's try using the old integration-space
mes <- ScaleData(mes)
mes<- RunPCA(mes, npcs = 100)
pdf(file='mes.PCs.pdf',width=10,height=8)
ElbowPlot(mes,ndims = 100)
PCHeatmap(mes,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
mes <- RunUMAP(mes, reduction = "pca", dims = 1:15)
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:15)
DefaultAssay(mes) <- 'integrated'
mes <- FindClusters(mes, resolution = 0.4)

# Visualization
pdf(file='mes.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(mes, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes) <- 'RNA'
p4 <- FeaturePlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(mes,features = c('Lgr5','Lgr6','Col13a1','Col14a1',
                                          'Acta2','Gucy1a1','Wnt11','Pi16'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(mes,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(mes,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

### Clean up ###
# Check cluster 8 & 9 which looks like multiplets
pdf(file='mes.cluster.8.9.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(mes,cells = WhichCells(mes,idents = c('8','9')))
FeaturePlot(mes,cells = WhichCells(mes),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(mes,cells = WhichCells(mes),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
VlnPlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(mes,features = c('percent.mt'))
VlnPlot(mes,features = c('nFeature_RNA'))
VlnPlot(mes,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
mes.clean <- subset(mes,idents = c('8','9'),invert = T)

# Re-Visualize
DefaultAssay(mes.clean) <- 'integrated' 
mes.clean <- ScaleData(mes.clean)
mes.clean<- RunPCA(mes.clean, npcs = 100)
pdf(file='mes.clean.PCs.pdf',width=10,height=8)
ElbowPlot(mes.clean,ndims = 100)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes.clean,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
mes.clean <- RunUMAP(mes.clean, reduction = "pca", dims = 1:14)
mes.clean <- FindNeighbors(mes.clean, reduction = "pca", dims = 1:14)
DefaultAssay(mes.clean) <- 'integrated'
mes.clean <- FindClusters(mes.clean, resolution = 0.6)

# Visualization
pdf(file='mes.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(mes.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes.clean) <- 'RNA'
p4 <- FeaturePlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(mes.clean,features = c('Lgr5','Lgr6','Col13a1','Col14a1',
                                   'Acta2','Gucy1a1','Wnt11','Pi16'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(mes.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(mes.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

# Identify Cell Types
DefaultAssay(mes.clean) <- 'RNA'
mes.clean.mark <- FindAllMarkers(mes.clean,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
mes.clean.mark$ratio <- mes.clean.mark$pct.1/mes.clean.mark$pct.2
save(mes.clean.mark,file='mes.clean.mark.Robj')
FeaturePlot(mes.clean,'Epcam',label=T,label.size = 10,pt.size = 1,order=T)
mark.10 <- FindMarkers(mes.clean,ident.1 = '10',ident.2 = '4')
mark.10$ratio <- mark.10$pct.1/mark.10$pct.2
mark.8 <- FindMarkers(mes.clean,ident.1 = '8',ident.2 = '5')
mes.clean <- RenameIdents(mes.clean,
                          '0'='Col13a1_Fib',
                          '1'='Myofibroblasts',
                          '2'='Col13a1_Fib',
                          '3'='Col14a1_Fib',
                          '4'='Mesothelium',
                          '5'='SMCs',
                          '6'='Multiplet',
                          '7'='Lgr5_Fib',
                          '8'='LowInfo',
                          '9'='Pericytes',
                          '10'='Wfdc21_Mesothelium')
mes.clean$CellType <- Idents(mes.clean)

pdf(file='mes.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(mes.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes.clean) <- 'RNA'
p4 <- FeaturePlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(mes.clean,features = c('Lgr5','Lgr6','Col13a1','Col14a1',
                                         'Acta2','Gucy1a1','Wnt11','Pi16'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(mes.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(mes.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

# Save and Clean up
save(mes.clean,file='mes.clean.labeled.2023-01-29.Robj')
mes.clean.metadata <- mes.clean@meta.data
save(mes.clean.metadata,file = 'mes.clean.metadata.2023-01-29.Robj')
gc()
rm(mes)
rm(mes.clean.mark)
rm(mes.clean)
rm(mes.clean.metadata)
gc()

##### Imm Subclustering and Cleaning #####
DefaultAssay(imm) <- 'integrated' # Let's try using the old integration-space
imm <- ScaleData(imm)
imm<- RunPCA(imm, npcs = 100)
pdf(file='imm.PCs.pdf',width=10,height=8)
ElbowPlot(imm,ndims = 100)
PCHeatmap(imm,cells=200,balanced=T,dims=1:9)
PCHeatmap(imm,cells=200,balanced=T,dims=10:18)
PCHeatmap(imm,cells=200,balanced=T,dims=19:27)
PCHeatmap(imm,cells=200,balanced=T,dims=28:36)
PCHeatmap(imm,cells=200,balanced=T,dims=37:45)
PCHeatmap(imm,cells=200,balanced=T,dims=46:54)
PCHeatmap(imm,cells=200,balanced=T,dims=55:63)
PCHeatmap(imm,cells=200,balanced=T,dims=64:72)
PCHeatmap(imm,cells=200,balanced=T,dims=73:81)
PCHeatmap(imm,cells=200,balanced=T,dims=82:90)
PCHeatmap(imm,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
imm <- RunUMAP(imm, reduction = "pca", dims = 1:50)
imm <- FindNeighbors(imm, reduction = "pca", dims = 1:50)
DefaultAssay(imm) <- 'integrated'
imm <- FindClusters(imm, resolution = 0.6)

# Visualization
pdf(file='imm.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(imm, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(imm) <- 'RNA'
p4 <- FeaturePlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(imm,features = c('Cd3e','Cd79a','Mrc1','C1qc',
                                         'Naaa','Jchain','Nkg7','Ly6c'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(imm,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(imm,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

### Clean up ###
# Check cluster 15 & 16 & 17 & 19 which look like multiplets
pdf(file='imm.cluster.15.16.17.19.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(imm,cells = WhichCells(imm,idents = c('15','16','17','19')))
FeaturePlot(imm,cells = WhichCells(imm),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(imm,cells = WhichCells(imm,idents = c('15','16','17','19')),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,order=T)
FeaturePlot(imm,cells = WhichCells(imm),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
FeaturePlot(imm,cells = WhichCells(imm,idents = c('15','16','17','19')),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
VlnPlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(imm,features = c('percent.mt'))
VlnPlot(imm,features = c('nFeature_RNA'))
VlnPlot(imm,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
imm.clean <- subset(imm,idents = c('15','16','17','19'),invert = T)

# Re-Visualize
DefaultAssay(imm.clean) <- 'integrated' 
imm.clean <- ScaleData(imm.clean)
imm.clean<- RunPCA(imm.clean, npcs = 100)
pdf(file='imm.clean.PCs.pdf',width=10,height=8)
ElbowPlot(imm.clean,ndims = 100)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=1:9)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=10:18)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=19:27)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=28:36)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=37:45)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=46:54)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=55:63)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=64:72)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=73:81)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=82:90)
PCHeatmap(imm.clean,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
imm.clean <- RunUMAP(imm.clean, reduction = "pca", dims = 1:35)
imm.clean <- FindNeighbors(imm.clean, reduction = "pca", dims = 1:35)
DefaultAssay(imm.clean) <- 'integrated'
imm.clean <- FindClusters(imm.clean, resolution = 0.6)

# Visualization
pdf(file='imm.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(imm.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(imm.clean) <- 'RNA'
p4 <- FeaturePlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(imm.clean,features = c('Cd3e','Cd79a','Mrc1','C1qc',
                                   'Naaa','Jchain','Nkg7','Ly6c'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(imm.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(imm.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

# Identify Cell Types
DefaultAssay(imm.clean) <- 'RNA'
imm.clean.mark <- FindAllMarkers(imm.clean,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
imm.clean.mark$ratio <- imm.clean.mark$pct.1/imm.clean.mark$pct.2
save(imm.clean.mark,file='imm.clean.mark.Robj')
FeaturePlot(imm.clean,'Pclaf',label=T,label.size = 10,pt.size = 1,order=T)
imm.clean <- RenameIdents(imm.clean,
                          '0'='B',
                          '1'='Mac_Alv',
                          '2'='Monocytes',
                          '3'='T',
                          '4'='NK',
                          '5'='Mac_Alv',
                          '6'='Monocytes',
                          '7'='T',
                          '8'='Neutrophils',
                          '9'='Mac_Inter',
                          '10'='DC',
                          '11'='ILCs',
                          '12'='Killer_T_Prolif',
                          '13'='Mac_Prolif',
                          '14'='Siglech+',
                          '15'='Plasma',
                          '16'='Megakaryocytes',
                          '17'='Multiplet')
imm.clean$CellType <- Idents(imm.clean)

pdf(file='imm.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(imm.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(imm.clean) <- 'RNA'
p4 <- FeaturePlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 4,order = T,pt.size = 1)
p5 <- FeaturePlot(imm.clean,features = c('Cd3e','Cd79a','Mrc1','C1qc',
                                         'Naaa','Jchain','Nkg7','Ly6c'),label = T,ncol = 4,order = T,pt.size = 1)
p6 <- FeaturePlot(imm.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T,ncol = 4,order = T,pt.size = 1)
p7 <- VlnPlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),log=T,pt.size=0,ncol = 2)
p8 <- VlnPlot(imm.clean,features = c('nCount_RNA','nFeature_RNA','percent.mt'),log= T,pt.size=0)
p1 + p2 + p3
p4
p5
p6
p7
p8
dev.off()

# Save and Clean up
save(imm.clean,file='imm.clean.labeled.2023-01-29.Robj')
imm.clean.metadata <- imm.clean@meta.data
save(imm.clean.metadata,file = 'imm.clean.metadata.2023-01-29.Robj')
gc()
rm(imm)
rm(imm.clean.mark)
rm(imm.clean)
rm(imm.clean.metadata)
gc()

#### Merge Clean Classes Back Together ####
load('imm.clean.labeled.2023-01-29.Robj')
load('mes.clean.labeled.2023-01-29.Robj')
load('endo.clean.labeled.2023-01-29.Robj')
load('epi.clean.labeled.2023-01-29.Robj')
gc()

# Merge
merge.clean <- merge(endo.clean,list(epi.clean,imm.clean,mes.clean))
rm(mes.clean)
rm(epi.clean)
rm(imm.clean)
rm(endo.clean)
gc()


# Clean
merge.clean
table(merge.clean$Dataset)
table(merge.clean$CellClass)
table(merge.clean$CellType)
table(merge.clean$CellClass,merge.clean$Dataset)
table(merge.clean$CellType,merge.clean$Dataset)
table(merge.clean$CellType,merge.clean$CellClass)

merge.clean <- subset(merge.clean,idents = c('LowInfo','Multiplet'),invert=T)
table(merge.clean$CellType,merge.clean$Dataset)
table(merge.clean$CellType,merge.clean$CellClass)

# Re-Integrate
# List
lung.list <- SplitObject(merge.clean,split.by = 'Dataset')
gc()

# normalize and identify variable features for each dataset independently
lung.list <- lapply(X = lung.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst") 
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = lung.list)

# integration
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, anchor.features = features)
gc()

lung.combined <- IntegrateData(anchorset = lung.anchors)
gc()

# Visualize
DefaultAssay(lung.combined) <- "integrated"
lung.combined <- ScaleData(lung.combined)
lung.combined <- RunPCA(lung.combined, npcs = 100)
pdf(file='lung.clean.combined.PCs.pdf',width=10,height=8)
ElbowPlot(lung.combined,ndims = 100)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(lung.combined,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
lung.combined <- RunUMAP(lung.combined, reduction = "pca", dims = 1:45)
lung.combined <- FindNeighbors(lung.combined, reduction = "pca", dims = 1:45)
DefaultAssay(lung.combined) <- "integrated"
lung.combined <- FindClusters(lung.combined, resolution = 0.6)

# Visualization
png(file='lung.clean.combined.UMAP.first.look.png',width=20,height=18,units = 'in',res=400)
p1 <- DimPlot(lung.combined, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(lung.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)
p4 <- DimPlot(lung.combined, reduction = "umap",  group.by = "CellType",label = TRUE, repel = TRUE)
p1 + p2 + p3 + p4
dev.off()


# Save for later
rm(merge.clean)
rm(lung.list)
rm(lung.anchors)
gc()
save(lung.combined,file = 'lung.clean.combined.2023-01-29.Robj')

#### Re-Annotate Class with fully clean integration ####
# Load from save
load('lung.clean.combined.2023-01-29.Robj')

# Class annotations
DefaultAssay(lung.combined) <- 'RNA'
png(file='lung.clean.combined.UMAP.class.png',width=10,height=10,units = 'in',res=300)
FeaturePlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,order = T)
dev.off()
png(file='lung.clean.combined.VLN.class.png',width=20,height=6,units = 'in',res=300)
VlnPlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 2,pt.size = 0)
dev.off()

lung.combined <- RenameIdents(lung.combined,
                              '0'='Endothelium',
                              '1'='Immune',
                              '2'='Immune',
                              '3'='Immune',
                              '4'='Epithelium',
                              '5'='Epithelium',
                              '6'='Immune',
                              '7'='Immune',
                              '8'='Epithelium',
                              '9'='Mesenchyme',
                              '10'='Epithelium',
                              '11'='Epithelium',
                              '12'='Immune',
                              '13'='Immune',
                              '14'='Endothelium',
                              '15'='Epithelium',
                              '16'='Immune',
                              '17'='Mesenchyme',
                              '18'='Endothelium',
                              '19'='Immune',
                              '20'='Immune',
                              '21'='Immune',
                              '22'='Immune',
                              '23'='Mesenchyme',
                              '24'='Epithelium',
                              '25'='Immune',
                              '26'='Epithelium',
                              '27'='Endothelium',
                              '28'='Immune')


# Stash new class labels
lung.combined$CellClass_Final <- Idents(lung.combined)

# Record landing point
Idents(lung.combined) <- lung.combined$CellClass_Final
pdf(file='lung.clean.combined.UMAP.class.labeled.pdf',width=12,height=10)
DimPlot(lung.combined,group.by = 'CellClass')+ggtitle('Cell Class - Final')
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Immune'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Epithelium'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Endothelium'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Mesenchyme'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
dev.off()

# Save metadata
lung.combined.metadata <- lung.combined@meta.data
save(lung.combined.metadata,file='lung.clean.combined.classed.metadata.2023-01-29.Robj')

##### Epi Subclustering and Identification #####
epi <- subset(lung.combined,idents = 'Epithelium')
gc()
DefaultAssay(epi) <- 'integrated' # Let's try using the old integration-space
epi <- ScaleData(epi)
epi<- RunPCA(epi, npcs = 100)
pdf(file='epi.round2.PCs.pdf',width=10,height=8)
ElbowPlot(epi,ndims = 100)
PCHeatmap(epi,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
DefaultAssay(epi) <- 'integrated'
epi <- RunUMAP(epi, reduction = "pca", dims = 1:8)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:8)
epi <- FindClusters(epi, resolution = 0.2)

# Visualization
pdf(file='epi.round2.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(epi, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi) <- 'RNA'
p4 <- FeaturePlot(epi,features = c('Sftpc','Scgb1a1','Dclk1','Ccdc153','Pdpn','Sox9'),label = T,ncol = 3)
p5 <- VlnPlot(epi,features = c('Sftpc','Scgb1a1','Dclk1','Ccdc153','Pdpn','Sox9'))
p1 + p2 + p3
p4
p5
dev.off()

# Identify Cell Types
DimPlot(epi,group.by = 'CellType')
epi <- RenameIdents(epi,
                    '0'='ATII',
                    '1'='ATI',
                    '2'='ATI',
                    '3'='ATII',
                    '4'='Tuft',
                    '5'='Secretory',
                    '6'='Ciliated',
                    '7'='BASC',
                    '8'='ATII-ATI')
epi$CellType_Final <- Idents(epi)

pdf(file='epi.final.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(epi, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi) <- 'RNA'
p4 <- FeaturePlot(epi,features = c('Sftpc','Scgb1a1','Dclk1','Ccdc153','Pdpn','Sox9'),label = T,ncol = 3)
p1 + p2 + p3
p4
dev.off()

# Save Epi Metadata
save(epi,file = 'epi.final.2023-01-29.Robj')
epi.final.metadata <- epi@meta.data
save(epi.final.metadata,file = 'epi.final.metadata.2023-01-29.Robj')
gc()
rm(epi)
gc()

##### Endo Subclustering and Identification #####
endo <- subset(lung.combined,idents = 'Endothelium')
DefaultAssay(endo) <- 'integrated' # Let's try using the old integration-space
endo <- ScaleData(endo)
endo<- RunPCA(endo, npcs = 100)
pdf(file='endo.round2.PCs.pdf',width=10,height=8)
ElbowPlot(endo,ndims = 100)
PCHeatmap(endo,cells=200,balanced=T,dims=1:9)
PCHeatmap(endo,cells=200,balanced=T,dims=10:18)
PCHeatmap(endo,cells=200,balanced=T,dims=19:27)
PCHeatmap(endo,cells=200,balanced=T,dims=28:36)
PCHeatmap(endo,cells=200,balanced=T,dims=37:45)
PCHeatmap(endo,cells=200,balanced=T,dims=46:54)
PCHeatmap(endo,cells=200,balanced=T,dims=55:63)
PCHeatmap(endo,cells=200,balanced=T,dims=64:72)
PCHeatmap(endo,cells=200,balanced=T,dims=73:81)
PCHeatmap(endo,cells=200,balanced=T,dims=82:90)
PCHeatmap(endo,cells=200,balanced=T,dims=91:99) 
dev.off()

# Embed and cluster
DefaultAssay(endo) <- 'integrated'
endo <- RunUMAP(endo, reduction = "pca", dims = c(1,3,4:6,9:10))
endo <- FindNeighbors(endo, reduction = "pca", dims = c(1,3,4:6,9:10))
endo <- FindClusters(endo, resolution = 0.2)

# Visualization
pdf(file='endo.round2.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(endo, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo) <- 'RNA'
p4 <- FeaturePlot(endo,features = c('Prx','Vcam1','Apln','Aplnr',
                                   'Ca4','Gja5','Efnb2','Nr2f2'),label = T,ncol = 3)
p5 <- VlnPlot(endo,features = c('Prx','Vcam1','Apln','Aplnr',
                               'Ca4','Gja5','Efnb2','Nr2f2'))
p1 + p2 + p3
p4
p5
dev.off()

# Identify Cell Types
DimPlot(endo, reduction = "umap", group.by = "CellType",shuffle = T,label=T)
endo <- RenameIdents(endo,
                     '0'='gCap',
                     '1'='gCap',
                     '2'='Arterial',#Efnb2
                     '3'='Venous',  #Vcam1 
                     '4'='aCap',
                     '5'='gCap',
                     '6'='Lymphatic'
)
endo$CellType_Final <- Idents(endo)

pdf(file='endo.Round2.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(endo, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo) <- 'RNA'
p5 <- VlnPlot(endo,features = c('Efnb2','Vcam1','Prx','Plvap','Prox1'),pt.size=0.1)
p1 + p2 + p3
p5
dev.off()

# Save and Clean up
save(endo,file = 'endo.final.2023-01-29.Robj')
endo.metadata <- endo@meta.data
save(endo.metadata,file = 'endo.final.metadata.2023-01-29.Robj')
gc()
rm(endo)
gc()


##### Mes Subclustering and Cleaning #####
mes <- subset(lung.combined,idents = 'Mesenchyme')
DefaultAssay(mes) <- 'integrated' # Let's try using the old integration-space
mes <- ScaleData(mes)
mes<- RunPCA(mes, npcs = 100)
pdf(file='mes.round2.PCs.pdf',width=10,height=8)
ElbowPlot(mes,ndims = 100)
PCHeatmap(mes,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
mes <- RunUMAP(mes, reduction = "pca", dims = 1:12)
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:12)
DefaultAssay(mes) <- 'integrated'
mes <- FindClusters(mes, resolution = 0.4)
DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)

# Visualization
pdf(file='mes.round2.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(mes, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes) <- 'RNA'
p4 <- FeaturePlot(mes,features = c('Lgr5','Lgr6','Col13a1','Col14a1',
                                   'Acta2','Gucy1a1','Wnt11','Msln'),label = T,ncol = 3)
p5 <- VlnPlot(mes,features = c('Lgr5','Lgr6','Col13a1','Col14a1',
                               'Acta2','Gucy1a1','Wnt11','Msln'),pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()


# Identify Cell Types
DefaultAssay(mes) <- 'RNA'
mes.mark <- FindAllMarkers(mes,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
mes.mark$ratio <- mes.mark$pct.1/mes.mark$pct.2
#save(mes.mark,file='mes.final.mark.Robj')
p1 <- DimPlot(mes,group.by = 'CellType',label=T)
p2 <- DimPlot(mes,label = T)
p1+p2
mes <- RenameIdents(mes,
                    '0'='Col13a1_Fib',
                    '1'='Myofibroblasts',
                    '2'='Col14a1_Fib',
                    '3'='Mesothelium',
                    '4'='Col13a1_Fib',
                    '5'='SMCs',
                    '6'='Lgr5_Fib',
                    '7'='Pericytes',
                    '8'='Wfdc21_Mesothelium')
mes$CellType_Final <- Idents(mes)
DimPlot(mes,label = T)
pdf(file='mes.final.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(mes, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes) <- 'RNA'
p4 <- FeaturePlot(mes,features = c('Col13a1','Col14a1','Gucy1a1','Acta2','Wnt11','Msln'),label = T,ncol = 3)
p5 <- VlnPlot(mes,features = c('Col13a1','Col14a1','Gucy1a1','Acta2','Wnt11','Msln'),pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
gc()
save(mes,file = 'mes.final.2023-01-29.Robj')
mes.metadata <- mes@meta.data
save(mes.metadata,file = 'mes.final.metadata.2023-01-29.Robj')
gc()
rm(mes)

##### Imm Subclustering and Cleaning #####
imm <- subset(lung.combined,idents = 'Immune')
DefaultAssay(imm) <- 'integrated' # Let's try using the old integration-space
imm <- ScaleData(imm)
imm<- RunPCA(imm, npcs = 100)
pdf(file='imm.round2.PCs.pdf',width=10,height=8)
ElbowPlot(imm,ndims = 100)
PCHeatmap(imm,cells=200,balanced=T,dims=1:9)
PCHeatmap(imm,cells=200,balanced=T,dims=10:18)
PCHeatmap(imm,cells=200,balanced=T,dims=19:27)
PCHeatmap(imm,cells=200,balanced=T,dims=28:36)
PCHeatmap(imm,cells=200,balanced=T,dims=37:45)
PCHeatmap(imm,cells=200,balanced=T,dims=46:54)
PCHeatmap(imm,cells=200,balanced=T,dims=55:63)
PCHeatmap(imm,cells=200,balanced=T,dims=64:72)
PCHeatmap(imm,cells=200,balanced=T,dims=73:81)
PCHeatmap(imm,cells=200,balanced=T,dims=82:90)
PCHeatmap(imm,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
imm <- RunUMAP(imm, reduction = "pca", dims = 1:25)
imm <- FindNeighbors(imm, reduction = "pca", dims = 1:25)
DefaultAssay(imm) <- 'integrated'
imm <- FindClusters(imm, resolution = 0.6)
DimPlot(imm,label = T)

# Visualization
pdf(file='imm.round2.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(imm, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE)
p1 +p2 + p3
dev.off()

# Identify Cell Types
DefaultAssay(imm) <- 'RNA'
imm.mark <- FindAllMarkers(imm,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
imm.mark$ratio <- imm.mark$pct.1/imm.mark$pct.2
#save(imm.mark,file='imm.final.mark.Robj')
p1 <- DimPlot(imm,group.by = 'CellType',label=T)
p2 <- DimPlot(imm,label = T)
p1+p2
imm <- RenameIdents(imm,
                    '0'='B',
                    '1'='T',
                    '2'='Mac_Alv',
                    '3'='Monocytes',
                    '4'='NK',
                    '5'='Mac_Alv',
                    '6'='Monocytes',
                    '7'='Neutrophils',
                    '8'='Mac_Inter',
                    '9'='DC',
                    '10'='ILCs',
                    '11'='Killer_T_Prolif',
                    '12'='Mac_Prolif',
                    '13'='Siglech+',
                    '14'='Plasma',
                    '15'='Multiplet',
                    '16'='Megakaryocytes')
imm$CellType_Final <- Idents(imm)

pdf(file='imm.Round2.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(imm, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 + p3
dev.off()
DimPlot(imm)
# Save and Clean up
save(imm,file = 'imm.final.2023-01-29.Robj')
imm.metadata <- imm@meta.data
save(imm.metadata,file = 'imm.final.metadata.2023-01-29.Robj')
rm(imm)
gc()

#### Map metadat back to integrated object ####
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/imm.final.metadata.2023-01-29.Robj")
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/endo.final.metadata.2023-01-29.Robj")
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/epi.final.metadata.2023-01-29.Robj")
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/mes.final.metadata.2023-01-29.Robj")

epi.meta <- data.frame(CellType_Final = epi.final.metadata$CellType_Final)
rownames(epi.meta) <- rownames(epi.final.metadata)

endo.meta <- data.frame(CellType_Final = endo.metadata$CellType_Final)
rownames(endo.meta) <- rownames(endo.metadata)

mes.meta <- data.frame(CellType_Final = mes.metadata$CellType_Final)
rownames(mes.meta) <- rownames(mes.metadata)

imm.meta <- data.frame(CellType_Final = imm.metadata$CellType_Final)
rownames(imm.meta) <- rownames(imm.metadata)

celltype.final <- rbind(epi.meta,endo.meta,imm.meta,mes.meta)

lung.combined <- AddMetaData(lung.combined,metadata = celltype.final)

# Final clean
Idents(lung.combined) <- lung.combined$CellType_Final
table(Idents(lung.combined),lung.combined$CellClass_Final)
lung.combined <- subset(lung.combined,idents = 'Multiplet',invert=T) # Remove 68 cells that are Multiplet
table(Idents(lung.combined),lung.combined$CellClass_Final)

# Final touches
table(Idents(lung.combined))
lung.combined <- RenameIdents(lung.combined,
                              'Killer_T_Prolif' = 'Lymphoid_Prolif')
table(Idents(lung.combined))
lung.combined$CellType_Final <- Idents(lung.combined)

lung.combined$CellType_Final <- factor(lung.combined$CellType_Final,
                                       levels = unique(lung.combined$CellType_Final)[stringr::str_order(unique(lung.combined$CellType_Final))])
# Visualization
png(file='lung.combined.final.UMAP.annotated.png',width=20,height=18,units = 'in',res=400)
p1 <- DimPlot(lung.combined, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(lung.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(lung.combined, reduction = "umap", group.by = "CellClass_Final")
p4 <- DimPlot(lung.combined, reduction = "umap",  group.by = "CellType_Final",label = TRUE, repel = TRUE)
p1 + p2 + p3 + p4
dev.off()

# Final save
save(lung.combined,file='lung.combined.clean.classed.annotated.final.2023-01-29.Robj')

##### Tuft BASC Niche ####
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/lung.combined.clean.classed.annotated.final.2023-01-29.Robj")
require(NICHES)
require(Seurat)
require(dplyr)
require(SeuratWrappers)
# lung.al <- RunALRA(lung.combined)
# save(lung.al,file = 'lung.combined.imputed.2023-03-29.Robj')
# table(Idents(lung.al))
scc <- RunNICHES(lung.combined,
                 assay = 'RNA',
                 species = 'rat',
                 LR.database = 'fantom5',
                 cell_types = 'CellType_Final',
                 meta.data.to.map = names(lung.combined@meta.data))
ctc <- scc$CellToCell
#VlnPlot(ctc,'nFeature_CellToCell',raster=F,log=T,pt.size = 0)+NoLegend()
Idents(ctc) <- ctc$ReceivingType
table(Idents(ctc))
sub <- subset(ctc,idents = c('BASC','Tuft'))
VlnPlot(sub,'nFeature_CellToCell',raster=F,log=T,pt.size = 0)+NoLegend()
sub <- subset(sub,nFeature_CellToCell > 25)
Idents(sub) <- sub$SendingType
table(Idents(sub))
sub <- ScaleData(sub)
sub <- FindVariableFeatures(sub)
sub <- RunPCA(sub)
ElbowPlot(sub,ndims = 50)
PCHeatmap(sub,dims=1:9,balanced=T,cells=300)
PCHeatmap(sub,dims=10:18,balanced=T,cells=300)
PCHeatmap(sub,dims=19:27,balanced=T,cells=300)
sub <- RunUMAP(sub,dims = c(1:2,4,6:7,9:12)) # Do not use degenerate PCs (driven by a single receptor, ie.)
DimPlot(sub)
FeaturePlot(sub,'nCount_CellToCell')
DimPlot(sub,group.by = 'Dataset.Joint')+NoLegend()
sub <- FindNeighbors(sub,dims = c(1:2,4,6:7,9:12),force.recalc = T)
sub <- FindClusters(sub,resolution = 0.4)
DimPlot(sub,label=T)
mark <- Seurat::FindAllMarkers(sub,slot = 'scale.data')
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio*mark$avg_diff
View(mark)
p1 <- DimPlot(sub,label=T,group.by = 'SendingType')+NoLegend()
p2 <- DimPlot(sub,label=T,group.by = 'ReceivingType')+NoLegend()
p3 <- DimPlot(sub,label=T)+NoLegend()
cowplot::plot_grid(p1,p2,p3)
FeaturePlot(sub,'Rspo1—Lgr5')
FeaturePlot(sub,'Rspo1—Lrp6',pt.size = 1,order = T)
FeaturePlot(sub,'Wnt5a—Fzd6',pt.size = 1,order = T)

VlnPlot(sub,'Rspo1—Lrp6',group.by = 'SendingType',sort=T)

index.remove <- grep('itg',mark$gene,ignore.case = T)
MOI <- mark[-index.remove,] %>% group_by(cluster) %>% top_n(10,power)
View(MOI)
DoHeatmap(sub,features = MOI$gene)
# try imputing the niches data??
sub.al <- SeuratWrappers::RunALRA(sub)
sub.al <- ScaleData(sub.al)
DoHeatmap(sub.al,features = MOI$gene)
