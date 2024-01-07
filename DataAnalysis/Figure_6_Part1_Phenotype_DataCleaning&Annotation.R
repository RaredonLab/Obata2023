# Figure 6 Part 1 Pneumonectomy Data Cleaning and Annotation

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell")
# Set Seed
set.seed(123)
# Load Packages
require(Seurat)
require(RColorBrewer)
require(ggplot2)
require(cowplot)

# Load Data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneum.data.cleaned.2023-03-01.Robj")
# Inspect
names(pneum@meta.data)
table(pneum$Sample) # Only has sample, no class or type information yet

#### Determine Cell Class ####
pneum <- ScaleData(pneum)
pneum <- FindVariableFeatures(pneum)
pneum <- RunPCA(pneum,npcs = 100)

# Look at PCs
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

# Try a PC selection and see how it looks in 2D space
pneum <- RunUMAP(pneum,dims = 1:53)
DimPlot(pneum,group.by = 'Sample')

# Look at how the four classes segregate
DefaultAssay(pneum) <- 'RNA'
FeaturePlot(pneum,features = c('Epcam','Ptprc','Col1a1','Cdh5'))
FeaturePlot(pneum,'Mmrn1')
FeaturePlot(pneum,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# We like this! The four classes segregate and we can see lymphatics and other rare populations
# Find neighbors and cluster
pneum <- FindNeighbors(pneum,dims = 1:53)
pneum <- FindClusters(pneum,resolution = 0.2)
DimPlot(pneum)
p1 <- DimPlot(pneum,label=T)
p2 <- FeaturePlot(pneum,features = c('Epcam','Ptprc','Col1a1','Cdh5'),label = T)
p3 <- VlnPlot(pneum,features = c('Epcam','Ptprc','Col1a1','Cdh5'),pt.size=0,ncol = 4)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,nrow=2)

# Label based on cell class:
cell.epi <- WhichCells(pneum, idents = c(4,5,12,15))
cell.mes <- WhichCells(pneum, idents = c(11,13,16,26))
cell.endo <- WhichCells(pneum, idents = c(0,19,23))
cell.immu <- WhichCells(pneum, idents = c(1:3,6:10,17:18,20:22,24:25))
cell.cycle <- WhichCells(pneum,idents = 14) # There is a distinct cell cycle cluster here, we will have to break up and redistribute 
pneum <- SetIdent(pneum, cells = cell.epi, value = 'Epithelium')
pneum <- SetIdent(pneum, cells = cell.mes, value = 'Mesenchyme')
pneum <- SetIdent(pneum, cells = cell.endo, value = 'Endothelium')
pneum <- SetIdent(pneum, cells = cell.immu, value = 'Immune')
pneum <- SetIdent(pneum, cells = cell.cycle, value = 'Cell_cycle')
DimPlot(pneum,label=T)

# Break up the cell cycle cluster:
cycle <- subset(pneum,idents = 'Cell_cycle')

# Determine Cell Class within the cell cycle cluster
cycle <- ScaleData(cycle)
cycle <- FindVariableFeatures(cycle)
cycle <- RunPCA(cycle,npcs = 100)

# Look at PCs
pdf(file='cycle.PCs.pdf',width=10,height=8)
ElbowPlot(cycle,ndims = 100)
PCHeatmap(cycle,cells=200,balanced=T,dims=1:9)
PCHeatmap(cycle,cells=200,balanced=T,dims=10:18)
PCHeatmap(cycle,cells=200,balanced=T,dims=19:27)
PCHeatmap(cycle,cells=200,balanced=T,dims=28:36)
PCHeatmap(cycle,cells=200,balanced=T,dims=37:45)
PCHeatmap(cycle,cells=200,balanced=T,dims=46:54)
PCHeatmap(cycle,cells=200,balanced=T,dims=55:63)
PCHeatmap(cycle,cells=200,balanced=T,dims=64:72)
PCHeatmap(cycle,cells=200,balanced=T,dims=73:81)
PCHeatmap(cycle,cells=200,balanced=T,dims=82:90)
PCHeatmap(cycle,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
cycle <- RunUMAP(cycle,dims = 1:20)
DimPlot(cycle,group.by = 'Sample')

# Look at how the four classes segregate
DefaultAssay(cycle) <- 'RNA'
FeaturePlot(cycle,features = c('Epcam','Ptprc','Col1a1','Cdh5'))
FeaturePlot(cycle,'Mmrn1')
FeaturePlot(cycle,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# We like this! The four classes segregate and we can see lymphatics and other rare populations
# Find neighbors and cluster
cycle <- FindNeighbors(cycle,dims = 1:20,force.recalc = T)
cycle <- FindClusters(cycle,resolution = 0.2)
DimPlot(cycle)
p1 <- DimPlot(cycle,label=T)
p2 <- FeaturePlot(cycle,features = c('Epcam','Ptprc','Col1a1','Cdh5'),label = T)
p3 <- VlnPlot(cycle,features = c('Epcam','Ptprc','Col1a1','Cdh5'),pt.size=0,ncol = 4)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,nrow=2)
# Label the cell cyle cluster based on cell class:
Idents(cycle) <- cycle$seurat_clusters
cell.epi <- WhichCells(cycle, idents = c(3))
cell.endo <- WhichCells(cycle, idents = c(1))
cell.immu <- WhichCells(cycle, idents = c(0,2,4,5,7))
cell.multi <- WhichCells(cycle,idents = 6) # Multiplets (and a tiny bit of mesenchyme we aren't going to worry about)
cycle <- SetIdent(cycle, cells = cell.epi, value = 'Epithelium')
cycle <- SetIdent(cycle, cells = cell.endo, value = 'Endothelium')
cycle <- SetIdent(cycle, cells = cell.immu, value = 'Immune')
cycle <- SetIdent(cycle, cells = cell.multi, value = 'Multiplet')
DimPlot(cycle,label=T)

#### Epithelium Cluster, Clean and Annotate ####
epi1 <- subset(pneum,idents='Epithelium')
epi2 <- subset(cycle,idents = 'Epithelium') 
epi <- merge(epi1,epi2)

# Cluster
epi <- ScaleData(epi)
epi <- FindVariableFeatures(epi)
epi <- RunPCA(epi,npcs = 100)

# Look at PCs
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

# Try a PC selection and see how it looks in 2D space
epi <- RunUMAP(epi,dims = 1:25)
DimPlot(epi,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(epi,features = c('Sox9','Abca3','Aqp5','Scgb1a1'))
FeaturePlot(epi,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# Looks good. Some sample separation, appears to be male specific?
# Find neighbors and cluster
epi <- FindNeighbors(epi,dims = 1:25)
epi <- FindClusters(epi,resolution = 0.6)
DimPlot(epi)
p1 <- DimPlot(epi,label=T)
p2 <- FeaturePlot(epi,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
#p2 <- FeaturePlot(epi,features = c('Sox9','Sftpc','Scgb1a1','Aqp5','Dclk1','Ccdc153'),label = T)
p3 <- VlnPlot(epi,features = c('Sox9','Sftpc','Scgb1a1','Aqp5','Dclk1','Ccdc153'),pt.size=0,ncol = 6)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,nrow=2)

# Plot to look at
png(filename = 'temp.png',width=16,height = 12,units = 'in',res=300)
cowplot::plot_grid(p4,p3,nrow=2)
dev.off()

# Remove multiplets
cell.multi <- WhichCells(epi, idents = c(13,15,16,18,19))
epi <- SetIdent(epi, cells = cell.multi, value = 'Multiplet')
epi.sub <- subset(epi,idents = 'Multiplet',invert=T)

# Cluster sub
epi.sub <- ScaleData(epi.sub)
epi.sub <- FindVariableFeatures(epi.sub)
epi.sub <- RunPCA(epi.sub,npcs = 100)

# Look at PCs
pdf(file='epi.sub.PCs.pdf',width=10,height=8)
ElbowPlot(epi.sub,ndims = 100)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi.sub,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
epi.sub <- RunUMAP(epi.sub,dims = 1:15)
DimPlot(epi.sub,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(epi.sub,features = c('Sox9','Abca3','Aqp5','Scgb1a1'))
FeaturePlot(epi.sub,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# Looks good. Some sample separation, appears to be male specific?
# Find neighbors and cluster
epi.sub <- FindNeighbors(epi.sub,dims = 1:15)
epi.sub <- FindClusters(epi.sub,resolution = 0.4)
DimPlot(epi.sub)
p1 <- DimPlot(epi.sub,label=T)
p5 <- VlnPlot(epi.sub,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 4)
p2 <- FeaturePlot(epi.sub,features = c('Sox9','Sftpc','Scgb1a1','Aqp5','Dclk1','Ccdc153'),label = T)
p3 <- VlnPlot(epi.sub,features = c('Sox9','Sftpc','Scgb1a1','Aqp5','Dclk1','Ccdc153'),pt.size=0,ncol = 6)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,p5,nrow=3)

# Plot to look at
png(filename = 'temp.png',width=16,height = 12,units = 'in',res=300)
cowplot::plot_grid(p4,p3,p5,nrow=3)
dev.off()

# Find Markers
mark.epi.sub <- FindAllMarkers(epi.sub,min.pct = 0.2,logfc.threshold = 0.25,only.pos = T)
mark.epi.sub$ratio <- mark.epi.sub$pct.1/mark.epi.sub$pct.2
mark.epi.sub$power <- mark.epi.sub$ratio*mark.epi.sub$avg_log2FC

#  Annotate
epi.sub <- RenameIdents(epi.sub,
                    '0'='ATII',
                    '1'='ATI',
                    '2'='ATI',
                    '3'='ATI',
                    '4'='Tuft',
                    '5'='Ciliated',
                    '6'='ATII-ATI',
                    '7'='ATII',
                    '8'='ATII',
                    '9'='Secretory',
                    '10'='BASC',
                    '11'='ATII',
                    '12'='Tuft',
                    '13'='Cell_cycle',
                    '14'='Multiplet')
DimPlot(epi.sub)
save(epi.sub,file = 'pneum.epi.sub.annotated.2023-06-14.Robj')
gc()

#### Endothelium Cluster, Clean and Annotate ####
end1 <- subset(pneum,idents='Endothelium')
end2 <- subset(cycle,idents = 'Endothelium') 
end <- merge(end1,end2)

# Cluster
end <- ScaleData(end)
end <- FindVariableFeatures(end)
end <- RunPCA(end,npcs = 100)

# Look at PCs
pdf(file='end.PCs.pdf',width=10,height=8)
ElbowPlot(end,ndims = 100)
PCHeatmap(end,cells=200,balanced=T,dims=1:9)
PCHeatmap(end,cells=200,balanced=T,dims=10:18)
PCHeatmap(end,cells=200,balanced=T,dims=19:27)
PCHeatmap(end,cells=200,balanced=T,dims=28:36)
PCHeatmap(end,cells=200,balanced=T,dims=37:45)
PCHeatmap(end,cells=200,balanced=T,dims=46:54)
PCHeatmap(end,cells=200,balanced=T,dims=55:63)
PCHeatmap(end,cells=200,balanced=T,dims=64:72)
PCHeatmap(end,cells=200,balanced=T,dims=73:81)
PCHeatmap(end,cells=200,balanced=T,dims=82:90)
PCHeatmap(end,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
end <- RunUMAP(end,dims = 1:25)
DimPlot(end,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(end,features = c('Prx','Lyve1','Mmrn1','Gja5'))
FeaturePlot(end,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))
# Looks good. Some sample separation. Definitely multiplets and low-info in here still
# Find neighbors and cluster
end <- FindNeighbors(end,dims = 1:25)
end <- FindClusters(end,resolution = 0.6)
DimPlot(end)
VlnPlot(end,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol=4,pt.size = 0)
p1 <- DimPlot(end,label=T)
p2 <- FeaturePlot(end,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
p3 <- VlnPlot(end,features = c('Epcam','Col1a1','Cdh5','Ptprc'),pt.size=0,ncol = 4)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,nrow=2)

# Plot to look at
png(filename = 'temp.png',width=16,height = 12,units = 'in',res=300)
cowplot::plot_grid(p4,p3,nrow=2)
dev.off()

# Remove multiplets and low info cells
cell.lowinfo <- WhichCells(end, idents = c(1,6))
cell.multi <- WhichCells(end, idents = c(13,14,16))
end <- SetIdent(end, cells = cell.lowinfo, value = 'LowInfo')
end <- SetIdent(end, cells = cell.multi, value = 'Multiplet')
end.sub <- subset(end,idents = c('Multiplet','LowInfo'),invert=T)

# Cluster sub
end.sub <- ScaleData(end.sub)
end.sub <- FindVariableFeatures(end.sub)
end.sub <- RunPCA(end.sub,npcs = 100)

# Look at PCs
pdf(file='end.sub.PCs.pdf',width=10,height=8)
ElbowPlot(end.sub,ndims = 100)
PCHeatmap(end.sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(end.sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(end.sub,cells=200,balanced=T,dims=19:27)
PCHeatmap(end.sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(end.sub,cells=200,balanced=T,dims=37:45)
PCHeatmap(end.sub,cells=200,balanced=T,dims=46:54)
PCHeatmap(end.sub,cells=200,balanced=T,dims=55:63)
PCHeatmap(end.sub,cells=200,balanced=T,dims=64:72)
PCHeatmap(end.sub,cells=200,balanced=T,dims=73:81)
PCHeatmap(end.sub,cells=200,balanced=T,dims=82:90)
PCHeatmap(end.sub,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
end.sub <- RunUMAP(end.sub,dims = 1:15)
DimPlot(end.sub,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(end.sub,features = c('Prx','Lyve1','Mmrn1','Gja5'))
FeaturePlot(end.sub,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# Looks good. Some sample separation, appears again to be male specific?
# Find neighbors and cluster
end.sub <- FindNeighbors(end.sub,dims = 1:15)
end.sub <- FindClusters(end.sub,resolution = 0.4)
DimPlot(end.sub)
p1 <- DimPlot(end.sub,label=T)
p5 <- VlnPlot(end.sub,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 4)
p2 <- FeaturePlot(end.sub,features = c('Prx','Lyve1','Mmrn1','Gja5'),label = T)
p3 <- VlnPlot(end.sub,features = c('Prx','Lyve1','Mmrn1','Gja5'),pt.size=0,ncol = 6)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,p5,nrow=3)

# Plot to look at
png(filename = 'temp.png',width=16,height = 12,units = 'in',res=300)
cowplot::plot_grid(p4,p3,p5,nrow=3)
dev.off()

# Find Markers
mark.end.sub <- FindAllMarkers(end.sub,min.pct = 0.2,logfc.threshold = 0.25,only.pos = T)
mark.end.sub$ratio <- mark.end.sub$pct.1/mark.end.sub$pct.2
mark.end.sub$power <- mark.end.sub$ratio*mark.end.sub$avg_log2FC

#  Annotate
end.sub <- RenameIdents(end.sub,
                        '0'='gCap',
                        '1'='Arterial',
                        '2'='gCap',
                        '3'='Venous',
                        '4'='aCap',
                        '5'='gCap',
                        '6'='Cell_cycle',
                        '7'='Lymphatic',
                        '8'='aCap',
                        '9'='Venous')
DimPlot(end.sub)
save(end.sub,file = 'pneum.end.sub.annotated.2023-06-14.Robj')
gc()
####  Mesenchyme Cluster, Clean and Annotate #### 
mes1 <- subset(pneum,idents='Mesenchyme')
#mes2 <- subset(cycle,idents = 'Mesenchyme') # No mesenchyme from the earlier cell cycle object
mes <- mes1

# Cluster
mes <- ScaleData(mes)
mes <- FindVariableFeatures(mes)
mes <- RunPCA(mes,npcs = 100)

# Look at PCs
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

# Try a PC selection and see how it looks in 2D space
mes <- RunUMAP(mes,dims = 1:20)
DimPlot(mes,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(mes,features = c('Msln','Col13a1','Col14a1','Acta2','Gucy1a1','Wnt11'))
FeaturePlot(mes,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))
# Looks good. Some sample separation. Definitely multiplets and low-info in here still
# Find neighbors and cluster
mes <- FindNeighbors(mes,dims = 1:20)
mes <- FindClusters(mes,resolution = 0.8)
DimPlot(mes)
VlnPlot(mes,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol=4,pt.size = 0)
FeaturePlot(mes,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),label=T)
p1 <- DimPlot(mes,label=T)
p2 <- FeaturePlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
p3 <- VlnPlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'),pt.size=0,ncol = 4)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,nrow=2)

# Plot to look at
png(filename = 'temp.png',width=16,height = 12,units = 'in',res=300)
cowplot::plot_grid(p4,p3,nrow=2)
dev.off()

# Remove multiplets and low info cells
cell.lowinfo <- WhichCells(mes, idents = c(0,16,22)) # Keep 17 and 18 for now...might be EMT?
cell.multi <- WhichCells(mes, idents = c(19,20,23))
mes <- SetIdent(mes, cells = cell.lowinfo, value = 'LowInfo')
mes <- SetIdent(mes, cells = cell.multi, value = 'Multiplet')
mes.sub <- subset(mes,idents = c('Multiplet','LowInfo'),invert=T)

# Cluster sub
mes.sub <- ScaleData(mes.sub)
mes.sub <- FindVariableFeatures(mes.sub)
mes.sub <- RunPCA(mes.sub,npcs = 100)

# Look at PCs
pdf(file='mes.sub.PCs.pdf',width=10,height=8)
ElbowPlot(mes.sub,ndims = 100)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=19:27)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=37:45)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=46:54)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=55:63)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=64:72)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=73:81)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=82:90)
PCHeatmap(mes.sub,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
mes.sub <- RunUMAP(mes.sub,dims = 1:20)
DimPlot(mes.sub,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(mes.sub,features = c('Msln','Col13a1','Col14a1','Acta2','Gucy1a1','Wnt11'))
FeaturePlot(mes.sub,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# Looks good. Some sample separation, appears again to be male specific?
# Find neighbors and cluster
mes.sub <- FindNeighbors(mes.sub,dims = 1:20)
mes.sub <- FindClusters(mes.sub,resolution = 0.4)
DimPlot(mes.sub)
p1 <- DimPlot(mes.sub,label=T)
p5 <- VlnPlot(mes.sub,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 4)
p2 <- FeaturePlot(mes.sub,features = c('Msln','Col13a1','Col14a1','Acta2','Gucy1a1','Wnt11'),label = T,ncol = 3)
p3 <- VlnPlot(mes.sub,features = c('Msln','Col13a1','Col14a1','Acta2','Gucy1a1','Wnt11'),pt.size=0,ncol = 6)
p6 <- VlnPlot(mes.sub,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol = 4)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,p5,p6,nrow=4)

# Plot to look at
png(filename = 'temp.png',width=16,height = 16,units = 'in',res=300)
cowplot::plot_grid(p4,p3,p5,p6,nrow=4)
dev.off()

# Find Markers
mark.mes.sub <- FindAllMarkers(mes.sub,min.pct = 0.2,logfc.threshold = 0.25,only.pos = T)
mark.mes.sub$ratio <- mark.mes.sub$pct.1/mark.mes.sub$pct.2
mark.mes.sub$power <- mark.mes.sub$ratio*mark.mes.sub$avg_log2FC

#  Annotate
mes.sub <- RenameIdents(mes.sub,
                        '0'='Myofibroblasts',
                        '1'='LowInfo',
                        '2'='Col14a1_Fib',
                        '3'='Col13a1_Fib',
                        '4'='Mesothelium',
                        '5'='Lgr5_Fib',
                        '6'='Mesothelium',
                        '7'='Myofibroblasts',
                        '8'='Pericytes',
                        '9'='Col14a1_Fib',
                        '10'='Mesothelium',
                        '11'='Frzb+/Smoc1+',
                        '12'='Multiplet',
                        '13'='Multiplet',
                        '14'='Multiplet',
                        '15'='Multiplet')
DimPlot(mes.sub)
save(mes.sub,file = 'pneum.mes.sub.annotated.2023-06-14.Robj')
gc()


####  Immune Cluster, Clean and Annotate #### 
imm1 <- subset(pneum,idents='Immune')
imm2 <- subset(cycle,idents = 'Immune') 
imm <- merge(imm1,imm2)

# Cluster
imm <- ScaleData(imm)
imm <- FindVariableFeatures(imm)
imm <- RunPCA(imm,npcs = 100)

# Look at PCs
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

# Try a PC selection and see how it looks in 2D space
imm <- RunUMAP(imm,dims = 1:50)
DimPlot(imm,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(imm,features = c('C1qc','Naaa','Cd3g','Nkg7','Cd79a','Jchain'))
FeaturePlot(imm,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))
FeaturePlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(imm,features = c('Sox9','Ager','Sftpc','Dclk1'))
# Looks good. Some sample separation. Definitely multiplets and low-info in here still
# Find neighbors and cluster
imm <- FindNeighbors(imm,dims = 1:50)
imm <- FindClusters(imm,resolution = 0.6)
DimPlot(imm)
VlnPlot(imm,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol=4,pt.size = 0)
FeaturePlot(imm,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),label=T)
p1 <- DimPlot(imm,label=T)
p2 <- FeaturePlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
p3 <- VlnPlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'),pt.size=0,ncol = 4)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
p5 <- VlnPlot(imm,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol=4,pt.size = 0)
cowplot::plot_grid(p4,p3,p5,nrow=3)

# Plot to look at
png(filename = 'temp.png',width=24,height = 18,units = 'in',res=300)
cowplot::plot_grid(p4,p3,p5,nrow=3)
dev.off()

# Remove multiplets and low info cells
cell.lowinfo <- WhichCells(imm, idents = c(23)) 
cell.multi <- WhichCells(imm, idents = c(25,29,30))
imm <- SetIdent(imm, cells = cell.lowinfo, value = 'LowInfo')
imm <- SetIdent(imm, cells = cell.multi, value = 'Multiplet')
imm.sub <- subset(imm,idents = c('Multiplet','LowInfo'),invert=T)

# Cluster sub
imm.sub <- ScaleData(imm.sub)
imm.sub <- FindVariableFeatures(imm.sub)
imm.sub <- RunPCA(imm.sub,npcs = 100)

# Look at PCs
pdf(file='imm.sub.PCs.pdf',width=10,height=8)
ElbowPlot(imm.sub,ndims = 100)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=1:9)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=10:18)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=19:27)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=28:36)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=37:45)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=46:54)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=55:63)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=64:72)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=73:81)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=82:90)
PCHeatmap(imm.sub,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
imm.sub <- RunUMAP(imm.sub,dims = 1:50)
DimPlot(imm.sub,group.by = 'Sample')

# Look at how a few key markers segregate
FeaturePlot(imm.sub,features = c('C1qc','Naaa','Cd3g','Nkg7','Cd79a','Jchain'))
FeaturePlot(imm.sub,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'))

# Looks good. Some sample separation, again appears again to be male specific
# Find neighbors and cluster
imm.sub <- FindNeighbors(imm.sub,dims = 1:50)
imm.sub <- FindClusters(imm.sub,resolution = 0.4)
DimPlot(imm.sub)
p1 <- DimPlot(imm.sub,label=T)
p5 <- VlnPlot(imm.sub,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 4,pt.size = 0)
p2 <- FeaturePlot(imm.sub,features = c('C1qc','Naaa','Cd3g','Nkg7','Cd79a','Jchain'),label = T,ncol = 3)
p3 <- VlnPlot(imm.sub,features = c('C1qc','Naaa','Cd3g','Nkg7','Cd79a','Jchain'),pt.size=0,ncol = 6)
p6 <- VlnPlot(imm.sub,features = c('nFeature_RNA','nCount_RNA','percent.mt','perc.spliced'),ncol = 4,pt.size = 0)
p4 <- cowplot::plot_grid(p1,p2,ncol = 2)
cowplot::plot_grid(p4,p3,p5,p6,nrow=4)

# Plot to look at
png(filename = 'temp.png',width=24,height = 24,units = 'in',res=300)
cowplot::plot_grid(p4,p3,p5,p6,nrow=4)
dev.off()

# Find Markers
mark.imm.sub <- FindAllMarkers(imm.sub,min.pct = 0.2,logfc.threshold = 0.25,only.pos = T)
mark.imm.sub$ratio <- mark.imm.sub$pct.1/mark.imm.sub$pct.2
mark.imm.sub$power <- mark.imm.sub$ratio*mark.imm.sub$avg_log2FC

#  Annotate
imm.sub <- RenameIdents(imm.sub,
                        '0'='B',
                        '1'='Mac_Alv',
                        '2'='T',
                        '3'='Monocytes',
                        '4'='Neutrophils',
                        '5'='NK',
                        '6'='Monocytes',
                        '7'='T',
                        '8'='NK',
                        '9'='B',
                        '10'='DC',
                        '11'='Mac_Inter',
                        '12'='Mac_Alv',
                        '13'='Monocytes',
                        '14'='T',
                        '15'='Mac_Alv',
                        '16'='Cell_cycle',
                        '17'='ILC',
                        '18'='Monocytes',
                        '19'='Siglech+',
                        '20'='Mac_Inter',
                        '21'='Rag1+',
                        '22'='Mac_Alv',
                        '23'='Mast_Il4+',
                        '24'='Plasma',
                        '25'='Mast')
DimPlot(imm.sub)
save(imm.sub,file = 'pneum.imm.sub.annotated.2023-06-14.Robj')
gc()

# Compile classes together into one object
epi.sub$CellType <- Idents(epi.sub)
end.sub$CellType <- Idents(end.sub)
mes.sub$CellType <- Idents(mes.sub)
imm.sub$CellType <- Idents(imm.sub)
epi.sub$CellClass <- "Epithelium"
end.sub$CellClass <- "Endothelium"
mes.sub$CellClass <- "Mesenchyme"
imm.sub$CellClass <- "Immune"
gc()
pneum.clean <- merge(epi.sub,list(end.sub,mes.sub,imm.sub))

# Add Timepoint metadata
Idents(pneum.clean) <- pneum.clean$Sample
table(Idents(pneum.clean))
cell.pn0 <- WhichCells(pneum.clean, idents = c('P0_f','P0_m')) 
cell.pn3 <- WhichCells(pneum.clean, idents = c('P3_f','P3_m')) 
cell.pn7 <- WhichCells(pneum.clean, idents = c('P7_f','P7_m')) 
cell.pn14 <- WhichCells(pneum.clean, idents = c('P14_f','P14_m')) 
pneum.clean <- SetIdent(pneum.clean, cells = cell.pn0, value = 'Day 0')
pneum.clean <- SetIdent(pneum.clean, cells = cell.pn3, value = 'Day 3')
pneum.clean <- SetIdent(pneum.clean, cells = cell.pn7, value = 'Day 7')
pneum.clean <- SetIdent(pneum.clean, cells = cell.pn14, value = 'Day 14')
pneum.clean$Timepoint <- Idents(pneum.clean)
pneum.clean$Timepoint <- factor(pneum.clean$Timepoint,levels=c('Day 0','Day 3','Day 7','Day 14'))
table(Idents(pneum.clean))

# Remove Multiplets and LowInfo
Idents(pneum.clean) <- pneum.clean$CellType
table(Idents(pneum.clean))
pneum.clean <- subset(pneum.clean,idents = c('LowInfo','Multiplet'),invert=T)
table(Idents(pneum.clean))

# Embed
pneum.clean <- ScaleData(pneum.clean)
pneum.clean <- FindVariableFeatures(pneum.clean)
pneum.clean <- RunPCA(pneum.clean,npcs = 100)

# Look at PCs
pdf(file='pneum.clean.PCs.pdf',width=10,height=8)
ElbowPlot(pneum.clean,ndims = 100)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=1:9)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=10:18)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=19:27)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=28:36)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=37:45)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=46:54)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=55:63)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=64:72)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=73:81)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=82:90)
PCHeatmap(pneum.clean,cells=200,balanced=T,dims=91:99)
dev.off()

# Try a PC selection and see how it looks in 2D space
pneum.clean <- RunUMAP(pneum.clean,dims = 1:50)
DimPlot(pneum.clean,group.by = 'Sample')
DimPlot(pneum.clean,group.by = 'CellClass')
DimPlot(pneum.clean,group.by = 'CellType',label=T)
p1 <- DimPlot(pneum.clean,group.by = 'CellType',label=T)
p2 <- DimPlot(pneum.clean,group.by = 'CellClass',label=T)
p3 <- DimPlot(pneum.clean,group.by = 'Timepoint',label=T)
p4 <- DimPlot(pneum.clean,group.by = 'Sample',label=T)
pdf(file='pneum.clean.UMAPS.pdf',width=24,height=12)
cowplot::plot_grid(p1,p2,p3,p4)
dev.off()

# Save
save(pneum.clean,file = 'pneum.clean.annotated.2023-06-14.Robj')
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneum.clean.annotated.2023-06-14.Robj")

#### Figure 6 Plotting for publication ####
# Order cell types
epi <- subset(pneum.clean@meta.data,CellClass == 'Epithelium')
ec <- subset(pneum.clean@meta.data,CellClass == 'Endothelium')
mes <- subset(pneum.clean@meta.data,CellClass =='Mesenchyme')
imm <- subset(pneum.clean@meta.data,CellClass == 'Immune')

celltypes.ordered <- c(sort(as.character(unique(epi$CellType))),
                       sort(as.character(unique(ec$CellType))),
                       sort(as.character(unique(imm$CellType))),
                       sort(as.character(unique(mes$CellType))))
celltypes.ordered <- c(celltypes.ordered[celltypes.ordered!='Cell_cycle'],'Cell_cycle')
pneum.clean$CellType <- factor(pneum.clean$CellType, levels = celltypes.ordered)
table(pneum.clean$CellType)
Idents(pneum.clean) <- pneum.clean$CellType
table(Idents(pneum.clean))

# Define colors
col.pal <- list()
col.pal$Class <- c(brewer.pal(4,'Set1'))
names(col.pal$Class) <- c('Epithelium','Endothelium','Mesenchyme','Immune')
col.pal$Timepoint <- c(brewer.pal(4,'Dark2'))
names(col.pal$Timepoint) <- c('Day 0','Day 3','Day 7','Day 14')
col.pal$CellType <- c('firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue',
                             'orangered','violetred','tomato','grey20','sandybrown',
                             'saddlebrown','royalblue','plum4','lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink',
                             'red2','paleturquoise1','palevioletred','orchid4','purple4','plum1','olivedrab2',
                             'slateblue','mediumvioletred','sienna','orange','seagreen',
                             'lightseagreen','mediumpurple4')
names(col.pal$CellType) <- celltypes.ordered

# Plot UMAPS (total pneumonectomy timecourse)
p1 <- DimPlot(pneum.clean,group.by = 'CellType',cols = col.pal$CellType,pt.size = 0.75,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p2 <- DimPlot(pneum.clean,group.by = 'CellClass',cols = col.pal$Class,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p3 <- DimPlot(pneum.clean,group.by = 'Timepoint',cols = col.pal$Timepoint,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p4 <- DimPlot(pneum.clean,group.by = 'Sample',shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p5 <- plot_grid(p2,p3,p4,ncol = 1,align = 'hv')
png(filename = 'pneumonectomy.rat.lung.UMAPs.combined.png',width = 16,height = 10,res = 400,units = 'in')
plot_grid(p1,p5,rel_widths = c(2.3,1))
dev.off()

# Quantifying Cell Proportion by Timepoint

# Format metadata for plotting
pneum.clean$Timepoint <- factor(pneum.clean$Timepoint,
                                levels = c('Day 0','Day 3','Day 7','Day 14'))
pneum.clean$Sample <- factor(pneum.clean$Sample,
                                levels = c('P0_f','P0_m','P3_f','P3_m',
                                           'P7_f','P7_m','P14_f','P14_m'))
# Cell type ratios
epi <- subset(pneum.clean,subset=CellClass=='Epithelium')
end <- subset(pneum.clean,subset=CellClass=='Endothelium')
imm <- subset(pneum.clean,subset=CellClass=='Immune')
mes <- subset(pneum.clean,subset=CellClass=='Mesenchyme')

# Global class proportion
table <- as.data.frame(prop.table(table(pneum.clean$Timepoint,pneum.clean$CellClass),margin = 1))
ggplot(data=table,
       aes(x=Var1,y=Freq,fill=Var2))+
  geom_col(position = 'fill')+
  scale_fill_manual(values=col.pal$Class)+
  theme_minimal()+
  ggtitle('Class Proportion') ## Immune cells go up and Epitheliual cells go down

# Epithelial proportion
table <- as.data.frame(prop.table(table(epi$Timepoint,epi$CellType),margin = 1))
ggplot(data=table,
       aes(x=Var1,y=Freq,fill=Var2))+geom_col(position = 'fill')+ggtitle('Epithelial Proportion') ## ATII-ATI cells increase and then decrease, ATII cells are up by D14
table <- as.data.frame(prop.table(prop.table(table(epi$Timepoint,epi$CellType),margin = 1),margin = 2))
ggplot(data=table,
       aes(x=Var2,y=Freq,fill=Var1))+geom_col(position = 'fill')
# Mesenchymal proportion
table <- as.data.frame(prop.table(table(mes$Timepoint,mes$CellType),margin = 1))
ggplot(data=table,
       aes(x=Var1,y=Freq,fill=Var2))+geom_col(position = 'fill')+ggtitle('Mesenchymal Proportion')  ## Expansion of mesothelium (and pericytes!) and decrease in myofibroblasts!
table <- as.data.frame(prop.table(prop.table(table(mes$Timepoint,mes$CellType),margin = 1),margin = 2))
ggplot(data=table,
       aes(x=Var2,y=Freq,fill=Var1))+geom_col(position = 'fill')
# Immune proportion
table <- as.data.frame(prop.table(table(imm$Timepoint,imm$CellType),margin = 1))
ggplot(data=table,
       aes(x=Var1,y=Freq,fill=Var2))+geom_col(position = 'fill')+ggtitle('Immune Proportion')  ## Increase in alveolar macrophages
table <- as.data.frame(prop.table(prop.table(table(imm$Timepoint,imm$CellType),margin = 1),margin = 2))
ggplot(data=table,
       aes(x=Var2,y=Freq,fill=Var1))+geom_col(position = 'fill')
# Endothelial proportion
table <- as.data.frame(prop.table(table(end$Timepoint,end$CellType),margin = 1))
ggplot(data=table,
       aes(x=Var1,y=Freq,fill=Var2))+geom_col(position = 'fill')+ggtitle('Endothelial Proportion')  ## Increase in cell cycle at day 7, drop in gCap (probably because more are entering cell cycle), constant aCap proportion
table <- as.data.frame(prop.table(prop.table(table(end$Timepoint,end$CellType),margin = 1),margin = 2))
ggplot(data=table,
       aes(x=Var2,y=Freq,fill=Var1))+geom_col(position = 'fill')

# Marker approach development
# Subset mesothelium
meso <- subset(pneum.clean,idents='Mesothelium')
Idents(meso) <- meso$Timepoint
mark.meso <- FindAllMarkers(meso,only.pos = F) ## Rspo1 down and Wnt4 up
# Subset ATI
ati <- subset(pneum.clean,idents='ATI')
Idents(ati) <- ati$Timepoint
mark.ati <- FindAllMarkers(ati,only.pos = F) ## very little Ati change
# Subset myofib and peri
myo <- subset(pneum.clean,idents=c('Myofibroblasts','Pericytes'))
Idents(myo) <- myo$Timepoint
mark.myo <- FindAllMarkers(myo,only.pos = F) # minimal change
# Subset mac alv
mac <- subset(pneum.clean,idents='Mac_Alv')
Idents(mac) <- mac$Timepoint
mark.mac <- FindAllMarkers(mac,only.pos = F) ## 
