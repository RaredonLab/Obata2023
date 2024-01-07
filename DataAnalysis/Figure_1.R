# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Homeostatic_Single_Cell")

# Packages
require(Seurat)
require(ggplot2)
require(RColorBrewer)
require(cowplot)
library(scales)
library(viridis)
set.seed(123)

# Load homeostatic atlas data
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/lung.combined.clean.classed.annotated.final.2023-01-29.Robj")

# Inspect
table(lung.combined$CellType)
table(lung.combined$CellClass)
Idents(lung.combined) <- 'CellClass'
table(Idents(lung.combined))

# Break into classes
epi <- subset(lung.combined, idents='Epithelium')
end <- subset(lung.combined, idents='Endothelium')
mes <- subset(lung.combined, idents='Mesenchyme')
imm <- subset(lung.combined, idents='Immune')

# Order cell types
celltypes.ordered <- c(sort(names(table(epi$CellType))),
                      sort(names(table(end$CellType))),
                      sort(names(table(mes$CellType))),
                      sort(names(table(imm$CellType))))

# Define colors
col.pal <- list()
col.pal$Class <- c(brewer.pal(4,'Set1'))
names(col.pal$Class) <- c('Epithelium','Endothelium','Mesenchyme','Immune')
col.pal$Dataset <- c(brewer.pal(3,'Dark2'))
names(col.pal$Dataset) <- c('10x_v2','10x_v3','DropSeq')
col.pal$Type <- c('firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue',
                             'orangered','violetred','tomato','grey20','sandybrown',
                             'saddlebrown','royalblue','plum4','lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink',
                             'red2','paleturquoise1','palevioletred','orchid4','purple4','plum1','olivedrab2',
                             'slateblue','mediumvioletred','sienna','orange','seagreen',
                             'lightseagreen','mediumpurple4')
names(col.pal$Type) <- celltypes.ordered
# Order the metadata in the object
lung.combined$CellType <- factor(lung.combined$CellType,levels = celltypes.ordered)

# Plot UMAPS (total tissue)
p1 <- DimPlot(lung.combined,group.by = 'CellType',cols = col.pal$Type,pt.size = 0.75,shuffle = T)+ggtitle(NULL)+
        theme(plot.title=element_text(hjust=0))+
        theme(plot.title = element_text(size = 30, face = "bold"))+
        guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
        NoAxes()
p2 <- DimPlot(lung.combined,group.by = 'CellClass',cols = col.pal$Class,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p3 <- DimPlot(lung.combined,group.by = 'Dataset',cols = col.pal$Dataset,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p4 <- DimPlot(lung.combined,group.by = 'Sample',shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
p5 <- plot_grid(p2,p3,p4,ncol = 1)
png(filename = 'homeostatic.adult.rat.lung.UMAPs.combined.png',width = 16,height = 10,res = 400,units = 'in')
plot_grid(p1,p5,rel_widths = c(2.6,1))
dev.off()

# Epi Heatmap
features <- list()
DefaultAssay(epi) <- 'RNA'
epi <- ScaleData(epi,features = rownames(epi))
Idents(epi) <- epi$CellType
# mark.epi <- FindAllMarkers(epi,only.pos = T,min.pct = 0.1)
# mark.epi$ratio <- mark.epi$pct.1/mark.epi$pct.2
# mark.epi$power <- mark.epi$ratio*mark.epi$avg_log2FC
# View(mark.epi)
# mark.basc <- FindMarkers(epi,ident.1 = 'BASC',ident.2 = 'Tuft')
# mark.basc$ratio <- mark.basc$pct.1/mark.basc$pct.2
# View(mark.basc)
epi$CellType <- factor(epi$CellType,levels = c('Tuft','BASC','ATII','ATII-ATI','ATI','Ciliated','Secretory'))
Idents(epi) <- epi$CellType
sum(is.na(epi$CellType))
features$epi <- c('Sox9','Dclk1','Pou2f3','Rgs13','Gng13','Avil','Krt8','Krt18','Il25','Plac8','Nrgn','Trpm5','Lgr5',
                  'Abca3','Sftpc','Sftpa1','Sftpd','Sftpb','Fabp5',
                  'Foxq1','Edn3','Pdzk1ip1','Adora2b','Plekhd1',
                  'Nkx2-1','Ager','Pdpn','Hopx','Vegfa','Sema3e',
                  'Ccdc153','Dnah6','Dnai1','Ccdc187','Spef2',
                  'Scgb1a1','Scgb3a2','Muc20','Lgi1','Bpifa5')
png(filename = 'homeostatic.adult.rat.epi.heatmap.png',width = 15,height = 8,res = 400,units = 'in')
DoHeatmap(epi,features$epi,group.colors = col.pal$Type[unique(epi$CellType)])
dev.off()

# Epi
DefaultAssay(epi) <- "integrated"
epi <- ScaleData(epi)
epi <- FindVariableFeatures(epi,selection.method = 'disp',assay = 'integrated',nfeatures = 2000)
epi <- RunPCA(epi,features = rownames(epi))
ElbowPlot(epi,ndims = 50)
DimHeatmap(epi,dims=1:9,balanced=T,cells=200)
DimHeatmap(epi,dims=10:18,balanced=T,cells=200)
DimHeatmap(epi,dims=19:27,balanced=T,cells=200)
epi <- RunUMAP(epi,dims = 1:20)
png(filename = 'homeostatic.adult.rat.epi.UMAP.png',width = 6,height = 5,res = 400,units = 'in')
DimPlot(epi,cols = col.pal$Type,pt.size = 1,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
dev.off()

# endo
Idents(end) <- end$CellType
DefaultAssay(end) <- "integrated"
end <- ScaleData(end)
end <- FindVariableFeatures(end,selection.method = 'disp',assay = 'integrated',nfeatures = 2000)
end <- RunPCA(end,features = rownames(end))
ElbowPlot(end,ndims = 50)
DimHeatmap(end,dims=1:9,balanced=T,cells=200)
DimHeatmap(end,dims=10:18,balanced=T,cells=200)
DimHeatmap(end,dims=19:27,balanced=T,cells=200)
end <- RunUMAP(end,dims = 3:9)
png(filename = 'homeostatic.adult.rat.end.UMAP.png',width = 6,height = 5,res = 400,units = 'in')
DimPlot(end,cols = col.pal$Type,pt.size = 1,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()+
  xlim(-7.5,10)+
  ylim(-10,10)
dev.off()

# mes
Idents(mes) <- mes$CellType
DefaultAssay(mes) <- "integrated"
mes <- ScaleData(mes)
mes <- FindVariableFeatures(mes,selection.method = 'disp',assay = 'integrated',nfeatures = 2000)
mes <- RunPCA(mes,features = rownames(mes))
ElbowPlot(mes,ndims = 50)
DimHeatmap(mes,dims=1:9,balanced=T,cells=200)
DimHeatmap(mes,dims=10:18,balanced=T,cells=200)
DimHeatmap(mes,dims=19:27,balanced=T,cells=200)
mes <- RunUMAP(mes,dims = 1:20)
png(filename = 'homeostatic.adult.rat.mes.UMAP.png',width = 6,height = 5,res = 400,units = 'in')
DimPlot(mes,cols = col.pal$Type,pt.size = 1,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
dev.off()

# imm
Idents(imm) <- imm$CellType
DefaultAssay(imm) <- "integrated"
imm <- ScaleData(imm)
imm <- FindVariableFeatures(imm,selection.method = 'disp',assay = 'integrated',nfeatures = 2000)
imm <- RunPCA(imm,features = rownames(imm))
ElbowPlot(imm,ndims = 50)
DimHeatmap(imm,dims=1:9,balanced=T,cells=200)
DimHeatmap(imm,dims=10:18,balanced=T,cells=200)
DimHeatmap(imm,dims=19:27,balanced=T,cells=200)
imm <- RunUMAP(imm,dims = 1:20)
png(filename = 'homeostatic.adult.rat.imm.UMAP.png',width = 6,height = 5,res = 400,units = 'in')
DimPlot(imm,cols = col.pal$Type,pt.size = 1,shuffle = T)+ggtitle(NULL)+
  theme(plot.title=element_text(hjust=0))+
  theme(plot.title = element_text(size = 30, face = "bold"))+
  guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
  NoAxes()
dev.off()

# Epi feature plots
DefaultAssay(epi) <- 'RNA'
pt.size = 0.1
p1 <- FeaturePlot(epi,'Sox9',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p2 <- FeaturePlot(epi,'Lgr5',max.cutoff = 2,pt.size = pt.size,order=T)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p3 <- FeaturePlot(epi,'Pou2f3',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p4 <- FeaturePlot(epi,'Dclk1',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p5 <- FeaturePlot(epi,'Il25',max.cutoff = 2,pt.size = pt.size,order=T)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p6 <- FeaturePlot(epi,'Avil',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p7 <- FeaturePlot(epi,'Rgs13',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p8 <- FeaturePlot(epi,'Gng13',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p9 <- FeaturePlot(epi,'Trpm5',max.cutoff = 2,pt.size = pt.size,order=T)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p10 <- FeaturePlot(epi,'Nrgn',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p11 <- FeaturePlot(epi,'Plac8',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
p12 <- FeaturePlot(epi,'Lrmp',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))+NoLegend()
png(filename = 'homeostatic.adult.rat.epi.FEATURES.png',width = 11,height = 8,res = 400,units = 'in')
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,nrow=3)
dev.off()
png(filename = 'homeostatic.adult.rat.epi.FEATURES.legend.png',width = 5,height = 5,res = 400,units = 'in')
FeaturePlot(epi,'Sox9',max.cutoff = 2,pt.size = pt.size,order=F)+
  DarkTheme()+NoAxes()+paletteer::scale_color_paletteer_c("viridis::viridis",limits=c(0,2))
dev.off()
