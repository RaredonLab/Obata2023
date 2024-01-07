# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Homeostatic_Single_Cell")

# Packages
require(Seurat)
require(ggplot2)
require(RColorBrewer)
require(cowplot)
require(dplyr)
library(scales)
library(viridis)
require(Connectome)
set.seed(123)

# Load homeostatic atlas data
load("~/T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2023-01-20/lung.combined.clean.classed.annotated.final.2023-01-29.Robj")

# Inspect
table(lung.combined$CellType)
table(lung.combined$CellClass)
Idents(lung.combined) <- 'CellClass'
table(Idents(lung.combined))

# Order cell types
sort(unique(lung.combined@meta.data[WhichCells(lung.combined,expression=CellClass=='Epithelium'),]$CellType))
            
celltypes.ordered <- c(sort(unique(lung.combined@meta.data[WhichCells(lung.combined,expression=CellClass=='Epithelium'),]$CellType)),
                       sort(unique(lung.combined@meta.data[WhichCells(lung.combined,expression=CellClass=='Endothelium'),]$CellType)),
                       sort(unique(lung.combined@meta.data[WhichCells(lung.combined,expression=CellClass=='Mesenchyme'),]$CellType)),
                       sort(unique(lung.combined@meta.data[WhichCells(lung.combined,expression=CellClass=='Immune'),]$CellType)))

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

# Set idents as cell types
Idents(lung.combined) <- lung.combined$CellType

# Run Connectome
connectome.genes <- union(Connectome::ncomms8866_rat$Ligand.ApprovedSymbol,Connectome::ncomms8866_rat$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(lung.combined)]
DefaultAssay(lung.combined) <- 'RNA'
lung.combined <- ScaleData(lung.combined,features = genes)
lung.con <- CreateConnectome(lung.combined,species = 'rat',min.cells.per.ident = 40,p.values = F,calculate.DOR = F)

# Run NICHES

table(lung.combined$CellClass)
table(lung.combined$CellType)
table(lung.combined$Sample)
table(lung.combined$Dataset)

#1. Split by sample
split <- SplitObject(lung.combined,split.by = 'Sample')

#2. Impute by sample
options(warn = 1)
for(i in 1:length(split)){
  gc()
  print(i)
  split[[i]] <- SeuratWrappers::RunALRA(split[[i]])
}

# Save imputed data
gc()
save(split,file='homeostatic.lung.split.imputed.Robj')

#3. Run NICHES, with all three outputs selected
set.seed(123)
scc.list <- list()
for(i in 1:length(split)){
  print(i)
  scc.list[[i]] <- NICHES::RunNICHES(split[[i]],
                             LR.database="fantom5",
                             species="rat",
                             assay="alra",
                             cell_types = "CellType",
                             meta.data.to.map = names(split[[i]]@meta.data),
                             SystemToCell = T,
                             CellToCell = T,
                             CellToSystem=T)
}
names(scc.list) <- names(split)     

#4. Merge outputs together for each analysis
temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToCell # Isolate CellToCell Signaling
}
cell.to.cell <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.cell

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$SystemToCell # Isolate SystemToCell Signaling
}
system.to.cell <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
system.to.cell

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToSystem # Isolate CellToSystem Signaling
}
cell.to.system <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.system

# Save point
save(cell.to.cell,file='cell.to.cell.Robj')
save(system.to.cell,file='system.to.cell.Robj')
save(cell.to.system,file='cell.to.system.Robj')

# Microenvironment
system.to.cell
system.to.cell <- ScaleData(system.to.cell)
table(Idents(system.to.cell))
mark.micro <- FindAllMarkers(system.to.cell,only.pos = T,min.pct = 0.04,logfc.threshold = 0.05)
mark.micro$ratio <- mark.micro$pct.1/mark.micro$pct.2
mark.micro$power <- mark.micro$ratio*mark.micro$avg_log2FC
View(mark.micro)

# Organize for plotting
custom.celltype.order <- c( "Tuft" ,'BASC', "ATII","ATII-ATI","ATI",  "Secretory",   "Ciliated",         
                            "aCap","Arterial","gCap","Lymphatic","Venous",
                   "Col13a1_Fib","Col14a1_Fib","Lgr5_Fib","Mesothelium","Myofibroblasts","Pericytes","SMCs","Wfdc21_Mesothelium",
                  "B","DC","ILCs","Killer_T_Prolif",
                   "Mac_Alv","Mac_Inter","Mac_Prolif","Megakaryocytes","Monocytes","Neutrophils",
                  "NK","Plasma","Siglech+","T")
col.pal$Type <- col.pal$Type[custom.celltype.order]
system.to.cell$CellClass <- factor(system.to.cell$CellClass,
                                      levels=c('Epithelium','Endothelium','Mesenchyme','Immune'))
system.to.cell$CellType <- factor(system.to.cell$CellType,
                                   levels=custom.celltype.order)

# Make small
small.sub <- subset(system.to.cell,cells=WhichCells(system.to.cell,downsample = 1000))

# Select features
MOI <- mark.micro %>% group_by(cluster) %>% top_n(7,ratio) %>% arrange(factor(cluster, levels = custom.celltype.order))
MOI$gene
MOI.custom <- c( "Rspo1—Lgr5","Il17b—Il17rb", "Il25—Il17rb",#"Sct—Sctr", # Tuft
                 'Hgf—Met', #BASC & ATII
                   "Wnt3a—Fzd2", # ATI
                  "Bmp2—Bmpr1b", #Proximal Epithelium
                  "Vegfa—Kdr", #aCap 
                 "Ptn—Ptprz1", #Arterial
                 #'Apln—Aplnr', # gCap
                 'Ntf3—Ntrk2',# gCap see https://www.science.org/doi/10.1126/scisignal.adf6696 induces proliferation!
                 #'Vegfc—Flt4', #Lymphatic
                 #'Efnb2—Ephb4', #Venous
                 'Cdh1—Cdh2', #Venous
                 'Pdgfa—Pdgfra',#Col13
                 'Rbp4—Stra6', # Col14 & Lgr5
                 #'Rbp4—Stra6',#Lgr5
                 'Rtn4—Rtn4r',#Meso
                 #'Ntf4—Ntrk3',#Myo
                 'Rspo1—Lgr6', #SMC & Myo
                 #'Shh—Smo',#Wfdc
                 'C3—Cd19',#B
                 'Cxcl2—Xcr1', #DC
                 'Ccl3—Ccr3', #ILC
                 'Cxcl2—Cxcr1', #Mac Alv
                 #'Csf2—Csf2rb', #Mac Int
                 'Il16—Cd4', #Mono
                 'Cxcl2—Cxcr2', # Neutro
                 'Il15—Il2rb',#NK
                 'B2m—Cd3d'#T
                 )

# Plot - Microenvironment heatmap
png('full.lung.microenvironment.heatmap.png',width=12,height = 7,units = 'in',res=600)
CustomHeatmap(object = small.sub,
              data.type = 'SystemToCell',
              primary = 'CellClass' ,
              secondary = 'CellType' ,
              tertiary = 'Dataset' ,
              quarternary = 'Sample' ,
              primary.cols = col.pal$Class,
              secondary.cols = col.pal$Type, # Need to be a named list of colors in the right order
              tertiary.cols = col.pal$Dataset,
              quarternary.cols = NULL,
              features = MOI.custom,
              labels = c('Cell Class','Microenvironment','Dataset','Sample'))
dev.off()



#5. Clean niches data
table(Idents(cell.to.cell))
table(Idents(cell.to.cell),cell.to.cell$Sample.Joint)
VlnPlot(cell.to.cell,
        features = 'nFeature_CellToCell',
        group.by = 'Sample.Joint',
        pt.size=0,log = T)
temp <- subset(cell.to.cell,nFeature_CellToCell>50)
temp
VlnPlot(temp,
        features = 'nFeature_CellToCell',
        group.by = 'Sample.Joint',
        pt.size=0,log = F)

# #6. Integrate by dataset (needed, in this instance)
# # split the dataset into a list of two seurat objects (stim and CTRL)
# niche.list <- SplitObject(temp, split.by = "Dataset.Sending")
# 
# # identify variable features for each dataset independently
# niche.list <- lapply(X = niche.list, FUN = function(x) {
#   x <- FindVariableFeatures(x)
# })
# 
# # select features that are repeatedly variable across datasets for integration run PCA on each
# # dataset using these features
# features <- SelectIntegrationFeatures(object.list = niche.list)
# niche.list <- lapply(X = niche.list, FUN = function(x) {
#   x <- ScaleData(x, features = features, verbose = FALSE)
#   x <- RunPCA(x, features = features, verbose = FALSE)
# })
# niche.anchors <- FindIntegrationAnchors(object.list = niche.list, anchor.features = features, reduction = "rpca")
# 
# # this command creates an 'integrated' data assay
# niche.combined <- IntegrateData(anchorset = niche.anchors)
# DefaultAssay(niche.combined) <- "integrated"
# gc()
# 
# #7. Visualization
# # Perform initial visualization
# niche.combined <- ScaleData(niche.combined)
# #niche.combined <- FindVariableFeatures(niche.combined)
# niche.combined <- RunPCA(niche.combined,npcs = 100)
# ElbowPlot(niche.combined,ndim=100)
# PCHeatmap(niche.combined,balanced=T,cells=200,dims=1:9)
# PCHeatmap(niche.combined,balanced=T,cells=200,dims=10:18)
# PCHeatmap(niche.combined,balanced=T,cells=200,dims=19:27)
# PCHeatmap(niche.combined,balanced=T,cells=200,dims=28:36)
# PCHeatmap(niche.combined,balanced=T,cells=200,dims=37:45)
# niche.combined <- RunUMAP(niche.combined,dims = 1:40)
# 
# # Plotting: Total connectomics
# p1 <- DimPlot(niche.combined,group.by = 'Dataset.Sending',shuffle = T,raster = F,cols = col.pal$Dataset)+NoAxes()+ggtitle('Dataset')
# p2 <- DimPlot(niche.combined,group.by = 'Sample.Sending',shuffle = T,raster = F)+NoLegend()+NoAxes()+ggtitle('Sample')
# p3 <- DimPlot(niche.combined,group.by = 'CellClass.Sending',shuffle = T,label = F,raster = F,cols = col.pal$Class)+NoAxes()+ggtitle('Sending Cell Class')
# p4 <- DimPlot(niche.combined,group.by = 'CellClass.Receiving',shuffle = T,label = F,raster = F,cols = col.pal$Class)+NoAxes()+ggtitle('Receiving Cell Class')
# png(filename = 'homeostatic.signaling.total.umap.png',width = 14,height = 12,res = 600,units = 'in')
# cowplot::plot_grid(p1,p2,p3,p4,align = T)
# dev.off()

######### Tuft BASC Niche ONLY
tuft.niche <- subset(temp,subset=CellType.Receiving==c('Tuft','BASC'))
# Integrate by dataset (trying on small subset here)
# split the dataset into a list of two seurat objects (stim and CTRL)
niche.list <- SplitObject(tuft.niche, split.by = "Dataset.Sending")

# identify variable features for each dataset independently
niche.list <- lapply(X = niche.list, FUN = function(x) {
  x <- FindVariableFeatures(x)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = niche.list)
niche.list <- lapply(X = niche.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
niche.anchors <- FindIntegrationAnchors(object.list = niche.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
niche.combined <- IntegrateData(anchorset = niche.anchors)
DefaultAssay(niche.combined) <- "integrated"
gc()

tuft.niche <- ScaleData(niche.combined)
tuft.niche <- FindVariableFeatures(tuft.niche,selection.method = 'disp')
tuft.niche <- RunPCA(tuft.niche)
ElbowPlot(tuft.niche,ndims = 50)
PCHeatmap(tuft.niche,balanced=T,cells=300,dims=1:9)
PCHeatmap(tuft.niche,balanced=T,cells=300,dims=10:18)
PCHeatmap(tuft.niche,balanced=T,cells=200,dims=19:27)
PCHeatmap(tuft.niche,balanced=T,cells=200,dims=28:36)
PCHeatmap(tuft.niche,balanced=T,cells=200,dims=37:45)
tuft.niche <- RunUMAP(tuft.niche,dims = c(1:5,7:15))
tuft.niche <- FindNeighbors(tuft.niche,dims = c(1:5,7:15))
tuft.niche <- FindClusters(tuft.niche,resolution = 0.2)
DimPlot(tuft.niche)
# Remove degenerate cluster 8
tuft.niche <- subset(tuft.niche,idents = '8',invert=T)

# FindMarkers
DefaultAssay(tuft.niche) <- 'CellToCell'
tuft.niche <- ScaleData(tuft.niche)
tuft.niche.mark <- FindAllMarkers(tuft.niche,only.pos=T,slot = 'scale.data')
tuft.niche.mark$ratio <- tuft.niche.mark$pct.1/tuft.niche.mark$pct.2
tuft.niche.mark$power <- tuft.niche.mark$ratio*tuft.niche.mark$avg_diff
View(tuft.niche.mark)

# Colors for feature plot with alpha
paletteFunc <- colorRampPalette(c('lightgrey','blue'));
palette     <- paletteFunc(100);
palette <- alpha(palette,0.50)

# Plot NICHES UMAPS for publication
tuft.niche$Signaling.Archetype <- tuft.niche$seurat_clusters
p1 <- DimPlot(tuft.niche,group.by = 'Signaling.Archetype',pt.size = 1.5,label=T,label.size = 7,shuffle = T)+ggtitle('Signaling Archetype')+NoLegend()+NoAxes()
p2 <- DimPlot(tuft.niche,group.by = 'CellType.Sending',pt.size = 1.5,cols = col.pal$Type,shuffle = T)+NoLegend()+NoAxes()+ggtitle('Sending Cell Type')
p3 <- DimPlot(tuft.niche,group.by = 'CellType.Receiving',pt.size = 1.5,cols = col.pal$Type,shuffle = T)+NoLegend()+NoAxes()+ggtitle('Receiving Cell Type')
p4 <- FeaturePlot(tuft.niche,'Rspo1—Lgr5',order = F,pt.size = 1.5,max.cutoff = 2,cols = palette)+NoLegend()+NoAxes()
p5 <- FeaturePlot(tuft.niche,'Il17b—Il17rb',order = F,pt.size = 1.5,max.cutoff = 2,cols = palette)+NoLegend()+NoAxes()
p6 <- FeaturePlot(tuft.niche,'Wnt3a—Fzd6',order = F,pt.size = 1.5,max.cutoff = 2,cols = palette)+NoLegend()+NoAxes()
png(filename = 'tuft.basc.signaling.umap.png',width = 10,height = 7,res = 600,units = 'in')
cowplot::plot_grid(p1,p2,p3,p4,p5,p6,align = T,nrow = 2)
dev.off()

# Select genes for tuft basc niches heatmap
MOI.to.plot <- tuft.niche.mark %>% group_by(cluster) %>% top_n(30,avg_diff)
MOI.to.plot <- MOI.to.plot[order(MOI.to.plot$cluster,-MOI.to.plot$power),]

# Setup tuft.niche levels and colors for plotting
tuft.niche$CellType.Sending <- factor(tuft.niche$CellType.Sending,levels = custom.celltype.order)
col.pal$Type <- col.pal$Type[custom.celltype.order]


# Plot - Tuft BASc Signaling Archetypes Heatmap
png('tuft.basc.niche.archetypes.heatmap.png',width=12,height = 7,units = 'in',res=600)
CustomHeatmap(object = tuft.niche,
              data.type = 'CellToCell',
              primary = 'Signaling.Archetype' ,
              secondary = 'CellType.Receiving' ,
              tertiary = 'CellType.Sending' ,
              quarternary = 'CellClass.Sending' ,
              primary.cols = NULL,
              secondary.cols = col.pal$Type, # Need to be a named list of colors in the right order
              tertiary.cols = col.pal$Type,
              quarternary.cols = col.pal$Class,
              features = unique(MOI.to.plot$gene),
              labels = c('Signaling Archetype','Receiving Cell Type','Sending Cell Type','Sending Cell Class'),
              selected.row.anotations = c('Il17b—Il17rb','Il25—Il17rb',
                                          'Ereg—Erbb3','Areg—Erbb3',
                                          'Sema4d—Plxnb2','Sema6a—Plxna2',#'Sema6a—Plxna4',
                                          'Wnt3a—Fzd6','Wnt7b—Lrp5','Wnt5a—Fzd5','Wnt5a—Fzd6','Wnt3a—Ryk',
                                          'Rspo2—Lgr5','Rspo3—Sdc4','Rspo1—Lrp6','Rspo1—Lrp6','Rspo1—Lgr5',
                                          'Dll1—Notch2','Dll4—Notch2','Dll3—Notch2',
                                          'Sct—Sctr','Dhh—Ptch1',
                                          'Tgfb1—Tgfbr1','Tgfb1—Tgfbr2','Ccl2—Ackr4','Ccl5—Sdc1',#'Ccl5—Sdc4',
                                          'Mdk—Sdc4','Ccl11—Ackr4','Fbln1—Itgb1',
                                          'Il18—Cd48','Tnf—Ltbr','Igf1—Insr','Tnfsf13—Tnfrsf13b'),
              selected.label.size = 10)
dev.off()

