# Set WD
setwd("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Tuft_Project_Data_Cleaning_2022-08-17")

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
library(monocle3)
library(ComplexHeatmap)
library(cowplot)
library(velocyto.R)

#### SetUp ####
# Functions (from Allie) # Some are not fully generalized #
LoadSeuratData <- function(file.names,sample.names){
  data <- list()
  for (i in 1:length(file.names)){
    message(sample.names[i])
    load(file = paste("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Pneumonectomy_Project/",file.names[i],sep='')) ## Script specific workaround
    data[[i]] <- CreateSeuratObject(counts = output[['Gene']])
    data[[i]][['GeneFull']] <- CreateAssayObject(output[['GeneFull']])
    data[[i]][['spliced']] <- CreateAssayObject(output[['spliced']])
    data[[i]][['unspliced']] <- CreateAssayObject(output[['unspliced']])
    data[[i]][["percent.mt"]] <- PercentageFeatureSet(data[[i]], pattern = "^Mt-")
    data[[i]]$perc.spliced <- 100*(colSums(data[[i]]@assays$spliced)/(colSums(data[[i]]@assays$spliced)+colSums(data[[i]]@assays$unspliced)))
    data[[i]]$Sample <- sample.names[i]
    gc()
  }
  return(data)
}
LoadMetaData <- function(seurat.data){
  data <- data.frame()
  for (i in 1:length(seurat.data)){
    data <- rbind(data,seurat.data[[i]]@meta.data)
    gc()
  }
  return(data)
}
QCViolinPlotter <- function(sample.name,min.counts,min.features,max.mt,tag){
  x <- which(sample.name == names(seurat.data))
  object <- seurat.data[[x]]
  filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Counts = ",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Features = ",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=50) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Max.Mt = ",max.mt)) +
    labs(caption = paste("Barcodes retained with these filters from starting total: ",ncol(filtered),"/",ncol(object))) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'),plot.caption = element_text(color='red'))
  png(paste(sample.name,'_filterselect_1',tag,'.png',sep=""),width = 1500,height=500)
  print(plot_grid(p,q,r,ncol=3))
  dev.off()
}
TriplePercentMtPlotter <- function(sample.name,min.counts,min.features,max.mt,tag){
  object <- subset(pneum.metadata,Sample %in% sample.name)
  lowpass <- subset(object,percent.mt < max.mt)
  highmt <- subset(object,percent.mt > max.mt)
  filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
  p <- ggplot() +
    geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
    geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
    geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(title = sample.name) +
    labs(subtitle = paste("Min.Counts = ",min.counts," / Min.Features = ",min.features," / Max.Mt = ",max.mt)) +
    labs(caption = paste("Barcodes retained with these filters from starting total: ",nrow(filtered),"/",nrow(object))) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          plot.subtitle = element_text(color='red'),
          plot.caption = element_text(color='red'),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  png(paste(sample.name,'_filterselect_2',tag,'.png',sep=""),width = 600,height=500)
  print(p)
  dev.off()
}
TriplePercentSplPlotter <- function(sample.name,min.counts,min.features,max.mt,tag){
  object <- subset(pneum.metadata,Sample %in% sample.name)
  lowpass <- subset(object,percent.mt < max.mt)
  highmt <- subset(object,percent.mt > max.mt)
  filtered <- subset(object,nFeature_RNA > min.features & nCount_RNA > min.counts & percent.mt < max.mt)
  p <- ggplot() +
    geom_point(data=highmt, aes(x=nCount_RNA, y=nFeature_RNA), color = 'gray') +
    geom_point(data=lowpass, aes(x=nCount_RNA, y=nFeature_RNA, color=perc.spliced)) +
    scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
    geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(title = sample.name) +
    labs(subtitle = paste("Min.Counts = ",min.counts," / Min.Features = ",min.features," / Max.Mt = ",max.mt)) +
    labs(caption = paste("Barcodes retained with these filters from starting total: ",nrow(filtered),"/",nrow(object))) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          plot.subtitle = element_text(color='red'),
          plot.caption = element_text(color='red'),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  png(paste(sample.name,'_filterselect_3',tag,'.png',sep=""),width = 600,height=500)
  print(p)
  dev.off()
}

# More Functions (from Sam)  # Some may not be fully generalized #
theProcess <- function(sample.name,min.counts,min.features,max.mt,tag){
  # Plot QC metrics
  QCViolinPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)
  TriplePercentMtPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)
  TriplePercentSplPlotter(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)
  
  # Pull object, apply subset, & dimensional reduction
  x <- which(sample.name == names(seurat.data))
  object <- seurat.data[[x]]
  object <- subset(object, nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object) # don't regress on anything for now, see if it matters
  object <- RunPCA(object, npcs = 50, verbose = F)
  png(paste(sample.name,'_elbow.png',tag,sep=""), width = 800, height = 500)
  print(ElbowPlot(object,ndims = 50)+labs(title = sample.name))
  dev.off()
  pdf(paste(sample.name,'_pcheatmap_1',tag,'.pdf',sep=""), width = 12, height = 16)
  print(DimHeatmap(object, dims = 1:15, cells = 200, balanced = T))
  dev.off()
  pdf(paste(sample.name,'_pcheatmap_2',tag,'.pdf',sep=""), width = 12, height = 16)
  print(DimHeatmap(object, dims = 16:30, cells = 200, balanced = T))
  dev.off()
  
  # Cluster, plot, & check QC
  pcs <- 30
  object <- FindNeighbors(object, dims = 1:pcs)
  object <- FindClusters(object, resolution = 0.5)
  object <- RunUMAP(object, reduction = "pca", dims = 1:pcs)
  png(paste(sample.name,'_umap_1',tag,'.png',sep=""), width = 500, height = 500)
  print(UMAPPlot(object,label = T) + labs(title = sample.name,subtitle = paste("PC's = ",pcs)) + NoLegend() + NoAxes())
  dev.off()
  png(paste(sample.name,'_fp_lin',tag,'.png',sep=""), width = 800, height = 800)
  print(FeaturePlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'), label = T,repel=T))
  dev.off()
  png(paste(sample.name,'_fp_qc',tag,'.png',sep=""), width = 800, height = 1200)
  print(FeaturePlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,label = T,repel=T))
  dev.off()
  png(paste(sample.name,'_vln_lin',tag,'.png',sep=""), width = 800, height = 800)
  print(VlnPlot(object, c('Cdh5','Col1a1','Epcam','Ptprc'),ncol=2,pt.size=0.1))
  dev.off()
  png(paste(sample.name,'_vln_qc',tag,'.png',sep=""), width = 800, height = 800)
  print(VlnPlot(object, c('percent.mt','nCount_RNA','nFeature_RNA','perc.spliced','Mki67','Top2a'),ncol=2,pt.size=0.1))
  dev.off()
  
  ## Backtrack to original Filtration plots
  nclusters <- length(levels(object$seurat_clusters))
  # ViolinPlots broken apart by cluster
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Counts = ",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Features = ",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=(max.mt+0.5)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Max.Mt = ",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'),plot.caption = element_text(color='red'))
  png(paste(sample.name,'_filterback_1',tag,'.png',sep=""),width = 1500,height=500)
  print(plot_grid(p,q,r,ncol=3))
  dev.off()
  
  # ViolinPlots not broken apart, jitter points colored by cluster
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.counts, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Counts = ",min.counts)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0,log=T,group.by='orig.ident',cols=1) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Min.Features = ",min.features)) + theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0,group.by='orig.ident',cols=1,y.max=(max.mt+0.5)) +
    geom_jitter(aes(color = factor(object$seurat_clusters)), size = 0.5) +
    scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
    geom_hline(yintercept=max.mt, linetype='dashed', color='red', size=1) +
    labs(subtitle = paste("Max.Mt = ",max.mt)) +
    theme(legend.position = 'none',axis.title.x = element_blank(),plot.subtitle = element_text(color='red'))
  png(paste(sample.name,'_filterback_2',tag,'.png',sep=""),width = 1500,height=500)
  print(plot_grid(p,q,r,ncol=3))
  dev.off()
  
  # Scatter Plot colored by cluster
  list <- subset(pneum.metadata,Sample %in% sample.name)
  filtered <- subset(list,nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  
  png(paste(sample.name,'_filterback_3',tag,'.png',sep=""),width = 600,height=500)
  print(
    ggplot() +
      geom_point(data=filtered, aes(x=nCount_RNA, y=nFeature_RNA, color=factor(object$seurat_clusters)), size = 0.5) +
      scale_colour_manual(name="color", values=hue_pal()(nclusters)) +
      geom_vline(xintercept=min.counts, linetype='dashed', color='red', size=1) +
      geom_hline(yintercept=min.features, linetype='dashed', color='red', size=1) +
      labs(title = sample.name) +
      labs(subtitle = paste("Min.Counts = ",min.counts," / Min.Features = ",min.features," / Max.Mt = ",max.mt)) +
      theme(plot.title = element_text(size = 24, hjust = 0.5),
            plot.subtitle = element_text(color='red'),
            plot.caption = element_text(color='red'),
            axis.title.x = element_text(size = 20,color='black'),
            axis.title.y = element_text(size = 20,color='black'),
            axis.text.x = element_text(size=12,color='black'),
            axis.text.y = element_text(size = 12,color='black'),
            axis.line = element_line(size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
  )
  dev.off()
  
  # Ribosomal Characterization
  object[['percent.ribo']] <- PercentageFeatureSet(object, pattern = "^Rp")
  png(paste(sample.name,'_fp_ribo',tag,'.png',sep=""), width = 500, height = 500)
  print(FeaturePlot(object, c('percent.ribo'), label = T))
  dev.off()
  png(paste(sample.name,'_vln_ribo',tag,'.png',sep=""), width = 500, height = 500)
  print(VlnPlot(object, c('percent.ribo'),pt.size=0.1))
  dev.off()
  
  # Epithelial Characterization
  png(paste(sample.name,'_fp_c5',tag,'.png',sep=""), width = 800, height = 800)
  print(FeaturePlot(object, c('Sftpc','Sftpb','Defb4','Lyz2'), label = T, repel = T))
  dev.off()
  png(paste(sample.name,'_fp_epi',tag,'.png',sep=""), width = 800, height = 800)
  print(FeaturePlot(object, c('Krt5','Sftpc','Hopx','Aqp5'), label = T, repel = T))
  dev.off()
  
  # Pass output to outside for further exploration
  object <<- object # Passes this outside of the function
}
titrationUMAP <- function(object,min.counts,min.features,max.mt){
  object@meta.data$Passing <- 'NoPass'
  cells.Pass <- WhichCells(object,expression = nCount_RNA > min.counts & nFeature_RNA > min.features & percent.mt < max.mt)
  meta.data.pass <- data.frame(Barcode = cells.Pass,Passing = 'Pass')
  rownames(meta.data.pass) <- cells.Pass
  if(length(cells.Pass) < ncol(object)){
    cells.NoPass <- WhichCells(object,cells = cells.Pass,invert = T)
    meta.data.nopass <- data.frame(Barcode = cells.NoPass,Passing = 'NoPass')
    rownames(meta.data.nopass) <- cells.NoPass
  }else{meta.data.nopass <- data.frame()}
  meta.data <- rbind(meta.data.pass,meta.data.nopass)
  object <- AddMetaData(object,metadata = meta.data)
  tag <- paste('min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
  png(paste(sample.name,'_TitrationDimPlot.',tag,'.png',sep=""), width = 500, height = 500)
  print(DimPlot(object,group.by = 'Passing',cols = c('#0529F4','#F42F05'),shuffle = T)
        + labs(title = sample.name,subtitle = paste('nCount_RNA >',min.counts,'& nFeature_RNA >',min.features,'& percent.mt <',max.mt)))
  dev.off()
}

# Get Filenames total
file.names.total <- list.files(file.path("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Pneumonectomy_Project"))

# Select samples to process in this script
file.names <- file.names.total


# Get Sample Names
sample.names <- strsplit(file.names, split = "[.]")
for (i in 1:length(sample.names)){
  sample.names[[i]] <- sample.names[[i]][1]
}
sample.names <- unlist(sample.names)
sample.names <- gsub("-", "_", sample.names)


# Define File Paths
file.paths <- paste("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Pneumonectomy_Project",file.names,sep='/')
file.paths

# Load data tag and save
seurat.data <- LoadSeuratData(file.names = file.names,sample.names = sample.names)
names(seurat.data) <- sample.names
#save(seurat.data,file = paste("Native_w_Allie_2022-05-26.seurat.objs.Robj",sep=""))

# Load metadata tag and save
pneum.metadata <- LoadMetaData(seurat.data = seurat.data)
pneum.metadata$Sample <- factor(pneum.metadata$Sample,
                                levels=c('2_7_rat','2_8_rat','2_12_rat','2_14_rat','2_16_rat','2_28_rat','3_1_rat','8_10_rat','rat1A','rat811',
                                         'TE1',
                                         'fRat','mRat',
                                         'P0_f','P0_m','P3_f','P3_m','P7_f','P7_m','P14_f','P14_m','P28_m','P56_m',
                                          'TCO1','TCO2','TCO3'))
table(pneum.metadata$Sample)
save(pneum.metadata,file = paste("pneum.metadata",Sys.Date(),".Robj",sep=""))

#### Plot Full Dataset Violins Before Filtration ####
png(filename = 'PNEUM_nCount_RNA_Unfiltered.png',width = 20,height = 5,res = 200,units = 'in')
ggplot(pneum.metadata, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE,size=0)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0,stroke = 0,position = position_jitterdodge())+
  ggtitle('nCount_RNA')
dev.off()

png(filename = 'PNEUM_nFeature_RNA_Unfiltered.png',width = 20,height = 5,res = 200,units = 'in')
ggplot(pneum.metadata, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE,size=0)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0,stroke = 0,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')
dev.off()

png(filename = 'PNEUM_percent.mt_Unfiltered.png',width = 20,height = 5,res = 200,units = 'in')
ggplot(pneum.metadata, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE,size=0)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0,stroke = 0,position = position_jitterdodge())+
  ggtitle('percent.mt')
dev.off()



#### Individual QC Plots Before Filtration ####
# Our standard Violin Plots of nUMI, nFeature_RNA, percent.mt - Individual PNG & PDF All
myplots <- vector('list',length(seurat.data))
for (i in 1:length(seurat.data)) {
  message(names(seurat.data[i]))
  object <- seurat.data[[i]]
  p <- VlnPlot(object, features = c("nCount_RNA"),pt.size=0.1,log=T) +
    theme(legend.position = 'none',axis.title.x = element_blank())
  q <- VlnPlot(object, features = c("nFeature_RNA"),pt.size=0.1,log=T) +
    theme(legend.position = 'none',axis.title.x = element_blank())
  r <- VlnPlot(object, features = c("percent.mt"),pt.size=0.1,y.max=50) +
    theme(legend.position = 'none',axis.title.x = element_blank())
  myplots[[i]] <- plot_grid(p,q,r,ncol=3)
  png(paste(names(seurat.data[i]),'_qc_1.png',sep=""),width = 1500,height=500)
  print(myplots[[i]])
  dev.off()
}
pdf(file = paste('PNEUM_QC_First_Extraction.pdf',sep='_'),width=24,height=12)
print(myplots)
dev.off()


#### # QC Plots from Luecken 2019, Fig. 2 - Replicas built in R ####
# Panel D - nFeature_RNA vs nUMI, colored by percent.mt
samples <- levels(pneum.metadata$Sample)
myplots <- vector('list',length(samples))
for (i in 1:length(samples)) {
  message(samples[i])
  object <- subset(pneum.metadata,Sample %in% samples[i])
  p <- ggplot(object, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
    scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
    labs(title = samples[i]) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  myplots[[i]] <- p
  png(paste(samples[i],'_qc_2.png',sep=""),width = 600,height=500)
  print(p)
  dev.off()
}
pdf(file = paste('PNEUM_QC_First_Extraction_2.pdf',sep='_'),width=12,height=10)
print(myplots)
dev.off()


# All cells together (percent.mt)
png('PNEUM_qc_2.png',width = 700,height=500)
ggplot(pneum.metadata, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
  scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
  labs(title = 'All Project Samples') +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# Same as Panel D (nFeature_RNA vs nUMI), BUT colored by perc.spliced
samples <- levels(pneum.metadata$Sample)
myplots <- vector('list',length(samples))
for (i in 1:length(samples)) {
  message(samples[i])
  object <- subset(pneum.metadata,Sample %in% samples[i])
  p <- ggplot(object, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=perc.spliced)) +
    scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
    labs(title = samples[i]) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.title.x = element_text(size = 20,color='black'),
          axis.title.y = element_text(size = 20,color='black'),
          axis.text.x = element_text(size=12,color='black'),
          axis.text.y = element_text(size = 12,color='black'),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  myplots[[i]] <- p
  png(paste(samples[i],'_qc_3.png',sep=""),width = 600,height=500)
  print(p)
  dev.off()
}
pdf(file = paste('PNEUM_QC_First_Extraction_3.pdf',sep='_'),width=12,height=10)
print(myplots)
dev.off()

# All cells together (perc.spliced)
png('PNEUM_qc_3.png',width = 700,height=500)
ggplot(pneum.metadata, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=perc.spliced)) +
  scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
  labs(title = 'All Project Samples') +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# All cells together (colored by sample)
png('PNEUM_qc_4.png',width = 700,height=500)
ggplot(pneum.metadata[sample(nrow(pneum.metadata)),], aes(x=nCount_RNA, y=nFeature_RNA)) + 
  geom_point(aes(color=Sample)) +
  labs(title = 'All Project Samples') +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()



#### Titrate Filtration and Filter Each Object Individually ####
# Desired visualizations:
#   Same 3 plots w appropriate cutoffs (pre-cut)
#   < Subset >
#   (Start with not regressing on anything (i.e. cell cycle) first)
#   < NormalizeData, FindVariableFeatures, ScaleData, RunPCA >
#   Elbow plot, DimHeatmap
#   < FindNeighbors, FindClusters, RunUMAP >

#### TCO1 ####
sample.name <- 'TCO1'
# Low Start
min.counts <- 1200
min.features <- 800
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(1000, 10000, by = 1000)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(1000, 3000, by = 250)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 3000
min.features <- 2000
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)




#### TCO2 ####
sample.name <- 'TCO2'
# Low Start
min.counts <- 500
min.features <- 200
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(1000, 10000, by = 1000)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(500, 3000, by = 250)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 2000
min.features <- 1000
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(2000, 5000, by = 500)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(1000, 2000, by = 250)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(25, 10, by = -5)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 3000
min.features <- 2000
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### TCO3 ####
sample.name <- 'TCO3'
# Low Start
min.counts <- 500
min.features <- 200
max.mt <- 75
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(1000, 5000, by = 500)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(1000, 3000, by = 250)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 3000
min.features <- 2000
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)



#### mRat ####
sample.name <- 'mRat'
# Low Start
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 300, by = 50)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 250, by = 50)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(25, 5, by = -5)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 150
min.features <- 150
max.mt <- 15
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)




#### fRat ####
sample.name <- 'fRat'
# Low Start
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 300, by = 50)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 250, by = 50)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(25, 5, by = -5)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 150
min.features <- 150
max.mt <- 15
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### P0_m ####
sample.name <- 'P0_m'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(40, 10, by = -5)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 500
min.features <- 300
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)



#### P0_f ####
sample.name <- 'P0_f'
# Low Start
min.counts <- 400
min.features <- 300
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(400, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(300, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(40, 10, by = -5)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 700
min.features <- 300
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### P3_f ####
sample.name <- 'P3_f'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 500
min.features <- 300
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### P3_m ####
sample.name <- 'P3_m'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 800
min.features <- 500
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### P7_f ####
sample.name <- 'P7_f'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 800
min.features <- 500
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### P7_m ####
sample.name <- 'P7_m'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 800
min.features <- 500
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### P14_f ####
sample.name <- 'P14_f'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 800
min.features <- 500
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### P14_m ####
sample.name <- 'P14_m'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 500
min.features <- 300
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### P28_m ####
sample.name <- 'P28_m'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 700
min.features <- 300
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### P56_m ####
sample.name <- 'P56_m'
# Low Start
min.counts <- 200
min.features <- 200
max.mt <- 40
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 1000, by = 100)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 1000, by = 100)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 500
min.features <- 300
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### 2_12_rat ####
sample.name <- '2_12_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### 2_14_rat ####
sample.name <- '2_14_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 125
min.features <- 100
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### 2_16_rat ####
sample.name <- '2_16_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 125
min.features <- 100
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)


#### 2_28_rat ####
sample.name <- '2_28_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 125
min.features <- 100
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### 2_7_rat ####
sample.name <- '2_7_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 200
min.features <- 100
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### 2_8_rat ####
sample.name <- '2_8_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### 3_1_rat ####
sample.name <- '3_1_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### 8_10_rat ####
sample.name <- '8_10_rat'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### rat1A ####
sample.name <- 'rat1A'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 200
min.features <- 100
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### rat811 ####
sample.name <- 'rat811'
# Low Start
min.counts <- 50
min.features <- 50
max.mt <- 50
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(50, 200, by = 25)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(50, 200, by = 25)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 10, by = -10)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 125
min.features <- 100
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

#### TE1 ####
sample.name <- 'TE1'
# Low Start
min.counts <- 100
min.features <- 80
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)

# Titration
count.test <- seq(100, 500, by = 50)
for (i in 1:length(count.test)){
  titrationUMAP(object=object,min.counts=count.test[i],min.features=min.features,max.mt=max.mt)
}
feature.test <- seq(100, 500, by = 50)
for (i in 1:length(feature.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=feature.test[i],max.mt=max.mt)
}
mt.test <- seq(50, 5, by = -5)
for (i in 1:length(mt.test)){
  titrationUMAP(object=object,min.counts=min.counts,min.features=min.features,max.mt=mt.test[i])
}

# Compromise
min.counts <- 400
min.features <- 250
max.mt <- 25
tag <- paste('.min.ct',min.counts,'min.ft',min.features,'max.mt',max.mt,sep='.')
theProcess(sample.name=sample.name,min.counts=min.counts,min.features=min.features,max.mt=max.mt,tag=tag)




#### Look at filtered data all together ####

# Subset based on above cutoffs
data.sub <- seurat.data
data.sub$TCO1 <- subset(data.sub$TCO1,nCount_RNA > 3000 & nFeature_RNA > 2000 & percent.mt < 25)
data.sub$TCO2 <- subset(data.sub$TCO2,nCount_RNA > 3000 & nFeature_RNA > 2000 & percent.mt < 25)
data.sub$TCO3 <- subset(data.sub$TCO3,nCount_RNA > 3000 & nFeature_RNA > 2000 & percent.mt < 25)
data.sub$mRat <- subset(data.sub$mRat,nCount_RNA > 150 & nFeature_RNA > 150 & percent.mt < 15)
data.sub$fRat <- subset(data.sub$fRat,nCount_RNA > 150 & nFeature_RNA > 150 & percent.mt < 15)
data.sub$P0_m <- subset(data.sub$P0_m,nCount_RNA > 500 & nFeature_RNA > 300 & percent.mt < 25)
data.sub$P0_f <- subset(data.sub$P0_f,nCount_RNA > 700 & nFeature_RNA > 300 & percent.mt < 25)
data.sub$P3_f <- subset(data.sub$P3_f,nCount_RNA > 500 & nFeature_RNA > 300 & percent.mt < 25)
data.sub$P3_m <- subset(data.sub$P3_m,nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 25)
data.sub$P7_f <- subset(data.sub$P7_f,nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 25)
data.sub$P7_m <- subset(data.sub$P7_m,nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 25)
data.sub$P14_f <- subset(data.sub$P14_f,nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 25)
data.sub$P14_m <- subset(data.sub$P14_m,nCount_RNA > 500 & nFeature_RNA > 300 & percent.mt < 25)
data.sub$P28_m <- subset(data.sub$P28_m,nCount_RNA > 700 & nFeature_RNA > 300 & percent.mt < 25)
data.sub$P56_m <- subset(data.sub$P56_m,nCount_RNA > 500 & nFeature_RNA > 300 & percent.mt < 25)
data.sub$`2_12_rat` <- subset(data.sub$`2_12_rat`,nCount_RNA > 100 & nFeature_RNA > 80 & percent.mt < 25)
data.sub$`2_14_rat` <- subset(data.sub$`2_14_rat`,nCount_RNA > 125 & nFeature_RNA > 100 & percent.mt < 25)
data.sub$`2_16_rat` <- subset(data.sub$`2_16_rat`,nCount_RNA > 125 & nFeature_RNA > 100 & percent.mt < 25)
data.sub$`2_28_rat` <- subset(data.sub$`2_28_rat`,nCount_RNA > 125 & nFeature_RNA > 100 & percent.mt < 25)
data.sub$`2_7_rat` <- subset(data.sub$`2_7_rat`,nCount_RNA > 200 & nFeature_RNA > 100 & percent.mt < 25)
data.sub$`2_8_rat` <- subset(data.sub$`2_8_rat`,nCount_RNA > 100 & nFeature_RNA > 80 & percent.mt < 25)
data.sub$`3_1_rat` <- subset(data.sub$`3_1_rat`,nCount_RNA > 100 & nFeature_RNA > 80 & percent.mt < 25)
data.sub$`8_10_rat` <- subset(data.sub$`8_10_rat`,nCount_RNA > 100 & nFeature_RNA > 80 & percent.mt < 25)
data.sub$rat1A <- subset(data.sub$rat1A,nCount_RNA > 200 & nFeature_RNA > 100 & percent.mt < 25)
data.sub$rat811 <- subset(data.sub$rat811,nCount_RNA > 125 & nFeature_RNA > 100 & percent.mt < 25)
data.sub$TE1 <- subset(data.sub$TE1,nCount_RNA > 400 & nFeature_RNA > 250 & percent.mt < 15)
gc()

# Merge
data.sub <- merge(data.sub[[1]],data.sub[2:length(data.sub)])
gc()

# Organize metadata
data.sub$Sample <- factor(data.sub$Sample,
                                levels=c('2_7_rat','2_8_rat','2_12_rat','2_14_rat','2_16_rat','2_28_rat','3_1_rat','8_10_rat','rat1A','rat811',
                                         'TE1',
                                         'fRat','mRat',
                                         'P0_f','P0_m','P3_f','P3_m','P7_f','P7_m','P14_f','P14_m','P28_m','P56_m',
                                         'TCO1','TCO2','TCO3'))
# Inspect
png(filename = 'PNEUM_nCount_RNA_Filtered.png',width = 10,height = 5,res = 200,units = 'in')
VlnPlot(data.sub,features = 'nCount_RNA',group.by = 'Sample') +
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  #geom_point(size=0.01,stroke = 0,position = position_jitterdodge())+
  ggtitle('nCount_RNA')
dev.off()

png(filename = 'PNEUM_nFeature_RNA_Filtered.png',width = 10,height = 5,res = 200,units = 'in')
VlnPlot(data.sub,features = 'nFeature_RNA',group.by = 'Sample') +
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  #geom_point(size=0.2,stroke = 0,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')
dev.off()

png(filename = 'PNEUM_percent.mt_Filtered.png',width = 10,height = 5,res = 200,units = 'in')
VlnPlot(data.sub,features = 'percent.mt',group.by = 'Sample') +
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  #geom_point(size=0.2,stroke = 0,position = position_jitterdodge())+
  ggtitle('percent.mt')
dev.off()


# All cells together (perc.spliced)
png('PNEUM_qc_3_filtered.png',width = 700,height=500)
ggplot(data.sub@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=perc.spliced)) +
  scale_color_viridis(option="magma",direction=-1,limits = c(50, 100), oob = scales::squish) +
  labs(title = 'All Organoid Samples Filtered') +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# All cells together (colored by sample)
png('PNEUM_qc_4_filtered.png',width = 700,height=500)
ggplot(data.sub@meta.data[sample(nrow(data.sub@meta.data)),], aes(x=nCount_RNA, y=nFeature_RNA)) + 
  geom_point(aes(color=Sample)) +
  labs(title = 'All Organoid Samples Filtered') +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# All cells together (percent.mt)
png('PNEUM_qc_5_filtered.png',width = 700,height=500)
ggplot(data.sub@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(aes(color=percent.mt)) +
  scale_color_viridis(limits = c(0, 50), oob = scales::squish) +
  labs(title = 'All Organoid Samples Filtered') +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.title.x = element_text(size = 20,color='black'),
        axis.title.y = element_text(size = 20,color='black'),
        axis.text.x = element_text(size=12,color='black'),
        axis.text.y = element_text(size = 12,color='black'),
        axis.line = element_line(size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

# Save
rm(seurat.data)
gc()
save(data.sub, file = "pneum.project.data.sub.2022-08-18.Robj")

######################################
###### UNADAPTED LEGACY CODE BELOW ######
#### First Look All Project Data ####
pneum.combined <- data.sub
pneum.combined <- NormalizeData(pneum.combined)
pneum.combined <- ScaleData(pneum.combined)
pneum.combined <- FindVariableFeatures(pneum.combined)
pneum.combined <- RunPCA(pneum.combined, npcs = 100)
pdf(file='pneum.combined.PCs.pdf',width=10,height=8)
ElbowPlot(pneum.combined,ndims = 100)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=1:9)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=10:18)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=19:27)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=28:36)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=37:45)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=46:54)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=55:63)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=64:72)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=73:81)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=82:90)
PCHeatmap(pneum.combined,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
pneum.combined <- RunUMAP(pneum.combined, reduction = "pca", dims = 1:48)
pneum.combined <- FindNeighbors(pneum.combined, reduction = "pca", dims = 1:48)
pneum.combined <- FindClusters(pneum.combined, resolution = 0.2)
DimPlot(pneum.combined, reduction = "umap", label = TRUE, repel = TRUE)


# Visualization
png(file='pneum.combined.UMAP.first.look.png',width=12,height=6,units = 'in',res=300)
p1 <- DimPlot(pneum.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p2 <- DimPlot(pneum.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
dev.off()





##### Integration #####
# List
lung.list <- list(pneum.lung.v2,pneum.lung.v3)
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

# this command creates an 'integrated' data assay
lung.combined <- IntegrateData(anchorset = lung.anchors)
gc()

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(lung.combined) <- "integrated"

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

# Embed and cluster [[EDIT 2022-06-30 MAKING SURE LYMPHATICS ARE SEPARATE]]
lung.combined <- RunUMAP(lung.combined, reduction = "pca", dims = 1:63)
lung.combined <- FindNeighbors(lung.combined, reduction = "pca", dims = 1:63)
DefaultAssay(lung.combined) <- "integrated"
lung.combined <- FindClusters(lung.combined, resolution = 0.9)
DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(lung.combined,'Prox1',label=T)
FeaturePlot(lung.combined,c('Cdh5','Prox1'),label=T,cells = WhichCells(lung.combined,idents = '35'))

# Visualization
png(file='lung.combined.UMAP.first.look.png',width=18,height=6,units = 'in',res=300)
p1 <- DimPlot(lung.combined, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(lung.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 + p3
dev.off()

# Class annotations
DefaultAssay(lung.combined) <- 'RNA'
lung.combined <- ScaleData(lung.combined) # Not necessary unless plotting scaled values, but doing just in case
png(file='lung.combined.UMAP.class.png',width=10,height=10,units = 'in',res=300)
FeaturePlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,order = T)
dev.off()
png(file='lung.combined.VLN.class.png',width=20,height=6,units = 'in',res=300)
VlnPlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 2)
dev.off()

lung.combined <- RenameIdents(lung.combined,
                              '0'='Epithelium',
                              '1'='Endothelium',
                              '2'='Immune',
                              '3'='Immune',
                              '4'='Epithelium',
                              '5'='Epithelium',
                              '6'='Immune',
                              '7'='Immune',
                              '8'='Immune',
                              '9'='Mesenchyme',
                              '10'='Immune',
                              '11'='Multiplet',
                              '12'='Endothelium',
                              '13'='Epithelium',
                              '14'='Immune',
                              '15'='Mesenchyme',
                              '16'='Immune',
                              '17'='Immune',
                              '18'='Endothelium',
                              '19'='Epithelium',
                              '20'='Endothelium',
                              '21'='Immune',
                              '22'='Immune',
                              '23'='Endothelium',
                              '24'='Epithelium',
                              '25'='Mesenchyme',
                              '26'='Immune',
                              '27'='Immune',
                              '28'='Low Info',
                              '29'='Immune',
                              '30'='Immune',
                              '31'='Epithelium',
                              '32'='Epithelium',
                              '33'='Mesenchyme',
                              '34'='Immune',
                              '35'='Endothelium',
                              '36'='Immune',
                              '37'='Immune',
                              '38'='Immune',
                              '39'='Immune',
                              '40'='Endothelium',
                              '41'='RBC',
                              '42'='Immune'
)

# Check cluster 11 which looks like multiplets
pdf(file='multiplet.cluster.11.on.first.integration.embedding.pdf',width = 14,height = 10)
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '11'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '11'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '11'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(lung.combined,features = c('percent.mt'))
VlnPlot(lung.combined,features = c('nFeature_RNA'))
VlnPlot(lung.combined,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 16 which looks like low info
pdf(file='multiplet.cluster.16.on.first.integration.embedding.pdf',width = 14,height = 10)
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '16'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '16'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '16'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(lung.combined,features = c('percent.mt'))
VlnPlot(lung.combined,features = c('nFeature_RNA'))
VlnPlot(lung.combined,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 28 which looks like low info
pdf(file='low.info.cluster.28.on.first.integration.embedding.pdf',width = 14,height = 10)
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '28'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '28'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '28'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(lung.combined,features = c('percent.mt'))
VlnPlot(lung.combined,features = c('nFeature_RNA'))
VlnPlot(lung.combined,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 36 which looks like multiplets
pdf(file='multiplet.cluster.36.on.first.integration.embedding.pdf',width = 14,height = 10)
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '36'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '36'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '36'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(lung.combined,features = c('percent.mt'))
VlnPlot(lung.combined,features = c('nFeature_RNA'))
VlnPlot(lung.combined,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 41 which looks like RBCs
DefaultAssay(lung.combined) <- 'RNA'
Idents(lung.combined) <-lung.combined$integrated_snn_res.0.9
temp.mark <- FindMarkers(lung.combined,ident.1 = '41',only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
View(temp.mark)
pdf(file='RBC.cluster.41.on.first.integration.embedding.pdf',width = 14,height = 10)
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '41'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = '41'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(lung.combined,'Hbb')
VlnPlot(lung.combined,features = c('percent.mt'))
VlnPlot(lung.combined,features = c('nFeature_RNA'))
VlnPlot(lung.combined,features = c('nCount_RNA'),log=T)
dev.off()


# Stash new class labels
lung.combined$CellClass <- Idents(lung.combined)

# Record landing point
Idents(lung.combined) <- lung.combined$CellClass
pdf(file='lung.combined.UMAP.class.labeled.pdf',width=12,height=10)
DimPlot(lung.combined,group.by = 'CellClass')+ggtitle('Cell Class - Initial')
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Immune'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Epithelium'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Endothelium'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Mesenchyme'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Low Info'),features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA','percent.mt','perc.spliced')) 
FeaturePlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Multiplet'),features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
dev.off()

# save
rm(lung.anchors)
rm(lung.list)
gc()
save(lung.combined,file='lung.v2.v3.combined.classed.for.BEFM.2022-06-10.Robj')

#### Class Separation ####
# Load Data from previous save
load("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Native_Data_Cleaning_2022-06-01/lung.v2.v3.combined.classed.for.BEFM.2022-06-10.Robj")

# Break out each cell class
epi <- subset(lung.combined,idents = 'Epithelium')
endo <- subset(lung.combined,idents = 'Endothelium')
mes <- subset(lung.combined,idents = 'Mesenchyme')
imm <- subset(lung.combined,idents = 'Immune')
gc()

# Remove total object to save working space
rm(lung.combined)
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
epi <- RunUMAP(epi, reduction = "pca", dims = 1:36)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:36)
epi <- FindClusters(epi, resolution = 0.2)

# Visualization
pdf(file='epi.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(epi, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi) <- 'RNA'
p4 <- FeaturePlot(epi,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(epi,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

### Clean up ###
# Check cluster 7 which looks like multiplets
pdf(file='epi.cluster.7.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(epi,cells = WhichCells(epi,idents = '7'))
FeaturePlot(epi,cells = WhichCells(epi),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi,idents = '7'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
FeaturePlot(epi,cells = WhichCells(epi,idents = '7'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(epi,features = c('percent.mt'))
VlnPlot(epi,features = c('nFeature_RNA'))
VlnPlot(epi,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 9 which looks like multiplets
pdf(file='epi.cluster.9.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
VlnPlot(epi,features=c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi,idents = '9'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi,idents = '9'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(epi,features = c('percent.mt'))
VlnPlot(epi,features = c('nFeature_RNA'))
VlnPlot(epi,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 10 which looks like multiplets
pdf(file='epi.cluster.10.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
VlnPlot(epi,features=c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi,idents = '10'),features = c('Epcam','Col1a1','Cdh5','Ptprc'))
FeaturePlot(epi,cells = WhichCells(epi,idents = '10'),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'))
VlnPlot(epi,features = c('percent.mt'))
VlnPlot(epi,features = c('nFeature_RNA'))
VlnPlot(epi,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
epi.clean <- subset(epi,idents = c('7','9','10'),invert = T)

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
epi.clean <- RunUMAP(epi.clean, reduction = "pca", dims = 1:30)
epi.clean <- FindNeighbors(epi.clean, reduction = "pca", dims = 1:30)
epi.clean <- FindClusters(epi.clean, resolution = 0.4)

# Visualization
pdf(file='epi.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(epi.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi.clean) <- 'RNA'
p4 <- FeaturePlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

### Clean up round 2 ###
# Check cluster 11 which looks like multiplets
pdf(file='epi.cluster.11.multiplet.on.second.integration.embedding.pdf',width = 10,height = 8)
FeaturePlot(epi.clean,label=T,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(epi.clean,features = c('percent.mt'))
VlnPlot(epi.clean,features = c('nFeature_RNA'))
VlnPlot(epi.clean,features = c('nCount_RNA'),log=T)
dev.off()

# Check cluster 10 which looks like multiplets
pdf(file='epi.cluster.10.multiplet.on.second.integration.embedding.pdf',width = 10,height = 8)
FeaturePlot(epi.clean,label=T,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(epi.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(epi.clean,features = c('percent.mt'))
VlnPlot(epi.clean,features = c('nFeature_RNA'))
VlnPlot(epi.clean,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
epi.clean.2 <- subset(epi.clean,idents = c('10','11'),invert = T)

# Re-Visualize round 2
DefaultAssay(epi.clean.2) <- 'integrated' 
epi.clean.2 <- ScaleData(epi.clean.2)
epi.clean.2<- RunPCA(epi.clean.2, npcs = 100)
pdf(file='epi.clean.2.PCs.pdf',width=10,height=8)
ElbowPlot(epi.clean.2,ndims = 100)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=1:9)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=10:18)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=19:27)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=28:36)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=37:45)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=46:54)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=55:63)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=64:72)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=73:81)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=82:90)
PCHeatmap(epi.clean.2,cells=200,balanced=T,dims=91:99)
dev.off()

# Embed and cluster
epi.clean.2 <- RunUMAP(epi.clean.2, reduction = "pca", dims = 1:29)
epi.clean.2 <- FindNeighbors(epi.clean.2, reduction = "pca", dims = 1:15)
epi.clean.2 <- FindClusters(epi.clean.2, resolution = 0.4)

# Visualization
pdf(file='epi.clean.2.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(epi.clean.2, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi.clean.2, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi.clean.2, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi.clean.2) <- 'RNA'
p4 <- FeaturePlot(epi.clean.2,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(epi.clean.2,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Identify Cell Types
DefaultAssay(epi.clean.2) <- 'RNA'
epi.clean.2.mark <- FindAllMarkers(epi.clean.2,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
epi.clean.2.mark$ratio <- epi.clean.2.mark$pct.1/epi.clean.2.mark$pct.2
save(epi.clean.2.mark,file='epi.clean.2.mark.Robj')
epi.clean.2 <- RenameIdents(epi.clean.2,
                          '0'='ATII',
                          '1'='ATI',
                          '2'='ATI',
                          '3'='Tuft',
                          '4'='Ciliated',
                          '5'='ATII',
                          '6'='Secretory',
                          '7'='ATII-ATI',
                          '8'='BASC'
)
epi.clean.2$CellType <- Idents(epi.clean.2)

pdf(file='epi.clean.2.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(epi.clean.2, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(epi.clean.2, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(epi.clean.2, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(epi.clean.2) <- 'RNA'
p4 <- FeaturePlot(epi.clean.2,features = c('Ager','Sftpc','Ccdc153','Scgb1a1'),label = T,ncol = 3)
p5 <- VlnPlot(epi.clean.2,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
save(epi.clean.2,file='epi.clean.2.labeled.2022-07-23.Robj')
epi.clean.2.metadata <- epi.clean.2@meta.data
save(epi.clean.2.metadata,file = 'epi.clean.2.metadata.2022-07-23.Robj')
gc()
rm(epi)
rm(epi.clean.2.mark)
rm(epi.clean.2)
rm(epi.clean.2.metadata)
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
endo <- RunUMAP(endo, reduction = "pca", dims = 1:40)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:40)
endo <- FindClusters(endo, resolution = 0.2)

# Visualization
pdf(file='endo.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(endo, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo) <- 'RNA'
p4 <- FeaturePlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

### Clean up ###
# Check cluster 6 & 7 which looks like multiplets
pdf(file='endo.cluster.6.7.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(endo,cells = WhichCells(endo,idents = c('6','7')))
FeaturePlot(endo,cells = WhichCells(endo),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(endo,cells = WhichCells(endo),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
VlnPlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(endo,features = c('percent.mt'))
VlnPlot(endo,features = c('nFeature_RNA'))
VlnPlot(endo,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
endo.clean <- subset(endo,idents = c('6','7'),invert = T)

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
endo.clean <- RunUMAP(endo.clean, reduction = "pca", dims = 1:25)
endo.clean <- FindNeighbors(endo.clean, reduction = "pca", dims = 1:25)
endo.clean <- FindClusters(endo.clean, resolution = 0.4)

# Visualization
pdf(file='endo.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(endo.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo.clean) <- 'RNA'
p4 <- FeaturePlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()


# Identify Cell Types
DefaultAssay(endo.clean) <- 'RNA'
endo.clean.mark <- FindAllMarkers(endo.clean,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
endo.clean.mark$ratio <- endo.clean.mark$pct.1/endo.clean.mark$pct.2
save(endo.clean.mark,file='endo.clean.mark.Robj')
endo.clean <- RenameIdents(endo.clean,
                             '0'='gCap',
                             '1'='Venous',
                             '2'='Arterial',
                             '3'='Venous',
                             '4'='aCap',
                             '5'='Arterial',
                             '6'='gCap',
                             '7'='Lymphatic'
)
endo.clean$CellType <- Idents(endo.clean)

pdf(file='endo.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(endo.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo.clean) <- 'RNA'
p4 <- FeaturePlot(endo.clean,features = c('Ager','Sftpc','Ccdc153','Scgb1a1'),label = T,ncol = 3)
p5 <- VlnPlot(endo.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
save(endo.clean,file='endo.clean.labeled.2022-07-23.Robj')
endo.clean.metadata <- endo.clean@meta.data
save(endo.clean.metadata,file = 'endo.clean.metadata.2022-07-23.Robj')
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
mes <- RunUMAP(mes, reduction = "pca", dims = 1:30)
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:30)
DefaultAssay(mes) <- 'integrated'
mes <- FindClusters(mes, resolution = 1.0)

# Visualization
pdf(file='mes.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(mes, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes) <- 'RNA'
p4 <- FeaturePlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

### Clean up ###
# Check cluster 11 & 12 & 13 which looks like multiplets
pdf(file='mes.cluster.11.12.13.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(mes,cells = WhichCells(mes,idents = c('11','12','13')))
FeaturePlot(mes,cells = WhichCells(mes),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(mes,cells = WhichCells(mes),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
VlnPlot(mes,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(mes,features = c('percent.mt'))
VlnPlot(mes,features = c('nFeature_RNA'))
VlnPlot(mes,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
mes.clean <- subset(mes,idents = c('11','12','13'),invert = T)

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
mes.clean <- RunUMAP(mes.clean, reduction = "pca", dims = 1:30)
mes.clean <- FindNeighbors(mes.clean, reduction = "pca", dims = 1:30)
mes.clean <- FindClusters(mes.clean, resolution = 0.4)

# Visualization
pdf(file='mes.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(mes.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes.clean) <- 'RNA'
p4 <- FeaturePlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Identify Cell Types
DefaultAssay(mes.clean) <- 'RNA'
mes.clean.mark <- FindAllMarkers(mes.clean,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
mes.clean.mark$ratio <- mes.clean.mark$pct.1/mes.clean.mark$pct.2
save(mes.clean.mark,file='mes.clean.mark.Robj')
mes.clean <- RenameIdents(mes.clean,
                          '0'='Fibroblasts',
                          '1'='Smooth_Muscle',
                          '2'='Alveolar_Fibroblasts',
                          '3'='Myofibroblasts',
                          '4'='Pi16+_Fibroblasts',
                          '5'='Mesothelium',
                          '6'='Notum+',
                          '7'='Pericyte'
                          
)
mes.clean$CellType <- Idents(mes.clean)

pdf(file='mes.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(mes.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes.clean) <- 'RNA'
p4 <- FeaturePlot(mes.clean,features = c('Ager','Sftpc','Ccdc153','Scgb1a1'),label = T,ncol = 3)
p5 <- VlnPlot(mes.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
save(mes.clean,file='mes.clean.labeled.2022-07-23.Robj')
mes.clean.metadata <- mes.clean@meta.data
save(mes.clean.metadata,file = 'mes.clean.metadata.2022-07-23.Robj')
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
imm <- RunUMAP(imm, reduction = "pca", dims = 1:85)
imm <- FindNeighbors(imm, reduction = "pca", dims = 1:85)
DefaultAssay(imm) <- 'integrated'
imm <- FindClusters(imm, resolution = 0.6)

# Visualization
pdf(file='imm.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(imm, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(imm) <- 'RNA'
p4 <- FeaturePlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

### Clean up ###
# Check cluster 16 & 17 & 18 & 19 which look like multiplets
pdf(file='imm.cluster.16.17.18.19.multiplet.on.first.integration.embedding.pdf',width = 10,height = 8)
DimPlot(imm,cells = WhichCells(imm,idents = c('16','17','18','19')))
FeaturePlot(imm,cells = WhichCells(imm),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(imm,cells = WhichCells(imm,idents = c('16','17','18','19')),features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T)
FeaturePlot(imm,cells = WhichCells(imm),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
FeaturePlot(imm,cells = WhichCells(imm,idents = c('16','17','18','19')),features = c('nCount_RNA','nFeature_RNA','percent.mt','perc.spliced'),label = T)
VlnPlot(imm,features = c('Epcam','Col1a1','Cdh5','Ptprc'))
VlnPlot(imm,features = c('percent.mt'))
VlnPlot(imm,features = c('nFeature_RNA'))
VlnPlot(imm,features = c('nCount_RNA'),log=T)
dev.off()

# Clean data
imm.clean <- subset(imm,idents = c('16','17','19'),invert = T) # keep 18 for now, does not look low info or like a multiplet

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
imm.clean <- RunUMAP(imm.clean, reduction = "pca", dims = 1:85)
imm.clean <- FindNeighbors(imm.clean, reduction = "pca", dims = 1:85)
imm.clean <- FindClusters(imm.clean, resolution = 0.2)

# Visualization
pdf(file='imm.clean.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(imm.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(imm.clean) <- 'RNA'
p4 <- FeaturePlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Identify Cell Types
DefaultAssay(imm.clean) <- 'RNA'
imm.clean.mark <- FindAllMarkers(imm.clean,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
imm.clean.mark$ratio <- imm.clean.mark$pct.1/imm.clean.mark$pct.2
save(imm.clean.mark,file='imm.clean.mark.Robj')
imm.clean <- RenameIdents(imm.clean,
                          '0'='T',
                          '1'='B',
                          '2'='Mac Alv',
                          '3'='Mono 1',
                          '4'='Mono 2',
                          '5'='NK',
                          '6'='Neutro',
                          '7'='Mac Inter',
                          '8'='Cell cycle',
                          '9'='ILC',
                          '10'='Kdr+',
                          '11'='Killer T',
                          '12'='Mast',
                          '13'='C6+'
)
imm.clean$CellType <- Idents(imm.clean)

pdf(file='imm.clean.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(imm.clean, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm.clean, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm.clean, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(imm.clean) <- 'RNA'
p4 <- FeaturePlot(imm.clean,features = c('Ager','Sftpc','Ccdc153','Scgb1a1'),label = T,ncol = 3)
p5 <- VlnPlot(imm.clean,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
save(imm.clean,file='imm.clean.labeled.2022-07-23.Robj')
imm.clean.metadata <- imm.clean@meta.data
save(imm.clean.metadata,file = 'imm.clean.metadata.2022-07-23.Robj')
gc()
rm(imm)
rm(imm.clean.mark)
rm(imm.clean)
rm(imm.clean.metadata)
gc()

#### Merge Clean Classes Back Together ####
load('imm.clean.labeled.2022-07-23.Robj')
load('mes.clean.labeled.2022-07-23.Robj')
load('endo.clean.labeled.2022-07-23.Robj')
load('epi.clean.2.labeled.2022-07-23.Robj')
gc()
merge.clean <- merge(endo.clean,list(epi.clean.2,imm.clean,mes.clean))
rm(mes.clean)
rm(epi.clean.2)
rm(imm.clean)
rm(endo.clean)
gc()

# Re-Integrate
# List
lung.list <- SplitObject(merge.clean,split.by = 'Dataset')
#rm(merge.clean)
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

# this command creates an 'integrated' data assay [Modified to allow 'perfect' downstream cell cycle regression here]
# Load cell cycle genes from Katie
 cc.genes.rat <- readRDS('cc.genes.rat.rds')
# cc.genes.present <- intersect((rownames(x = GetAssayData(merge.clean,assay = "RNA"))),c(cc.genes.rat$s.genes,cc.genes.rat$g2m.genes))
# cc.genes.present
# features.to.integrate <- union(lung.anchors@anchor.features,cc.genes.present)
# features.to.integrate
# lung.combined <- IntegrateData(anchorset = lung.anchors, 
#                         dims = 1:30,
#                         features.to.integrate = features.to.integrate)

object <- IntegrateData(anchorset = lung.anchors)
gc()

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(lung.combined) <- "integrated"

# Embed the integrated object, regressing on cell cycle
# Assign cell cycle scores (on integrated slot, per https://github.com/satijalab/seurat/issues/2148)
lung.combined <- CellCycleScoring(lung.combined, 
                                  s.features = cc.genes.rat$s.genes, 
                                  g2m.features = cc.genes.rat$g2m.genes)
# Scale integrated slot, regressing on cell cycle genes
lung.combined <- ScaleData(lung.combined,vars.to.regress = c("S.Score", "G2M.Score"))
# Run the standard workflow for visualization and clustering
lung.combined <- RunPCA(lung.combined, npcs = 100)
pdf(file='lung.combined.regressed.PCs.pdf',width=10,height=8)
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
lung.combined <- RunUMAP(lung.combined, reduction = "pca", dims = 1:65)
lung.combined <- FindNeighbors(lung.combined, reduction = "pca", dims = 1:65,k.param = 5,force.recalc = T)
DefaultAssay(lung.combined) <- "integrated"
lung.combined <- FindClusters(lung.combined, resolution = 1.4)
DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)
lung.combined$stash <- Idents(lung.combined)

# Visualization
png(file='lung.clean.combined.regressed.UMAP.first.look.png',width=20,height=18,units = 'in',res=400)
p1 <- DimPlot(lung.combined, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(lung.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(lung.combined, reduction = "umap", label = TRUE, repel = TRUE)
p4 <- DimPlot(lung.combined, reduction = "umap",  group.by = "CellType",label = TRUE, repel = TRUE)
p1 + p2 + p3 + p4
dev.off()

# Show distribution of previouosly annotated 'Cell cycle' cluster
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = '26'))
Idents(lung.combined) <- lung.combined$CellType
DimPlot(lung.combined,cells = WhichCells(lung.combined,idents = 'Cell cycle'))
Idents(lung.combined) <- lung.combined$stash

# Save for later
save(lung.combined,file = 'lung.clean.combined.2022-07-23.Robj')

#### Annotate and Split by Class ####
# Load from save
load('lung.clean.combined.2022-07-23.Robj')

# Class annotations
DefaultAssay(lung.combined) <- 'RNA'
png(file='lung.clean.combined.UMAP.class.png',width=10,height=10,units = 'in',res=300)
FeaturePlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,order = T)
dev.off()
png(file='lung.clean.combined.VLN.class.png',width=20,height=6,units = 'in',res=300)
VlnPlot(lung.combined,features = c('Epcam','Col1a1','Cdh5','Ptprc'),ncol = 2)
dev.off()

lung.combined <- RenameIdents(lung.combined,
                              '0'='Immune',
                              '1'='Epithelium',
                              '2'='Immune',
                              '3'='Epithelium',
                              '4'='Epithelium',
                              '5'='Immune',
                              '6'='Endothelium',
                              '7'='Immune',
                              '8'='Immune',
                              '9'='Immune',
                              '10'='Endothelium',
                              '11'='Epithelium',
                              '12'='Endothelium',
                              '13'='Mesenchyme',
                              '14'='Immune',
                              '15'='Immune',
                              '16'='Mesenchyme',
                              '17'='Immune',
                              '18'='Endothelium',
                              '19'='Endothelium',
                              '20'='Epithelium',
                              '21'='Epithelium',
                              '22'='Epithelium',
                              '23'='Immune',
                              '24'='Endothelium',
                              '25'='Epithelium',
                              '26'='Immune',
                              '27'='Immune',
                              '28'='Immune',
                              '29'='Mesenchyme',
                              '30'='Mesenchyme',
                              '31'='Immune',
                              '32'='Immune',
                              '33'='Immune',
                              '34'='Epithelium',
                              '35'='Immune',
                              '36'='Mesenchyme',
                              '37'='Mesenchyme',
                              '38'='Endothelium',
                              '39'='Immune',
                              '40'='Immune',
                              '41'='Immune'
)


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
gc()
lung.combined$integrated_snn_res.1.2 <- NULL
lung.combined$integrated_snn_res.1 <- NULL
lung.combined$integrated_snn_res.0.9 <- NULL
lung.combined$integrated_snn_res.0.8 <- NULL
lung.combined$integrated_snn_res.0.6 <- NULL
lung.combined$integrated_snn_res.0.4 <- NULL
lung.combined$integrated_snn_res.0.3 <- NULL
lung.combined$integrated_snn_res.0.2 <- NULL
gc()
lung.combined.metadata <- lung.combined@meta.data
save(lung.combined.metadata,file='lung.clean.combined.classed.metadata.2022-07-24.Robj')

##### Epi Subclustering and Identification #####
epi <- subset(lung.combined,idents = 'Epithelium')
gc()
DefaultAssay(epi) <- 'integrated' # Let's try using the old integration-space
epi <- ScaleData(epi,vars.to.regress = c("S.Score", "G2M.Score"))
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
epi <- RunUMAP(epi, reduction = "pca", dims = 1:32)
epi <- FindNeighbors(epi, reduction = "pca", dims = 1:32)
DefaultAssay(epi) <- 'integrated'
epi <- FindClusters(epi, resolution = 0.4)

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
epi <- RenameIdents(epi,
                    '0'='ATII',
                    '1'='ATI',
                    '2'='ATI',
                    '3'='Tuft',
                    '4'='ATI',
                    '5'='Ciliated',
                    '6'='Secretory',
                    '7'='ATII',
                    '8'='ATII-ATI',
                    '9'='BASC'
)
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
epi.final.metadata <- epi@meta.data
save(epi.final.metadata,file = 'epi.final.metadata.2022-07-24.Robj')
gc()
rm(epi)
gc()

##### Endo Subclustering and Identification #####
endo <- subset(lung.combined,idents = 'Endothelium')
DefaultAssay(endo) <- 'integrated' # Let's try using the old integration-space
endo <- ScaleData(endo,vars.to.regress = c("S.Score", "G2M.Score"))
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
endo <- RunUMAP(endo, reduction = "pca", dims = 1:40)
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:40,force.recalc = T)
endo <- FindClusters(endo, resolution = 0.2)

# Visualization
pdf(file='endo.round2.UMAP.first.look.pdf',width=18,height=6)
p1 <- DimPlot(endo, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(endo, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(endo, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(endo) <- 'RNA'
p4 <- FeaturePlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc'),label = T,ncol = 3)
p5 <- VlnPlot(endo,features = c('Epcam','Col1a1','Cdh5','Ptprc','nCount_RNA','nFeature_RNA'),log=T,pt.size=0.1)
p1 + p2 + p3
p4
p5
dev.off()

# Identify Cell Types
DefaultAssay(endo) <- 'RNA'
endo.mark <- FindAllMarkers(endo,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
endo.mark$ratio <- endo.mark$pct.1/endo.mark$pct.2
save(endo.mark,file='endo.final.mark.Robj')
endo <- RenameIdents(endo,
                           '0'='gCap',
                           '1'='gCap',
                           '2'='aCap',
                           '3'='Arterial', #Efnb2
                           '4'='Venous', #Vcam1
                           '5'='Lymphatic'
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
endo.metadata <- endo@meta.data
save(endo.metadata,file = 'endo.final.metadata.2022-07-24.Robj')
gc()
rm(endo)
gc()


##### Mes Subclustering and Cleaning #####
mes <- subset(lung.combined,idents = 'Mesenchyme')
DefaultAssay(mes) <- 'integrated' # Let's try using the old integration-space
mes <- ScaleData(mes,vars.to.regress = c("S.Score", "G2M.Score"))
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
mes <- RunUMAP(mes, reduction = "pca", dims = 1:25)
mes <- FindNeighbors(mes, reduction = "pca", dims = 1:25)
DefaultAssay(mes) <- 'integrated'
mes <- FindClusters(mes, resolution = 1.5)
DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)

# Visualization
pdf(file='mes.round2.UMAP.first.look.pdf',width=18,height=6)
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


# Identify Cell Types
DefaultAssay(mes) <- 'RNA'
mes.mark <- FindAllMarkers(mes,only.pos = T,min.pct = 0.5,logfc.threshold = 0.5)
mes.mark$ratio <- mes.mark$pct.1/mes.mark$pct.2
save(mes.mark,file='mes.final.mark.Robj')
mes <- RenameIdents(mes,
                    '0'='Pi16+_Fibroblasts',
                    '1'='Mesothelium',
                    '2'='Itga8+_Fibroblasts',
                    '3'='Itga8+_Fibroblasts',
                    '4'='Myofibroblasts',
                    '5'='Myofibroblasts',
                    '6'='Itga8+_Fibroblasts',
                    '7'='Myofibroblasts',
                    '8'='Myofibroblasts',
                    '9'='Itga8+_Fibroblasts',
                    '10'='Notum+_Fibroblasts',
                    '11'='Smooth_Muscle',
                    '12'='Myofibroblasts',
                    '13'='Pericytes'
                          
)
mes$CellType_Final <- Idents(mes)

pdf(file='mes.final.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(mes, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(mes, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(mes, reduction = "umap", label = TRUE, repel = TRUE)
DefaultAssay(mes) <- 'RNA'
p4 <- FeaturePlot(mes,features = c('Col13a1','Col14a1','Gucy1a1','Acta2','Wnt11','Msln'),label = T,ncol = 3)
p5 <- VlnPlot(mes,features = c('Col13a1','Col14a1','Gucy1a1','Acta2','Wnt11','Msln'),pt.size=0.1)p1 + p2 + p3
p4
p5
dev.off()

# Save and Clean up
mes.metadata <- mes@meta.data
save(mes.metadata,file = 'mes.final.metadata.2022-07-24.Robj')
gc()
rm(mes)

##### Imm Subclustering and Cleaning #####
imm <- subset(lung.combined,idents = 'Immune')
DefaultAssay(imm) <- 'integrated' # Let's try using the old integration-space
imm <- ScaleData(imm,vars.to.regress = c("S.Score", "G2M.Score"))
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
imm <- RunUMAP(imm, reduction = "pca", dims = 1:85)
imm <- FindNeighbors(imm, reduction = "pca", dims = 1:85)
DefaultAssay(imm) <- 'integrated'
imm <- FindClusters(imm, resolution = 0.6)

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
save(imm.mark,file='imm.final.mark.Robj')
imm <- RenameIdents(imm,
                    '0'='B',
                    '1'='Mac Alv',
                    '2'='T',
                    '3'='Mono 1',
                    '4'='NK',
                    '5'='Mono 2',
                    '6'='Neutro',
                    '7'='Mac Inter',
                    '8'='Mono 3',
                    '9'='T',
                    '10'='Cell cycle',
                    '11'='ILC',
                    '12'='pDC',
                    '13'='Killer T',
                    '14'='T',
                    '15'='Mast',
                    '16'='Ccl24+'
)
imm$CellType_Final <- Idents(imm)

pdf(file='imm.Round2.UMAP.labeled.pdf',width=18,height=6)
p1 <- DimPlot(imm, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(imm, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(imm, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2 + p3
dev.off()

# Save and Clean up
imm.metadata <- imm@meta.data
save(imm.metadata,file = 'imm.final.metadata.2022-07-24.Robj')
rm(imm)
gc()

#### Map metadat back to integrated object ####
load("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Native_Data_Cleaning_2022-06-01/imm.final.metadata.2022-07-24.Robj")
load("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Native_Data_Cleaning_2022-06-01/mes.final.metadata.2022-07-24.Robj")
load("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Native_Data_Cleaning_2022-06-01/endo.final.metadata.2022-07-24.Robj")
load("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff/Native_Data_Cleaning_2022-06-01/epi.final.metadata.2022-07-24.Robj")

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

# Visualization
png(file='lung.combined.final.UMAP.annotated.png',width=20,height=18,units = 'in',res=400)
p1 <- DimPlot(lung.combined, reduction = "umap", group.by = "Dataset",shuffle = T)
p2 <- DimPlot(lung.combined, reduction = "umap", group.by = "Sample",shuffle = T)
p3 <- DimPlot(lung.combined, reduction = "umap", group.by = "CellClass_Final")
p4 <- DimPlot(lung.combined, reduction = "umap",  group.by = "CellType_Final",label = TRUE, repel = TRUE)
p1 + p2 + p3 + p4
dev.off()

# Final save
save(lung.combined,file='lung.combined.clean.classed.annotated.final.2022-07-24.Robj')



