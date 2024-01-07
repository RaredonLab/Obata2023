# Set WD
setwd('/gpfs/loomis/scratch60/kaminski/mbr29/STARsolo-Output/DGEs_30000_Rank_Cutoff')
# OR
setwd("/Volumes/Samsung_T5/DGEs_30000_Rank_Cutoff")

# Packages
require(Seurat)
require(ggplot2)

# Get Filenames
file.names <- list.files(file.path("Pneumonectomy_Project"))
file.names

# Get Sample Names
sample.names <- strsplit(file.names, split = "[.]")
for (i in 1:length(sample.names)){
  sample.names[[i]] <- sample.names[[i]][1]
}
sample.names <- unlist(sample.names)
sample.names

# Define File Paths
file.paths <- paste("Pneumonectomy_Project",file.names,sep='/')
file.paths

# Load Metadata Function(don't retain data)
LoadMetaData <- function(file.paths,sample.names){
  data <- data.frame()
  for (i in 1:length(file.paths)){
    load(file = file.paths[i])
    seurat.data <- CreateSeuratObject(counts = output[['Gene']])
    seurat.data[['GeneFull']] <- CreateAssayObject(output[['GeneFull']])
    seurat.data[['spliced']] <- CreateAssayObject(output[['spliced']])
    seurat.data[['unspliced']] <- CreateAssayObject(output[['unspliced']])
    seurat.data[["percent.mt"]] <- PercentageFeatureSet(seurat.data, pattern = "^Mt-")
    seurat.data$perc.spliced <- 100*(colSums(seurat.data@assays$spliced)/(colSums(seurat.data@assays$spliced)+colSums(seurat.data@assays$unspliced)))
    seurat.data$Sample <- sample.names[i] # Tag with sample metadata
    data <- rbind(data,seurat.data@meta.data)
    rm(output)
    rm(seurat.data)
    gc()
  }
  return(data)
}

# Run function on complete dataset
pneum.metadata <- LoadMetaData(file.paths = file.paths,sample.names = sample.names)
# Order the sample levels
pneum.metadata$Sample <- factor(pneum.metadata$Sample,levels=c("2-12_rat", "2-14_rat", "2-16_rat" ,"2-28_rat" ,"2-7_rat" , "2-8_rat" , "3-1_rat" , "8-10_rat","rat1A" ,   "rat811", 
                                                                                 "TE1"   , "fRat"  ,   "mRat"    , "P0-f"  ,   "P0-m"  ,  "P3-f"  ,   "P3-m" ,    "P7-f"  ,   "P7-m"   ,"P14-f"   , "P14-m"  ,"P28-m"    ,   "P56-m" ))
save(pneum.metadata,file = 'pneum.metadata.Robj')

#### Plot Full Dataset Violins Before Filtration ####
png(filename = 'Pneum_nCount_RNA_Unfiltered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(pneum.metadata, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')
dev.off()

png(filename = 'Pneum_nFeature_RNA_Unfiltered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(pneum.metadata, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')
dev.off()

png(filename = 'Pneum_percent.mt_Unfiltered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(pneum.metadata, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')
dev.off()

#Initialize empty dataframe for total project filtered metadata
pneum.metadata.filtered <- data.frame()

#### Pre-Process DropSeq Data ####
dataset.members <- c(  "2-12_rat" ,"2-14_rat", "2-16_rat" ,
                      "2-28_rat", "2-7_rat" , "2-8_rat" , "3-1_rat" , "8-10_rat",
                      "rat1A"  ,  "rat811" ) 
dataset.name <- 'Pneum_DropSeq'

dataset <- subset(pneum.metadata,Sample %in% dataset.members)
dataset$Dataset <- dataset.name

# Propose thresholds
min.Counts <- 200
min.Features <- 200
max.Mt <- 10

# Plot QC metrics
pdf(file = paste(dataset.name,'QC_First_Extraction.pdf',sep='_'),width=12,height=8)
ggplot(dataset, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Implement Secondary Filtering
dataset.sub <- subset(dataset,percent.mt < max.Mt & nFeature_RNA > min.Features & nCount_RNA > min.Counts)

# Plot QC metrics after filtering
pdf(file = paste(dataset.name,'QC_After_Filtration.pdf',sep='_'),width=12,height=8)
ggplot(dataset.sub, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Store for later
pneum.metadata.filtered <- rbind(pneum.metadata.filtered,dataset.sub)

#### Pre-Process v2 Lung Data ####
dataset.members <- c( "fRat"   ,  "mRat"    ) 
dataset.name <- 'Pneum_Lung_v2'

dataset <- subset(pneum.metadata,Sample %in% dataset.members)
dataset$Dataset <- dataset.name

# Propose thresholds
min.Counts <- 200
min.Features <- 200
max.Mt <- 10

# Plot QC metrics
pdf(file = paste(dataset.name,'QC_First_Extraction.pdf',sep='_'),width=12,height=8)
ggplot(dataset, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Implement Secondary Filtering
dataset.sub <- subset(dataset,percent.mt < max.Mt & nFeature_RNA > min.Features & nCount_RNA > min.Counts)

# Plot QC metrics after filtering
pdf(file = paste(dataset.name,'QC_After_Filtration.pdf',sep='_'),width=12,height=8)
ggplot(dataset.sub, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Store for later
pneum.metadata.filtered <- rbind(pneum.metadata.filtered,dataset.sub)


#### Pre-Process v2 Trachea Data ####
dataset.members <- c( "TE1" ) 
dataset.name <- 'Pneum_Trachea_v2'

dataset <- subset(pneum.metadata,Sample %in% dataset.members)
dataset$Dataset <- dataset.name

# Propose thresholds
min.Counts <- 200
min.Features <- 200
max.Mt <- 10

# Plot QC metrics
pdf(file = paste(dataset.name,'QC_First_Extraction.pdf',sep='_'),width=12,height=8)
ggplot(dataset, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Implement Secondary Filtering
dataset.sub <- subset(dataset,percent.mt < max.Mt & nFeature_RNA > min.Features & nCount_RNA > min.Counts)

# Plot QC metrics after filtering
pdf(file = paste(dataset.name,'QC_After_Filtration.pdf',sep='_'),width=12,height=8)
ggplot(dataset.sub, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Store for later
pneum.metadata.filtered <- rbind(pneum.metadata.filtered,dataset.sub)


#### Pre-Process v3 Lung Data ####
dataset.members <- c( "P0-f"   ,  "P0-m"  ,   "P14-f"   , "P14-m"  ,  "P28-m"  ,  "P3-f"  ,   "P3-m"   ,  "P56-m"  ,  "P7-f"  ,   "P7-m"  ) 
dataset.name <- 'Pneum_Lung_v3'

dataset <- subset(pneum.metadata,Sample %in% dataset.members)
dataset$Dataset <- dataset.name

# Propose thresholds
min.Counts <- 1000
min.Features <- 500
max.Mt <- 20

# Plot QC metrics
pdf(file = paste(dataset.name,'QC_First_Extraction.pdf',sep='_'),width=12,height=8)
ggplot(dataset, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Implement Secondary Filtering
dataset.sub <- subset(dataset,percent.mt < max.Mt & nFeature_RNA > min.Features & nCount_RNA > min.Counts)

# Plot QC metrics after filtering
pdf(file = paste(dataset.name,'QC_After_Filtration.pdf',sep='_'),width=12,height=8)
ggplot(dataset.sub, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nCount_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')+geom_hline(yintercept=min.Counts, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=nFeature_RNA, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')+geom_hline(yintercept=min.Features, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
ggplot(dataset.sub, aes(x=Dataset, y=percent.mt, fill=Dataset)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')+geom_hline(yintercept=max.Mt, linetype='dashed', color='red', size=1)
dev.off()

# Store for later
pneum.metadata.filtered <- rbind(pneum.metadata.filtered,dataset.sub)

#### Plot Full Dataset Violins After Filtration ####

# Order the sample levels
pneum.metadata.filtered$Sample <- factor(pneum.metadata.filtered$Sample,levels=c("2-12_rat", "2-14_rat", "2-16_rat" ,"2-28_rat" ,"2-7_rat" , "2-8_rat" , "3-1_rat" , "8-10_rat","rat1A" ,   "rat811", 
                                                              "TE1"   , "fRat"  ,   "mRat"    , "P0-f"  ,   "P0-m"  ,  "P3-f"  ,   "P3-m" ,    "P7-f"  ,   "P7-m"   ,"P14-f"   , "P14-m"  ,"P28-m"    ,   "P56-m" ))


png(filename = 'Pneum_nCount_RNA_Filtered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(pneum.metadata.filtered, aes(x=Sample, y=nCount_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nCount_RNA')
dev.off()

png(filename = 'Pneum_nFeature_RNA_Filtered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(pneum.metadata.filtered, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_violin(trim=FALSE)+
  scale_y_log10()+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('nFeature_RNA')
dev.off()

png(filename = 'Pneum_percent.mt_Filtered.png',width = 20,height = 10,res = 200,units = 'in')
ggplot(pneum.metadata.filtered, aes(x=Sample, y=percent.mt, fill=Sample)) +
  geom_violin(trim=FALSE)+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(size=0.01,position = position_jitterdodge())+
  ggtitle('percent.mt')
dev.off()
