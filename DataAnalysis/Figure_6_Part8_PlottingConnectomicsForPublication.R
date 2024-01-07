# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell")
# Set Seed
set.seed(123)
# Load Packages
require(Seurat)
require(RColorBrewer)
require(ggplot2)
require(cowplot)
library(tidyverse)
library(viridis)
require(ggsignif)
require(ggrepel)

# Load colors
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneumonectomy.col.pal.Robj")

# Connectomics
# Load NICHES data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/cell.to.cell.imputed.by.timepoint.Robj")

# Load connectomic significance table "combined.findings"
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/combined.findings.imputed.by.timepoint.2023-06-25.Robj")

# Load connectomic significance findings (thresholded) "combined.findings.niches.thresh"
#load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/combined.findings.niches.thresh.2023-06-28.Robj")
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/combined.findings.niches.thresh.2023-06-28.Robj")

# Format metadata for plotting
cell.to.cell.imputed.by.timepoint$Timepoint <- factor(cell.to.cell.imputed.by.timepoint$Timepoint.Sending,
                                                      levels = c('Day 0','Day 3','Day 7','Day 14'))
cell.to.cell.imputed.by.timepoint$Sample <- factor(cell.to.cell.imputed.by.timepoint$Sample.Sending,
                                                   levels = c('P0_f','P0_m','P3_f','P3_m',
                                                              'P7_f','P7_m','P14_f','P14_m'))
cell.to.cell.imputed.by.timepoint$Day <- NA
cell.to.cell.imputed.by.timepoint$Day[cell.to.cell.imputed.by.timepoint$Timepoint=='Day 0'] <- 0
cell.to.cell.imputed.by.timepoint$Day[cell.to.cell.imputed.by.timepoint$Timepoint=='Day 3'] <- 3
cell.to.cell.imputed.by.timepoint$Day[cell.to.cell.imputed.by.timepoint$Timepoint=='Day 7'] <- 7
cell.to.cell.imputed.by.timepoint$Day[cell.to.cell.imputed.by.timepoint$Timepoint=='Day 14'] <- 14

#### Cluster mechanisms at each timepoint???? ####
Idents(cell.to.cell.imputed.by.timepoint) <- cell.to.cell.imputed.by.timepoint$VectorType
day0 <- subset(cell.to.cell.imputed.by.timepoint,subset=Timepoint=='Day 0')
day3 <- subset(cell.to.cell.imputed.by.timepoint,subset=Timepoint=='Day 3')
day7 <- subset(cell.to.cell.imputed.by.timepoint,subset=Timepoint=='Day 7')
day14 <- subset(cell.to.cell.imputed.by.timepoint,subset=Timepoint=='Day 14')
day0 <- ScaleData(day0,features = rownames(day0))
day3 <- ScaleData(day3,features = rownames(day3))
day7 <- ScaleData(day7,features = rownames(day7))
day14 <- ScaleData(day14,features = rownames(day14))
day0.avg <- AverageExpression(day0,slot='scale.data',assays = 'CellToCell')
day3.avg <- AverageExpression(day3,slot='scale.data',assays = 'CellToCell')
day7.avg <- AverageExpression(day7,slot='scale.data',assays = 'CellToCell')
day14.avg <- AverageExpression(day14,slot='scale.data',assays = 'CellToCell')
day0.avg <- day0.avg$CellToCell
day3.avg <- day3.avg$CellToCell
day7.avg <- day7.avg$CellToCell
day14.avg <- day14.avg$CellToCell
day0.mech <- CreateSeuratObject(counts = t(day0.avg),assay = 'MechByVec')
day0.mech <- ScaleData(day0.mech,features = rownames(day0.mech))
day0.mech <- RunUMAP(day0.mech,features = rownames(day0.mech),slot = 'scale.data')
day0.mech <- FindNeighbors(day0.mech,features = rownames(day0.mech),dims = NULL)
day0.mech <- FindClusters(day0.mech,res=0.6)
DimPlot(day0.mech)


#### Signal fractions over time #####
cell.classes <- c('Epithelium','Endothelium','Mesenchyme','Immune')
cell.types <- unique(combined.findings.niches.thresh$SendingType)
cell.types <- cell.types[-30] # Remove Frzb+/Smoc1+ celltype, for now
timepoints <- c('Day 0','Day 3','Day 7','Day 14')

# What fraction is mesothelial?
timepoint.fractions.sending <- data.frame()
for(i in 1:length(timepoints)){
  temp <- combined.findings.niches.thresh[combined.findings.niches.thresh$Day==timepoints[i],]
  celltype.fractions <- list()
  for(j in 1:length(cell.types)){
    celltype.fractions[[j]] <- nrow(temp[temp$SendingType==cell.types[j],])/nrow(temp)
  }
  temp <- data.frame(Fraction = unlist(celltype.fractions)*100,
                     Day = timepoints[i],
                     CellType = cell.types)
  timepoint.fractions.sending <- rbind(timepoint.fractions.sending,temp)
}
timepoint.fractions.receiving <- data.frame()
for(i in 1:length(timepoints)){
  temp <- combined.findings.niches.thresh[combined.findings.niches.thresh$Day==timepoints[i],]
  celltype.fractions <- list()
  for(j in 1:length(cell.types)){
    celltype.fractions[[j]] <- nrow(temp[temp$ReceivingType==cell.types[j],])/nrow(temp)
  }
  temp <- data.frame(Fraction = unlist(celltype.fractions)*100,
                     Day = timepoints[i],
                     CellType = cell.types)
  timepoint.fractions.receiving <- rbind(timepoint.fractions.receiving,temp)
}
# Set up metadata ordering for fractional information
timepoint.fractions.sending$Day <- factor(timepoint.fractions.sending$Day,
                                          levels = c('Day 0','Day 3','Day 7','Day 14'))
timepoint.fractions.receiving$Day <- factor(timepoint.fractions.receiving$Day,
                                          levels = c('Day 0','Day 3','Day 7','Day 14'))

# Add labels for selected cells
cells.to.label <- c('Mesothelium','Tuft')
timepoint.fractions.sending$Label <- NA
timepoint.fractions.sending[timepoint.fractions.sending$CellType%in%cells.to.label,]$Label <- timepoint.fractions.sending[timepoint.fractions.sending$CellType%in%cells.to.label,]$CellType
timepoint.fractions.receiving$Label <- NA
timepoint.fractions.receiving[timepoint.fractions.receiving$CellType%in%cells.to.label,]$Label <- timepoint.fractions.receiving[timepoint.fractions.receiving$CellType%in%cells.to.label,]$CellType
# make plots
ggplot(data=timepoint.fractions.sending,
       aes(x=Day,y=Fraction,color=CellType,size=Fraction,label=CellType))+
  geom_point()+
  geom_label_repel(force        = 0.5,
                   nudge_x      = -0.25,
                   direction    = "y",
                   hjust        = 1,
                   segment.size = 0.2,
                   size=3,
                   max.overlaps = 10)+
  scale_color_manual(values = col.pal$CellType)+
  ggtitle('Fraction of Perturbed Signaling - Outgoing')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  NoLegend()

p1 <- ggplot(data=timepoint.fractions.sending,
             aes(x=Day,y=Fraction,color=CellType,size=Fraction,label=CellType))+
  geom_point()+
  geom_label_repel(force        = 0.5,
                   nudge_x      = -0.25,
                   direction    = "y",
                   hjust        = 1,
                   segment.size = 0.2,
                   size=3,
                   max.overlaps = 10)+
  scale_color_manual(values = col.pal$CellType)+
  ggtitle('Fraction of Perturbed Signaling - Outgoing')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  NoLegend()
p2 <- ggplot(data=timepoint.fractions.receiving,
             aes(x=Day,y=Fraction,color=CellType,size=Fraction,label=CellType))+
  geom_point()+
  geom_label_repel(force        = 0.5,
                   nudge_x      = -0.25,
                   direction    = "y",
                   hjust        = 1,
                   segment.size = 0.2,
                   size=3,
                   max.overlaps = 4)+
  scale_color_manual(values = col.pal$CellType)+
  ggtitle('Fraction of Perturbed Signaling - Incoming')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  NoLegend()

png(filename = 'Signaling.Fractions.over.time.png',width = 8,height = 12,units = 'in',res=300)
plot_grid(p1,p2,nrow=2)
dev.off()

#### Heatmaps showing biggest changes, both up and down, for certain signaling families across whole dataset for each timepoint ####

# Pull out significant findings for one signaling family
temp <- combined.findings.niches.thresh[grep('Wnt|Rspo',combined.findings.niches.thresh$Mechanism),]
# Break into positive and negative
temp.up <- temp[temp$delta=='Positive',]
temp.down <- temp[temp$delta=='Negative',]
# Summarize as number of unique ligands per vectorype per timepoint that are differentially perturbed
temp.up <- temp.up %>% group_by(Day,VectorType) %>% summarize(ligand.count = length(unique(Ligand)))
temp.down <- temp.down %>% group_by(Day,VectorType) %>% summarize(ligand.count = length(unique(Ligand)))

# Summarize as number of significant hits per vectortype per day
# temp.up <- as.data.frame(table(temp.up$VectorType,temp.up$Day))
# temp.down <- as.data.frame(table(temp.down$VectorType,temp.down$Day))
# Add Sending and receiving cell type metadata
temp.up$SendingType <- stringr::str_split_fixed(temp.up$VectorType,'—',n=2)[,1]
temp.up$ReceivingType <- stringr::str_split_fixed(temp.up$VectorType,'—',n=2)[,2]
temp.down$SendingType <- stringr::str_split_fixed(temp.down$VectorType,'—',n=2)[,1]
temp.down$ReceivingType <- stringr::str_split_fixed(temp.down$VectorType,'—',n=2)[,2]
# Add signaling family metadata
temp.up$Family <- 'WNT Family'
temp.down$Family <- 'WNT Family'

# Order cell type levels for plotting
temp.up$SendingType <- factor(temp.up$SendingType,levels = names(col.pal$CellType))
temp.up$ReceivingType <- factor(temp.up$ReceivingType,levels = names(col.pal$CellType))
temp.down$SendingType <- factor(temp.down$SendingType,levels = names(col.pal$CellType))
temp.down$ReceivingType <- factor(temp.down$ReceivingType,levels = names(col.pal$CellType))

# Fill in missing values
day0.up <- temp.up %>% subset(Day=='Day 0') %>% complete(SendingType,ReceivingType)
day3.up <- temp.up %>% subset(Day=='Day 3') %>% complete(SendingType,ReceivingType)
day7.up <- temp.up %>% subset(Day=='Day 7') %>% complete(SendingType,ReceivingType)
day14.up <- temp.up %>% subset(Day=='Day 14') %>% complete(SendingType,ReceivingType)
day0.down <- temp.down %>% subset(Day=='Day 0') %>% complete(SendingType,ReceivingType)
day3.down <- temp.down %>% subset(Day=='Day 3') %>% complete(SendingType,ReceivingType)
day7.down <- temp.down %>% subset(Day=='Day 7') %>% complete(SendingType,ReceivingType)
day14.down <- temp.down %>% subset(Day=='Day 14') %>% complete(SendingType,ReceivingType)

# Convert NA to zero
day0.up[is.na(day0.up$ligand.count),]$ligand.count <- 0
day3.up[is.na(day3.up$ligand.count),]$ligand.count <- 0
day7.up[is.na(day7.up$ligand.count),]$ligand.count <- 0
day14.up[is.na(day14.up$ligand.count),]$ligand.count <- 0
day0.down[is.na(day0.down$ligand.count),]$ligand.count <- 0
day3.down[is.na(day3.down$ligand.count),]$ligand.count <- 0
day7.down[is.na(day7.down$ligand.count),]$ligand.count <- 0
day14.down[is.na(day14.down$ligand.count),]$ligand.count <- 0

# Scale
day0.down$scale <- scale(day0.down$ligand.count)
day3.down$scale <- scale(day3.down$ligand.count)
day7.down$scale <- scale(day7.down$ligand.count)
day14.down$scale <- scale(day14.down$ligand.count)
day0.up$scale <- scale(day0.up$ligand.count)
day3.up$scale <- scale(day3.up$ligand.count)
day7.up$scale <- scale(day7.up$ligand.count)
day14.up$scale <- scale(day14.up$ligand.count)

# Define color limits
total.data.up <- rbind(day0.up,day3.up,day7.up,day14.up)
total.data.down <- rbind(day0.down,day3.down,day7.down,day14.down)
limits.use.up <- c(min(total.data.up$ligand.count),max(total.data.up$ligand.count))
limits.use.down <- c(min(total.data.down$ligand.count),max(total.data.down$ligand.count))

# Wnt signaling over time
p1 <- ggplot(day0.up,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top") +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.up)+NoLegend()+ggtitle('Day 0')
p2 <- ggplot(day3.up,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top")  +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.up)+NoLegend()+ggtitle('Day 3')
p3 <- ggplot(day7.up,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top")  +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.up)+NoLegend()+ggtitle('Day 7')
p4 <- ggplot(day14.up,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top")   +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.up)+ggtitle('Day 14')
p5 <- ggplot(day0.down,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top") +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.down)+NoLegend()+ggtitle('Day 0')
p6 <- ggplot(day3.down,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top")  +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.down)+NoLegend()+ggtitle('Day 3')
p7 <- ggplot(day7.down,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top")  +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.down)+NoLegend()+ggtitle('Day 7')
p8 <- ggplot(day14.down,
             aes(x=ReceivingType,y=SendingType,fill=ligand.count))+
  geom_tile()+scale_x_discrete(position = "top")   +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))+
  scale_fill_viridis_c(limits=limits.use.down)+ggtitle('Day 14')
png(filename = 'Wnt.Signaling.Sig.Num.Heatmaps.png',width = 22,height = 10,unit='in',res=300)
plot_grid(p1,p2,p3,p4,
          p5,p6,p7,p8,
          nrow=2,align = 'hv',axis = 'lr')
dev.off()

# Do not separate up from down...?
# Specify pattern
temp <- combined.findings.niches.thresh[grep('Wnt|Rspo',combined.findings.niches.thresh$Mechanism),] # Just a patterns worth
# Summarize as number of mechanisms that are differentially perturbed
#temp <- temp %>% group_by(Day,VectorType) %>% summarize(mech.count = length(unique(Mechanism)))
#temp$scaled.delta <- scale(temp$mean.log2FC)
temp <- temp %>% group_by(Day,VectorType) %>% summarize(mech.sum = sum(abs(mean.log2FC))) # as a sum of log fold change delta?

# Add Sending and receiving cell type metadata
temp$SendingType <- stringr::str_split_fixed(temp$VectorType,'—',n=2)[,1]
temp$ReceivingType <- stringr::str_split_fixed(temp$VectorType,'—',n=2)[,2]

# Order cell type levels for plotting
temp$SendingType <- factor(temp$SendingType,levels = names(col.pal$CellType))
temp$ReceivingType <- factor(temp$ReceivingType,levels = names(col.pal$CellType))

# Fill in missing values
day0 <- temp %>% subset(Day=='Day 0') %>% complete(SendingType,ReceivingType)
day3 <- temp %>% subset(Day=='Day 3') %>% complete(SendingType,ReceivingType)
day7 <- temp %>% subset(Day=='Day 7') %>% complete(SendingType,ReceivingType)
day14 <- temp %>% subset(Day=='Day 14') %>% complete(SendingType,ReceivingType)


# Convert NA to zero
day0[is.na(day0$mech.sum),]$mech.sum <- 0
day3[is.na(day3$mech.sum),]$mech.sum <- 0
day7[is.na(day7$mech.sum),]$mech.sum <- 0
day14[is.na(day14$mech.sum),]$mech.sum <- 0


# Scale
day0$scale <- scale(day0$mech.sum)
day3$scale <- scale(day3$mech.sum)
day7$scale <- scale(day7$mech.sum)
day14$scale <- scale(day14$mech.sum)

# Define color limits
total.data <- rbind(day0,day3,day7,day14)
total.data.down <- rbind(day0.down,day3.down,day7.down,day14.down)
limits.use <- c(min(total.data$mech.sum),max(total.data$mech.sum))

# Wnt signaling over time
p1 <- ggplot(day0,
             aes(x=ReceivingType,y=SendingType,fill=mech.sum))+
  geom_tile()+scale_x_discrete(position = "top") +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, hjust=0),axis.text.x.top = element_text(vjust = 0.5))+
  scale_fill_viridis_c()+ggtitle('Day 0')
p2 <- ggplot(day3,
             aes(x=ReceivingType,y=SendingType,fill=mech.sum))+
  geom_tile()+scale_x_discrete(position = "top")  +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, hjust=0),axis.text.x.top = element_text(vjust = 0.5))+
  scale_fill_viridis_c()+ggtitle('Day 3')
p3 <- ggplot(day7,
             aes(x=ReceivingType,y=SendingType,fill=mech.sum))+
  geom_tile()+scale_x_discrete(position = "top")  +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, hjust=0),axis.text.x.top = element_text(vjust = 0.5))+
  scale_fill_viridis_c()+ggtitle('Day 7')
p4 <- ggplot(day14,
             aes(x=ReceivingType,y=SendingType,fill=mech.sum))+
  geom_tile()+scale_x_discrete(position = "top")   +scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 90, hjust=0),axis.text.x.top = element_text(vjust = 0.5))+
  scale_fill_viridis_c()+ggtitle('Day 14')

png(filename = 'Wnt.Signaling.Sum.Delta.png',width = 22,height = 5,unit='in',res=300)
plot_grid(p1,p2,p3,p4,
          nrow=1,align = 'hv',axis = 'lr')
dev.off()

#### Volcano plots showing driving features to/from each cell type, by class, at each timepoint ####
cell.classes <- c('Epithelium','Endothelium','Mesenchyme','Immune')
cell.types <- unique(combined.findings.niches.thresh$SendingType)
cell.types <- cell.types[-30] # Remove Frzb+/Smoc1+ celltype, for now
timepoints <- c('Day 0','Day 3','Day 7','Day 14')


# Color work
colors <- col.pal$CellType

# Global parameters
p.thresh <- 0.0001
log2fc.thresh <- 0.25
cells.to.exclude <- c('Ciliated','Secretory','Frzb+/Smoc1+')
pattern.to.exclude <- c('Col|Itg|Hsp|Scgb|Sftp|Lama|Adam|Psen|Vim|Cdh')

# Influence arrangement
plot.list.celltypes <- list()
for(k in 1:length(cell.types)){
  plot.list.class <- list()
  for(i in 1:length(cell.classes)){
    plot.list.time <- list()
    for(j in 1:length(timepoints)){
      to.plot <- subset(combined.findings.niches.thresh,
                        combined.findings.niches.thresh$SendingType==cell.types[k] &
                          combined.findings.niches.thresh$CellClass.Receiving==cell.classes[i]&
                          combined.findings.niches.thresh$Day==timepoints[j])
      if(length(grep(pattern.to.exclude,to.plot$Mechanism))>0){
      to.plot <- to.plot[-grep(pattern.to.exclude,to.plot$Mechanism),] # this line removes mechanisms in advance to prevent cluttering the plot
      }
      if(nrow(to.plot)>0){
      to.plot$diffexpressed <- "NO"
      to.plot$diffexpressed[to.plot$mean.log2FC > log2fc.thresh & to.plot$comb.p.val.adj < p.thresh] <- "UP"
      to.plot$diffexpressed[to.plot$mean.log2FC < -log2fc.thresh & to.plot$comb.p.val.adj < p.thresh] <- "DOWN"
      to.exclude <- to.plot$Mechanism[grep(pattern.to.exclude,to.plot$Mechanism)]
      to.plot$delabel <- NA
      to.plot$delabel[to.plot$diffexpressed != "NO" & !(to.plot$Mechanism %in% to.exclude) & !(to.plot$ReceivingType %in% cells.to.exclude)] <- to.plot$Mechanism[to.plot$diffexpressed != "NO" & !(to.plot$Mechanism %in% to.exclude) & !(to.plot$ReceivingType %in% cells.to.exclude)]
      to.plot$ReceivingType <- factor(to.plot$ReceivingType,
                                      levels = names(colors))
      x.max <- max(to.plot$mean.log2FC)
      x.min <- min(to.plot$mean.log2FC)
      if(x.max>abs(x.min)){limit <- x.max}else{limit <- abs(x.min)}
      plot.list.time[[j]] <- ggplot(data=to.plot, aes(x=mean.log2FC, y=-log10(comb.p.val.adj), col=ReceivingType,label=delabel))+
                        geom_point()+
                        geom_text_repel(max.overlaps = 30)+
                        geom_vline(xintercept=c(-log2fc.thresh, log2fc.thresh), col="red") +
                        geom_hline(yintercept=-log10(p.thresh), col="red")+
                        xlim(-limit,limit)+
                        theme_minimal()+
                        ggtitle(paste(timepoints[j],cell.types[k],'to',cell.classes[i]))+
                        scale_color_manual(values = colors)
      }else(plot.list.time[[j]] <- NA)
    }
    names(plot.list.time) <- timepoints
    plot.list.class[[i]] <- plot.list.time
  }
  names(plot.list.class) <- cell.classes
  plot.list.celltypes[[k]] <- plot.list.class
}
names(plot.list.celltypes) <- cell.types
save(plot.list.celltypes,file = 'volcano.plot.list.celltypes.influence.Robj')

for(k in 1:length(cell.types)){
  print(k)
plot.list <- list(plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Epithelium']],nrow=1,align='hv'),
                              plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Endothelium']],nrow=1,align='hv'),
                              plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Mesenchyme']],nrow=1,align='hv'),
                              plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Immune']],nrow=1,align='hv'))
png(filename = paste(cell.types[k],'.volcano.tester.influence.png',sep=''),width=36,height = 24,units='in',res=300)
print(plot_grid(plotlist = plot.list,nrow=4,align = 'hv'))
dev.off()
}

# Niche arrangement
plot.list.celltypes <- list()
for(k in 1:length(cell.types)){
  plot.list.class <- list()
  for(i in 1:length(cell.classes)){
    plot.list.time <- list()
    for(j in 1:length(timepoints)){
      to.plot <- subset(combined.findings.niches.thresh,
                        combined.findings.niches.thresh$ReceivingType==cell.types[k] &
                          combined.findings.niches.thresh$CellClass.Sending==cell.classes[i]&
                          combined.findings.niches.thresh$Day==timepoints[j])
      if(length(grep(pattern.to.exclude,to.plot$Mechanism))>0){
        to.plot <- to.plot[-grep(pattern.to.exclude,to.plot$Mechanism),] # this line removes mechanisms in advance to prevent cluttering the plot
      }
      if(nrow(to.plot)>0){
        to.plot$diffexpressed <- "NO"
        to.plot$diffexpressed[to.plot$mean.log2FC > log2fc.thresh & to.plot$comb.p.val.adj < p.thresh] <- "UP"
        to.plot$diffexpressed[to.plot$mean.log2FC < -log2fc.thresh & to.plot$comb.p.val.adj < p.thresh] <- "DOWN"
        to.exclude <- to.plot$Mechanism[grep(pattern.to.exclude,to.plot$Mechanism)]
        to.plot$delabel <- NA
        to.plot$delabel[to.plot$diffexpressed != "NO" & !(to.plot$Mechanism %in% to.exclude) & !(to.plot$SendingType %in% cells.to.exclude)] <- to.plot$Mechanism[to.plot$diffexpressed != "NO" & !(to.plot$Mechanism %in% to.exclude) & !(to.plot$SendingType %in% cells.to.exclude)]
        to.plot$SendingType <- factor(to.plot$SendingType,
                                        levels = names(colors))
        x.max <- max(to.plot$mean.log2FC)
        x.min <- min(to.plot$mean.log2FC)
        if(x.max>abs(x.min)){limit <- x.max}else{limit <- abs(x.min)}
        plot.list.time[[j]] <- ggplot(data=to.plot, aes(x=mean.log2FC, y=-log10(comb.p.val.adj), col=SendingType,label=delabel))+
          geom_point()+
          geom_text_repel(max.overlaps = 30)+
          geom_vline(xintercept=c(-log2fc.thresh, log2fc.thresh), col="red") +
          geom_hline(yintercept=-log10(p.thresh), col="red")+
          xlim(-limit,limit)+
          theme_minimal()+
          ggtitle(paste(timepoints[j],cell.types[k],'Niche Signaling from',cell.classes[i]))+
          scale_color_manual(values = colors)
      }else(plot.list.time[[j]] <- NA)
    }
    names(plot.list.time) <- timepoints
    plot.list.class[[i]] <- plot.list.time
  }
  names(plot.list.class) <- cell.classes
  plot.list.celltypes[[k]] <- plot.list.class
}
names(plot.list.celltypes) <- cell.types
save(plot.list.celltypes,file = 'volcano.plot.list.celltypes.niche.Robj')

for(k in 1:length(cell.types)){
  print(k)
  plot.list <- list(plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Epithelium']],nrow=1,align='hv'),
                    plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Endothelium']],nrow=1,align='hv'),
                    plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Mesenchyme']],nrow=1,align='hv'),
                    plot_grid(plotlist = plot.list.celltypes[[cell.types[k]]][['Immune']],nrow=1,align='hv'))
  png(filename = paste(cell.types[k],'.volcano.tester.niche.png',sep=''),width=36,height = 24,units='in',res=300)
  print(plot_grid(plotlist = plot.list,nrow=4,align = 'hv'))
  dev.off()
}

# Mechanism arrangment (colored by sending celltype)
mech.list <- sort(unique(combined.findings.niches.thresh$Mechanism))
plot.list.mechanisms <- list()
for(k in 1:length(mech.list)){
  #plot.list.class <- list()
  #for(i in 1:length(cell.classes)){
    plot.list.time <- list()
    for(j in 1:length(timepoints)){
      to.plot <- subset(combined.findings.niches.thresh,
                        combined.findings.niches.thresh$Mechanism==mech.list[k] &
                          #combined.findings.niches.thresh$CellClass.Sending==cell.classes[i]&
                          combined.findings.niches.thresh$Day==timepoints[j])
      if(length(grep(pattern.to.exclude,to.plot$Mechanism))>0){
        to.plot <- to.plot[-grep(pattern.to.exclude,to.plot$Mechanism),] # this line removes mechanisms in advance to prevent cluttering the plot
      }
      if(nrow(to.plot)>0){
        to.plot$diffexpressed <- "NO"
        to.plot$diffexpressed[to.plot$mean.log2FC > log2fc.thresh & to.plot$comb.p.val.adj < p.thresh] <- "UP"
        to.plot$diffexpressed[to.plot$mean.log2FC < -log2fc.thresh & to.plot$comb.p.val.adj < p.thresh] <- "DOWN"
        to.exclude <- to.plot$Mechanism[grep(pattern.to.exclude,to.plot$Mechanism)]
        to.plot$delabel <- NA
        to.plot$delabel[to.plot$diffexpressed != "NO" & !(to.plot$Mechanism %in% to.exclude) & !(to.plot$SendingType %in% cells.to.exclude)] <- to.plot$VectorType[to.plot$diffexpressed != "NO" & !(to.plot$Mechanism %in% to.exclude) & !(to.plot$SendingType %in% cells.to.exclude)]
        to.plot$SendingType <- factor(to.plot$SendingType,
                                      levels = names(colors))
        x.max <- max(to.plot$mean.log2FC)
        x.min <- min(to.plot$mean.log2FC)
        if(x.max>abs(x.min)){limit <- x.max}else{limit <- abs(x.min)}
        plot.list.time[[j]] <- ggplot(data=to.plot, aes(x=mean.log2FC, y=-log10(comb.p.val.adj), col=ReceivingType,label=delabel))+
          geom_point()+
          geom_text_repel(max.overlaps = 30)+ # Could think about modifying to have two color labels? https://andrewwhitby.com/2017/09/18/multi-color-text-ggplot2/
          geom_vline(xintercept=c(-log2fc.thresh, log2fc.thresh), col="red") +
          geom_hline(yintercept=-log10(p.thresh), col="red")+
          xlim(-limit,limit)+
          theme_minimal()+
          ggtitle(paste(timepoints[j],mech.list[k],'Niche Signaling'))+
          scale_color_manual(values = colors)
        
      }else(plot.list.time[[j]] <- NA)
    }
    names(plot.list.time) <- timepoints
    #plot.list.class[[i]] <- plot.list.time
  #}
  #names(plot.list.class) <- cell.classes
  plot.list.mechanisms[[k]] <- plot.list.time
}
names(plot.list.mechanisms) <- mech.list
save(plot.list.mechanisms,file = 'volcano.plot.list.mechanisms.niche.Robj')

for(k in 1:length(mech.list)){
  print(k)
  plot <- plot_grid(plot.list.mechanisms[[mech.list[k]]][['Day 0']],
                              plot.list.mechanisms[[mech.list[k]]][['Day 3']],
                              plot.list.mechanisms[[mech.list[k]]][['Day 7']],
                              plot.list.mechanisms[[mech.list[k]]][['Day 14']],nrow=1,align='hv')
  png(filename = paste(mech.list[k],'.volcano.tester.receiving.color.png',sep=''),width=20,height = 4,units='in',res=300)
  print(plot)
  dev.off()
}


#### Trend Plotting ####
# Define TrendPlot plotting function
TrendPlot <- function(MOI,
                      emphasis.trends,
                      emphasis.colors,
                      logFC.thresh = 0.1,
                      cells.exclude = c('Secretory','Ciliated','Plasma','Frzb+/Smoc1+','Cell_cycle','Siglech+'),
                      buffer=0.1){
  
  meta.data <- cell.to.cell.imputed.by.timepoint@meta.data
  data <- data.frame(MOI = cell.to.cell.imputed.by.timepoint@assays$CellToCell@data[MOI,])
  data <- cbind(data,meta.data)
  data$vector.mech <- paste(data$VectorType,MOI)
  combined.findings.niches.thresh$vector.mech <- paste(combined.findings.niches.thresh$VectorType,
                                                      combined.findings.niches.thresh$Mechanism)
  # Limit to only significant trends
  data <- data[data$vector.mech%in%combined.findings.niches.thresh$vector.mech,]
  
  # Threshold based on mean log FC
  combined.findings.limited <- combined.findings.niches.thresh[abs(combined.findings.niches.thresh$mean.log2FC)>logFC.thresh,]
  
  # Exclude cells not of interest
  combined.findings.limited <- combined.findings.limited[!(combined.findings.limited$SendingType%in%cells.exclude),]
  combined.findings.limited <- combined.findings.limited[!(combined.findings.limited$ReceivingType%in%cells.exclude),]
  
  # Limit data to only limited trends
  limited.trends <- unique(data$VectorType)[unique(data$VectorType) %in% combined.findings.limited[combined.findings.limited$Mechanism==MOI,]$VectorType]
  data <- data[data$VectorType%in%limited.trends,]
  
  summary.out <- data%>%group_by(VectorType,Timepoint.Sending)%>%summarise(mean.connectivity=mean(MOI))
  summary.out$mech <- MOI
  summary.out <<- summary.out
  
  combined.findings.limited.MOI <<- combined.findings.limited[combined.findings.limited$Mechanism==MOI,]
  
  # Set y-axis limits based on the significant vectors
  #y.max <- max(data[data$VectorType%in%significant.trends,]$MOI)
  y.max <- max(summary.out$mean.connectivity)+buffer
  
  # Make everything grey to start
  vector.types <- unique(data$VectorType)
  line.colors <- rep('lightgrey',length(vector.types))
  names(line.colors) <- vector.types
  
  # Choose things to emphasize
  rankings <- data %>% group_by(VectorType,Day) %>% summarise(mean.MOI = mean(MOI))
  rankings <- rankings[order(rankings$mean.MOI,decreasing = T),]
  unique(rankings$VectorType[1:10])
  #emphasis.trends <- c('Mesothelium—Tuft','Mesothelium—BASC','Mesothelium—Mesothelium','Mesothelium—Arterial')
  
  # Emphasize certain vectors with color
  line.colors[emphasis.trends] <- emphasis.colors
  
  # Order the data to put significant lines on top
  data$VectorType <- factor(data$VectorType,
                            levels = c(vector.types[!(vector.types%in%emphasis.trends)],emphasis.trends))
  
  # Create plot
  to.plot <- ggplot(data=data,
                    aes(x=Day,y=MOI,group=VectorType,color=VectorType))+
    #geom_point()+
    geom_smooth(data = data,
                mapping = aes(group = VectorType,color=VectorType),method = "loess")+
    scale_color_manual(values = line.colors,breaks = emphasis.trends)+
    scale_x_continuous(breaks = c(0,3,7,14),expand=c(0,0))+
    labs(color='Trends of Interest')+
    ylim(0,y.max)+
    ggtitle(MOI)+
    ylab('Connectivity')+
    xlab(NULL)+
    theme_classic()+
    theme(legend.position = c(0.7, 0.9))
  
  png(filename = paste(MOI,'Trends.png'),width = 3.75,height = 3.75,units = 'in',res=300)
  print(to.plot)
  dev.off()
  return(to.plot)
  
}

# Wnt4—Fzd6
MOI <- 'Wnt4—Fzd6'
emphasis.trends <- c('Mesothelium—Tuft','Mesothelium—BASC','Mesothelium—Mesothelium','Mesothelium—Arterial')
emphasis.colors <- c("#021580","#A034F0",'brown','darkgreen')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Bmp5—Bmpr2
MOI <- 'Bmp5—Bmpr2'
emphasis.trends <- c('Mesothelium—Tuft','Col13a1_Fib—Tuft','Mesothelium—gCap','Col13a1_Fib—gCap')
emphasis.colors <- c("#021580","#A034F0",'brown','darkgreen')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Fgf18—Fgfr2
MOI <- 'Fgf18—Fgfr2'
emphasis.trends <- c('Mesothelium—Tuft','Col13a1_Fib—Tuft','Myofibroblasts—Tuft','Col13a1_Fib—gCap')
emphasis.colors <- c("#021580","#A034F0",'brown','darkgreen')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Rspo1—Lgr5
MOI <- 'Rspo1—Lgr5'
emphasis.trends <- c('Mesothelium—Tuft','Mac_Inter—Tuft','Mac_Inter—BASC')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Bmp4—Bmpr2
MOI <- 'Bmp4—Bmpr2'
emphasis.trends <- c('Mesothelium—gCap','Mesothelium—aCap')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Bmp4—Acvr2a
MOI <- 'Bmp4—Acvr2a'
emphasis.trends <- c('Mesothelium—Tuft','aCap—Tuft','Mesothelium—Myofibroblasts','Mesothelium—Col13a1_Fib')
emphasis.colors <- c("#021580",'brown','darkgreen',"#A034F0")
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Bmp4—Bmpr1a
MOI <- 'Bmp4—Bmpr1a'
emphasis.trends <- c('gCap—Tuft','Myofibroblasts—Tuft','Mesothelium—Myofibroblasts','aCap—Myofibroblasts')
emphasis.colors <- c("#021580",'brown','darkgreen',"#A034F0")
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=0.1)


# Angptl4—Tie1
MOI <- 'Angptl4—Tie1'
emphasis.trends <- c('Mesothelium—gCap','Mesothelium—aCap','Monocytes—aCap')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Igfbp4—Lrp6
MOI <- 'Igfbp4—Lrp6'
emphasis.trends <- c('Mesothelium—Tuft','Mesothelium—BASC')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer = 0.25)

# Cxcl6—Cxcr2
MOI <- 'Cxcl6—Cxcr2'
emphasis.trends <- c('Mesothelium—Neutrophils')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Il17b-Il17rb
MOI <- 'Il17b—Il17rb'
emphasis.trends <- c('Mesothelium—Tuft','Mesothelium—ILC')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=0.1)

# Wnt5a—Fzd1
MOI <- 'Wnt5a—Fzd1'
emphasis.trends <- c('Mesothelium—Tuft','Mesothelium—Lgr5_Fib','Myofibroblasts—Tuft','Col13a1_Fib—Lgr5_Fib')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Cthrc1—Fzd3
MOI <- 'Cthrc1—Fzd3'
emphasis.trends <- c('Myofibroblasts—Tuft','Lymphatic—Tuft','ATII—Mesothelium','ATII—Tuft')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=0.1)

# Bgn—Tlr4
MOI <- 'Bgn—Tlr4'
emphasis.trends <- c('Myofibroblasts—aCap','Mesothelium—Mac_Alv','Col14a1_Fib—aCap','Col14a1_Fib—Mac_Inter')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=2.75)

# Rgmb—Bmpr2
MOI <- 'Rgmb—Bmpr2'
emphasis.trends <- c('ATII—aCap','Mesothelium—Tuft','ATII—Myofibroblasts','ATII—gCap')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=0.75)

# Il6—Il6st
MOI <- 'Il6—Il6st'
emphasis.trends <- c('Col13a1_Fib—Lymphatic','Col13a1_Fib—Myofibroblasts','ILC—Venous','Mac_Inter—Venous')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Thbs1—Lrp5
MOI <- 'Thbs1—Lrp5'
emphasis.trends <- c('Mac_Inter—ATII','Mac_Inter—Tuft','Venous—Tuft','Venous—Myofibroblasts')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=0.5)

# Inha—Acvr1
MOI <- 'Inha—Acvr1'
emphasis.trends <- c('Tuft—Mesothelium','Tuft—Col13a1_Fib','Mesothelium—Mesothelium','Mesothelium—Myofibroblasts')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)
# Inhba—Acvr1b
MOI <- 'Inhba—Acvr1b'
emphasis.trends <- c('Lgr5_Fib—Tuft','Venous—Lgr5_Fib','Venous—aCap')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)
# Kitlg—Kit
MOI <- 'Kitlg—Kit'
emphasis.trends <- c('Mesothelium—gCap','Mesothelium—Venous')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)
# Il11—Il11ra
MOI <- 'Il11—Il11ra1'
emphasis.trends <- c('Mesothelium—Tuft','Mesothelium—Col13a1_Fib','Mesothelium—Lymphatic','Mesothelium—Mac_Alv')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors)

# Pgf—Flt1 
MOI <- 'Pgf—Flt1'
emphasis.trends <- c('Tuft—gCap','Lymphatic—aCap','Lymphatic—Venous','Myofibroblasts—Arterial')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Igf1—Igf1r 
MOI <- 'Igf1—Igf1r'
emphasis.trends <- c('Mac_Alv—aCap','Col13a1_Fib—aCap','Lymphatic—aCap','Myofibroblasts—aCap')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Gpc3—Cd81 https://pubmed.ncbi.nlm.nih.gov/23665349/
MOI <- 'Gpc3—Cd81'
emphasis.trends <- c('Col13a1_Fib—Col13a1_Fib','Col13a1_Fib—Mesothelium','BASC—Mesothelium','Mesothelium—ATII')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Fgf1—Fgfr2 
MOI <- 'Fgf1—Fgfr2'
emphasis.trends <- c('ATII—Tuft','ATI—ATII','ATII—ATII','ATI—Tuft')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Fgf9—Fgfr2 
MOI <- 'Fgf9—Fgfr2'
emphasis.trends <- c('Tuft—BASC','Tuft—ATI')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Wnt11—Fzd6 
MOI <- 'Wnt11—Fzd7'
emphasis.trends <- c('Mesothelium—Tuft')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

# Tgfb2—Tgfbr1 
MOI <- 'Tgfb2—Tgfbr1'
emphasis.trends <- c('Mac_Inter—Mac_Inter','Mac_Inter—Mesothelium','Mac_Inter—Mac_Alv','Mac_Inter—Lymphatic')
emphasis.colors <- c("#021580",'brown','darkgreen','#A034F0')
TrendPlot(MOI,emphasis.trends,emphasis.colors,buffer=1)

#### System to Cell Heatmap ####
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/cell.to.system.Robj")
VlnPlot(cell.to.system,'Rspo1—Lgr5',idents = 'Mesothelium',group.by = 'Timepoint')
cell.to.system$Timepoint <- factor(cell.to.system$Timepoint,levels=c('Day 0','Day 3','Day 7','Day 14'))
VlnPlot(cell.to.system,'Rspo1—Lgr5',idents = 'Mesothelium',group.by = 'Timepoint')
sub <- subset(cell.to.system,subset=CellType=='Mesothelium')
sub <- ScaleData(sub)
Idents(sub) <- sub[['Timepoint']]
sub.mark <- FindAllMarkers(sub,only.pos = F,logfc.threshold = 0.01)
sub.mark$ratio <- sub.mark$pct.1/sub.mark$pct.2

acceptible.mechs <- combined.findings.niches.thresh[combined.findings.niches.thresh$SendingType=='Mesothelium',]
acceptible.mechs <- unique(acceptible.mechs$Mechanism)
sub.mark.accept <- sub.mark[sub.mark$gene%in%acceptible.mechs,]
#MOI <- sub.mark.accept %>% group_by(cluster) %>% top_n(50,ratio)
png(filename = 'mesothelial.influence.timecourse.png',width = 6,height = 8,units = 'in',res=300)
DoHeatmap(sub,
          cells = sample(colnames(sub),ncol(sub)),
          features = c('Bmp4—Bmpr1a','Rspo1—Lgr5','Dll1—Notch1','Dkk2—Lrp6','Inhba—Acvr1b','Rgmb—Bmpr2','Csf1—Csf1r','Pgf—Nrp1',
                       'Wnt4—Fzd6','Igfbp4—Lrp6','Angptl4—Tie1','Cxcl6—Cxcr1','Cxcl6—Cxcr2','Il17b—Il17rb','Il11—Il11ra1',
                       'Bmp2—Bmpr2','Bmp2—Bmpr1a','Bmp5—Bmpr2','Fgf12—Fgfr2','Bgn—Tlr4','Fgf9—Fgfr2','Tgfb2—Tgfbr1','Kitlg—Kit',
                       'Wnt5a—Fzd1','Wnt5a—Lrp5','Fgf18—Fgfr1','Wnt11—Fzd7'),
          #features = MOI$gene,
          group.by = 'Timepoint',group.colors = col.pal$Timepoint) + scale_fill_gradient2(low = 'grey',mid = 'white',high = 'blue')

dev.off()

# Experimental embedding
male <- subset(cell.to.system,subset=Sample%in%c('P0_m','P3_m','P7_m','P14_m'))
female <- subset(cell.to.system,subset=Sample%in%c('P0_m','P3_m','P7_m','P14_m'),invert=T)
male <- FindVariableFeatures(male,nfeatures = 2000)
female <- FindVariableFeatures(female,nfeatures = 2000)
var.features.use <- intersect(female@assays$CellToSystem@var.features,male@assays$CellToSystem@var.features)
cell.to.system <- RunUMAP(cell.to.system,reduction = NULL,features = var.features.use)
DimPlot(cell.to.system,group.by='CellType',cols = col.pal$CellType)
VlnPlot(cell.to.system,'nFeature_CellToSystem')

Idents(male) <- male[['Sample']]
Idents(female) <- female[['Sample']]
features.male.3 <- FindMarkers(male,ident.1 = 'P3_m',ident.2 = 'P0_m',only.pos=F,min.pct = 0.04,logfc.threshold = 0.1)
features.female.3 <- FindMarkers(female,ident.1 = 'P3_f',ident.2 = 'P0_f',only.pos=F,min.pct = 0.04,logfc.threshold = 0.1)
features.male.7 <- FindMarkers(male,ident.1 = 'P7_m',ident.2 = 'P0_m',only.pos=F,min.pct = 0.04,logfc.threshold = 0.1)
features.female.7 <- FindMarkers(female,ident.1 = 'P7_f',ident.2 = 'P0_f',only.pos=F,min.pct = 0.04,logfc.threshold = 0.1)
features.male.14 <- FindMarkers(male,ident.1 = 'P14_m',ident.2 = 'P0_m',only.pos=F,min.pct = 0.04,logfc.threshold = 0.1)
features.female.14 <- FindMarkers(female,ident.1 = 'P14_f',ident.2 = 'P0_f',only.pos=F,min.pct = 0.04,logfc.threshold = 0.1)
var.features.use <- unique(c(intersect(rownames(features.male.3),rownames(features.female.3)),
                    intersect(rownames(features.male.7),rownames(features.female.7)),
                    intersect(rownames(features.male.14),rownames(features.female.14))))
cell.to.system <- RunUMAP(cell.to.system,reduction = NULL,features = var.features.use,min.dist = 0.8,n.neighbors = 40)
DimPlot(cell.to.system,group.by='CellType',cols = col.pal$CellType,pt.size=1)
mes <- subset(cell.to.system,subset=CellClass=='Mesenchyme')
mes <- FindVariableFeatures(mes)
mes <- RunUMAP(mes,reduction = NULL,features = mes@assays$CellToSystem@var.features)
DimPlot(mes,group.by='CellType',cols = col.pal$CellType,pt.size=1)

mes.list <- SplitObject(mes, split.by = "Sample")
mes.list <- lapply(X = mes.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "disp", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = mes.list)
mes.anchors <- FindIntegrationAnchors(object.list = mes.list, anchor.features = features)
mes.combined <- IntegrateData(anchorset = mes.anchors)
DefaultAssay(mes.combined) <- "integrated"
mes.combined <- ScaleData(mes.combined, verbose = FALSE)
mes.combined <- FindVariableFeatures(mes.combined,selection.method = 'disp')
mes.combined <- RunUMAP(mes.combined,reduction = NULL,features = mes.combined@assays$integrated@var.features)
DimPlot(mes.combined,group.by='CellType',cols = col.pal$CellType,pt.size=1)

# Cluster the mesenchyme influence just by mesenchymal ligands and tuft basc receptors
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/combined.findings.phenotype.2023-06-28.Robj")
combined.findings.phenotype <- combined.findings
receptors.consider <- combined.findings.phenotype[combined.findings.phenotype$CellType==c('Tuft','BASC'),]
receptors.consider <- receptors.consider[receptors.consider$consensus==TRUE,]
#receptors.consider <- receptors.consider[receptors.consider$comb.p.val.adj<0.05,]
ligands.consider <- combined.findings.phenotype[combined.findings.phenotype$CellType%in%c('Mesothelium','Col13a1_Fib','Col14a1_Fib','Frzb+/Smoc1+','Lgr5_Fib','Myofibroblasts','Pericytes'),]
ligands.consider <- ligands.consider[ligands.consider$consensus==TRUE,]
#ligands.consider <- ligands.consider[ligands.consider$comb.p.val.adj<0.05,]
var.features.use <- rownames(mes)
var.features.use <- var.features.use[str_split_fixed(var.features.use,pattern = '—',n=2)[,1] %in% ligands.consider$Gene]
var.features.use <- var.features.use[str_split_fixed(var.features.use,pattern = '—',n=2)[,2] %in% receptors.consider$Gene]
var.features.use
mes <- RunUMAP(mes,reduction = NULL,features = var.features.use)
DimPlot(mes,group.by='CellType',cols = col.pal$CellType,pt.size=1)
