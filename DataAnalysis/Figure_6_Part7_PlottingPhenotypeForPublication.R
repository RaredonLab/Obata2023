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

#### Phenotype ####
# Load data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneum.clean.annotated.2023-06-14.Robj")

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

# Format metadata for plotting
pneum.clean$Timepoint <- factor(pneum.clean$Timepoint,
                                levels = c('Day 0','Day 3','Day 7','Day 14'))
pneum.clean$Sample <- factor(pneum.clean$Sample,
                             levels = c('P0_f','P0_m','P3_f','P3_m',
                                        'P7_f','P7_m','P14_f','P14_m'))

# Define colors
col.pal <- list()
col.pal$Class <- c(brewer.pal(4,'Set1'))
names(col.pal$Class) <- c('Epithelium','Endothelium','Mesenchyme','Immune')
col.pal$Timepoint <- c(brewer.pal(4,'Dark2'))
names(col.pal$Timepoint) <- c('Day 0','Day 3','Day 7','Day 14')
col.pal$Sample <- c(brewer.pal(8,'Paired'))
names(col.pal$Sample) <- c('P0_f','P0_m','P3_f','P3_m',
                             'P7_f','P7_m','P14_f','P14_m')
col.pal$CellType <- c('firebrick','steelblue','springgreen','purple','salmon','skyblue','navyblue',
                                 'orangered','violetred','tomato','grey20','sandybrown',
                                 'saddlebrown','royalblue','plum4','lightgoldenrod','lawngreen','forestgreen','dimgray','deeppink',
                                 'red2','paleturquoise1','palevioletred','orchid4','purple4','plum1','olivedrab2',
                                 'slateblue','mediumvioletred','sienna','orange','seagreen',
                                 'lightseagreen','mediumpurple4')
names(col.pal$CellType) <- celltypes.ordered
#save(col.pal,file= 'pneumonectomy.col.pal.Robj')

# Plot UMAPS (total pneumonectomy timecourse)
p1 <- DimPlot(pneum.clean,group.by = 'CellType',cols = col.pal$CellType,pt.size = 0.75,shuffle = T)+
  ggtitle(NULL)+
 theme(plot.title=element_text(hjust=0))+
 theme(plot.title = element_text(size = 30, face = "bold"))+
 guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
 NoAxes()+NoLegend()+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p2 <- DimPlot(pneum.clean,group.by = 'CellClass',cols = col.pal$Class,shuffle = T)+
  ggtitle(NULL)+
 theme(plot.title=element_text(hjust=0))+
 theme(plot.title = element_text(size = 30, face = "bold"))+
 guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
 NoAxes()+NoLegend()+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p3 <- DimPlot(pneum.clean,group.by = 'Timepoint',cols = col.pal$Timepoint,shuffle = T)+
  ggtitle(NULL)+
 theme(plot.title=element_text(hjust=0))+
 theme(plot.title = element_text(size = 18, face = "bold"))+
 guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
 NoAxes()+NoLegend()+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p4 <- DimPlot(pneum.clean,group.by = 'Sample',cols = col.pal$Sample,shuffle = T)+
  ggtitle(NULL)+
 theme(plot.title=element_text(hjust=0))+
 theme(plot.title = element_text(size = 18, face = "bold"))+
 guides(color=guide_legend(ncol =1,override.aes = list(size=5)))+
 NoAxes()+NoLegend()+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
p5 <- plot_grid(p2,p3,p4,ncol = 1,align = 'hv')
png(filename = 'pneumonectomy.rat.lung.UMAPs.combined.png',width = 10,height = 0.7*10,res = 400,units = 'in')
plot_grid(p1,p5,rel_widths = c(2.75,1))
dev.off()

# Cell type ratios
epi <- droplevels(subset(pneum.clean,subset=CellClass=='Epithelium'))
end <- droplevels(subset(pneum.clean,subset=CellClass=='Endothelium'))
imm <- droplevels(subset(pneum.clean,subset=CellClass=='Immune'))
mes <- droplevels(subset(pneum.clean,subset=CellClass=='Mesenchyme'))
epi$CellType <- droplevels(epi$CellType)
end$CellType <- droplevels(end$CellType)
mes$CellType <- droplevels(mes$CellType)
imm$CellType <- droplevels(imm$CellType)

# Global class proportion
# Fraction of tissue classified as each Cell Class in each Timepoint
table <- as.data.frame(prop.table(table(pneum.clean$Timepoint,pneum.clean$CellClass),margin = 1))
p1 <- ggplot(data=table,
      aes(x=Var1,y=Freq,fill=Var2))+
         geom_col(position = 'fill')+
         scale_fill_manual(values=col.pal$Class)+
         theme_minimal()+
  ggtitle('Cell Class')+
  ylab('Cell Class Contribution to Tissue')+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Relative contribution of each Timepoint to each Cell Class
table <- as.data.frame(prop.table(prop.table(table(pneum.clean$Timepoint,pneum.clean$CellClass),margin = 1),margin = 2))
p2 <- ggplot(data=table,
       aes(x=Var2,y=Freq,fill=Var1))+
          geom_col(position = 'fill')+
          scale_fill_manual(values=col.pal$Timepoint)+
  ggtitle(NULL)+
  ylab('Timepoint Contribution to Class')+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Epithelial proportion
# Fraction of epithelium classified as each identity in each timepoint
table <- as.data.frame(prop.table(table(epi$Timepoint,epi$CellType),margin = 1))
p3 <- ggplot(data=table,
      aes(x=Var1,y=Freq,fill=Var2))+
        geom_col(position = 'fill')+
        scale_fill_manual(values=col.pal$CellType)+
  ggtitle('Epithelium')+
  ylab('Cell Type Contribution to Class')+
  xlab(NULL)+theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Relative contribution of each Timepoint to each identity
table <- as.data.frame(prop.table(prop.table(table(epi$Timepoint,epi$CellType),margin = 1),margin = 2))
p4 <- ggplot(data=table,
      aes(x=Var2,y=Freq,fill=Var1))+
        geom_col(position = 'fill')+
        scale_fill_manual(values=col.pal$Timepoint)+
  ggtitle(NULL)+
  ylab('Timepoint Contribution to Cell Type')+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Mesenchymal proportion
# Fraction of mesenchyme classified as each identity in each timepoint
table <- as.data.frame(prop.table(table(mes$Timepoint,mes$CellType),margin = 1))
p5 <- ggplot(data=table,
             aes(x=Var1,y=Freq,fill=Var2))+
                geom_col(position = 'fill')+
                scale_fill_manual(values=col.pal$CellType)+
  ggtitle('Mesenchyme')+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Relative contribution of each Timepoint to each identity
table <- as.data.frame(prop.table(prop.table(table(mes$Timepoint,mes$CellType),margin = 1),margin = 2))
p6 <- ggplot(data=table,
             aes(x=Var2,y=Freq,fill=Var1))+
          geom_col(position = 'fill')+
          scale_fill_manual(values=col.pal$Timepoint)+
  ggtitle(NULL)+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Immune proportion
# Fraction of immune cells classified as each identity in each timepoint
table <- as.data.frame(prop.table(table(imm$Timepoint,imm$CellType),margin = 1))
p7 <- ggplot(data=table,
             aes(x=Var1,y=Freq,fill=Var2))+
                geom_col(position = 'fill')+
                scale_fill_manual(values=col.pal$CellType)+
  ggtitle('Immune')+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Relative contribution of each Timepoint to each identity
table <- as.data.frame(prop.table(prop.table(table(imm$Timepoint,imm$CellType),margin = 1),margin = 2))
p8 <- ggplot(data=table,
             aes(x=Var2,y=Freq,fill=Var1))+
              geom_col(position = 'fill')+
              scale_fill_manual(values=col.pal$Timepoint)+
  ggtitle(NULL)+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Endothelial proportion
# Fraction of mesenchyme classified as each identity in each timepoint
table <- as.data.frame(prop.table(table(end$Timepoint,end$CellType),margin = 1))
p9 <- ggplot(data=table,
             aes(x=Var1,y=Freq,fill=Var2))+
                geom_col(position = 'fill')+
                scale_fill_manual(values=col.pal$CellType)+
  ggtitle('Endothelium')+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Relative contribution of each Timepoint to each identity
table <- as.data.frame(prop.table(prop.table(table(end$Timepoint,end$CellType),margin = 1),margin = 2))
p10 <- ggplot(data=table,
             aes(x=Var2,y=Freq,fill=Var1))+
                geom_col(position = 'fill')+
                scale_fill_manual(values=col.pal$Timepoint)+
  ggtitle(NULL)+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal()+
  NoLegend()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Assemble into a composite plot
png(filename = 'pneumonectomy.cell.type.ratios.png',width = 15,height = 7,units = 'in',res=300)
plot_grid(p1,p3,p5,p7,p9,p2,p4,p6,p8,p10,ncol = 5,align = 'hv')
dev.off()

# Create a composite legend
epi.types <- as.character(unique(epi$CellType))
end.types <- as.character(unique(end$CellType))
mes.types <- as.character(unique(mes$CellType))
imm.types <- as.character(unique(imm$CellType))
epi.types <- epi.types[epi.types!='Cell_cycle']
end.types <- end.types[end.types!='Cell_cycle']
mes.types <- mes.types[mes.types!='Cell_cycle']
imm.types <- imm.types[imm.types!='Cell_cycle']
cell.cycle <- 'Cell_cycle'

# Custom built legend
# Options are “bottomright”, “bottom”, “bottomleft”, “left”, “topleft”, “top”, “topright”, “right”, “center”
# Initial empty plot
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)

legend("topleft", 
       legend = names(col.pal$CellType[epi.types]), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$CellType[epi.types])

legend("left", 
       legend = names(col.pal$CellType[end.types]), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$CellType[end.types])


legend("bottomleft", 
       legend = names(col.pal$CellType[mes.types]), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$CellType[mes.types])


legend("top", 
       legend = names(col.pal$CellType[imm.types]), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$CellType[imm.types])

legend("bottom", 
       legend = names(col.pal$CellType[cell.cycle]), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$CellType[cell.cycle],)

legend("topright", 
       legend = names(col.pal$Class), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$Class,)

legend("right", 
       legend = names(col.pal$Timepoint), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$Timepoint)

legend("bottomright", 
       legend = names(col.pal$Sample), 
       pch=16, pt.cex=3, cex=1, bty='n',
       col = col.pal$Sample)

# Legend titles
# mtext("Sample", at=0.9, cex=2)
# mtext("Timepoint", at=0.9, cex=2)
# mtext("Cell Class", at=0.5, cex=2)
# mtext("Common State", at=0.075, cex=2)
# mtext("Immune Cell Types", at=0.075, cex=2)
# mtext("Mesenchymal Cell Types", at=0.075, cex=2)
# mtext("Endothelial Cell Types", at=0.075, cex=2)
# mtext("Epithelial Cell Types", at=0.075, cex=2)

# Store legend as a plot
legend.plot <- recordPlot()
png(filename = 'pnuemonectomy.legend.png',width = 8,height = 8,units='in',res=300)
legend.plot
dev.off()

# # Follow up questions Andre
# # Wt1?
# #https://febs.onlinelibrary.wiley.com/doi/epdf/10.1111/febs.16091
# 
# VlnPlot(pneum.clean,'Wt1',pt.size = 0)
# VlnPlot(pneum.clean,'Wt1',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Ccn2',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Yap1',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Mst1',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Tead1',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Tead2',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Tead3',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Tead4',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Brd4',split.by='Timepoint',pt.size = 0.1)
# 
# 
# VlnPlot(pneum.clean,'Nf2',split.by='Timepoint',pt.size = 0.1) # decrease of this would cause yap translocation to the nucleus: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9317057/
# VlnPlot(pneum.clean,'Ccn2',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Lats2',split.by='Timepoint',pt.size = 0.1) # decrease of this would cause yap translocation to the nucleus: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9317057/
# VlnPlot(pneum.clean,'Lats1',split.by='Timepoint',pt.size = 0.1) 
# VlnPlot(pneum.clean,'Igf1',split.by='Timepoint',pt.size = 0.1) # increase of this would promote yap translocation to nucleus
# VlnPlot(pneum.clean,'Pdk1',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Pik3ca',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Akt1',split.by='Timepoint',pt.size = 0.1)
# VlnPlot(pneum.clean,'Taok1',split.by='Timepoint',pt.size = 0.1) # decrease of this would cause yap translocation to the nucleus: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9317057/
# VlnPlot(pneum.clean,'Stat3',split.by='Timepoint',pt.size = 0.1)
# 
# VlnPlot(pneum.clean,'Strip1',split.by='Timepoint',pt.size = 0.1) #STRIPAK
# 
# mes <- subset(pneum.clean,idents = 'Mesothelium')
# Idents(mes) <- mes$Timepoint
# output <- data.frame(AverageExpression(mes,features = 'Wt1',assays = 'RNA'))
# View(output)
# 
# VlnPlot(mes,'Yap1',split.by='Timepoint',pt.size = 0.1)
# output <- data.frame(AverageExpression(mes,features = 'Yap1',assays = 'RNA'))
# output <- data.frame(AverageExpression(mes,features = 'Tead1',assays = 'RNA'))
# output <- data.frame(AverageExpression(mes,features = 'Nf2',assays = 'RNA')) # Merlin
# output <- data.frame(AverageExpression(mes,features = 'Ccn2',assays = 'RNA')) # Ctgf
# output <- data.frame(AverageExpression(mes,features = 'Taok1',assays = 'RNA')) # depletes if yaptaz activated
# output <- data.frame(AverageExpression(mes,features = 'Lats1',assays = 'RNA'))
# 
# View(output)
# 
# names <- c('Yap1','Tead1','Tead2','Tead3','Tead4','Taok1','Ccn2','Nf2','Lats2','Lats1','Igf1')
# hippo.on <- c(1,-1,-1,-1,-1,1,1,1,1,1,-1)
# names(hippo.on) <- names
# hippo.off <- -hippo.on
# 
# Idents(mes) <- 'Timepoint'
# timepoint.vals <- data.frame(AverageExpression(mes,features = rownames(hippo.stat),assays = 'RNA'))
# control.vals <- timepoint.vals[,1]
# timepoint.vals$Day0 <- timepoint.vals$RNA.Day.0 - timepoint.vals$RNA.Day.3
# timepoint.vals$Day3 <- timepoint.vals$RNA.Day.3 - timepoint.vals$RNA.Day.0
# timepoint.vals$Day7 <- timepoint.vals$RNA.Day.7 - timepoint.vals$RNA.Day.0
# timepoint.vals$Day14 <- timepoint.vals$RNA.Day.14 - timepoint.vals$RNA.Day.0
# View(timepoint.vals)
# timepoint.deltas <- timepoint.vals[,c('Day0','Day3','Day7','Day14')]
# timepoint.norms <- timepoint.vals[,c('RNA.Day.0','RNA.Day.3','RNA.Day.7','RNA.Day.14')]
# 
# ranges <- matrixStats::rowRanges(as.matrix(timepoint.norms))
# for (i in 1:length(hippo.on)){
#   if(hippo.on[i]<0){hippo.on[i]<-ranges[names(hippo.on[i]),1]}else{hippo.on[i]<-ranges[names(hippo.on[i]),2]}
# }
# 
# for (i in 1:length(hippo.off)){
#   if(hippo.off[i]<0){hippo.off[i]<-ranges[names(hippo.off[i]),1]}else{hippo.off[i]<-ranges[names(hippo.off[i]),2]}
# }
# 
# hippo.stat <- data.frame(hippo.on,hippo.off)
# 
# data <- cbind(hippo.stat,timepoint.norms)
# 
# View(data)
# 
# 
# # transform to a distance matrix
# d <- dist(t(data))
# # reduce dimension
# fit <- cmdscale(d,eig=TRUE, k=1) # k is the number of dim
# # plot solution
# ggplot(data.frame(fit$points),
#        aes(x=fit$points[,1],
#            y=1,
#            label=rownames(fit$points)))+
#   geom_point()+
#   geom_text_repel() # that is fucking cool. publish this it is dope.
# 
# # what if we try with all mesothelial cells?
# # define gene space
# names <- c('Yap1','Tead1','Tead2','Tead3','Tead4','Taok1',
#            'Ccn2','Nf2','Lats2','Lats1','Igf1')
# hippo.on <- c(1,-1,-1,-1,-1,1,
#               1,1,1,1,-1) 
# # 1 means it should up UP if HIPPO in ON (yap is in cytoplasm and is not actively promoting transcription 
# # -1 means it should be DOWN if HIPPO is ON (yap is in cytoplasm and is not actively promoting transcription)
# names(hippo.on) <- names
# hippo.off <- -hippo.on
# # pull data values for that gene space
# temp.data <- mes[names,]
# meta.data <- temp.data@meta.data
# # compute ranges
# ranges <- matrixStats::rowRanges(as.matrix(temp.data@assays$RNA@data))
# for (i in 1:length(hippo.on)){
#   if(hippo.on[i]<0){hippo.on[i]<-ranges[names(hippo.on[i]),1]}else{hippo.on[i]<-ranges[names(hippo.on[i]),2]}
# }
# 
# for (i in 1:length(hippo.off)){
#   if(hippo.off[i]<0){hippo.off[i]<-ranges[names(hippo.off[i]),1]}else{hippo.off[i]<-ranges[names(hippo.off[i]),2]}
# }
# 
# hippo.stat <- data.frame(hippo.on,hippo.off)
# 
# data <- cbind(hippo.stat,as.data.frame(temp.data@assays$RNA@data))
# 
# View(data)
# d <- dist(t(data))
# # reduce dimension
# fit <- cmdscale(d,eig=TRUE, k=1) # k is the number of dim
# # make plotting data frame
# to.plot <- as.data.frame(fit$points)
# to.plot$Sample <- NA
# to.plot[rownames(meta.data),]$Sample <- as.character(meta.data$orig.ident)
# to.plot['hippo.on',]$Sample <- 'Hippo-ON'
# to.plot['hippo.off',]$Sample <- 'Hippo-OFF'
# to.plot$Timepoint <- NA
# to.plot[rownames(meta.data),]$Timepoint <- as.character(meta.data$Timepoint)
# to.plot['hippo.on',]$Timepoint <- 'Hippo-ON'
# to.plot['hippo.off',]$Timepoint <- 'Hippo-OFF'
# 
# # plot solution
# to.plot$Timepoint <- factor(to.plot$Timepoint,levels = c('Day 0','Day 3','Day 7','Day 14'))
# ggplot(to.plot,
#        aes(x=V1,y=1,color=to.plot$Timepoint))+
#   geom_point()
# 
# ggplot(to.plot, 
#        aes(V1, fill = Timepoint, colour = Timepoint)) +
#   geom_density(alpha = 1) +
#   xlim(to.plot['hippo.off',]$V1,to.plot['hippo.on',]$V1)
# 
# # Using Small multiple
# ggplot(data=to.plot[-c(1:2),], 
#        aes(x=V1, group=Timepoint, fill=Timepoint)) +
#   geom_density(adjust=1.5) +
#   theme_classic() +
#   facet_wrap(~Timepoint) +
#   theme(
#     legend.position="none",
#     panel.spacing = unit(0.1, "lines"),
#     axis.ticks.x=element_blank()
#   )
# 
# to.plot$scale <- scale(to.plot$V1,center = T)[,1]
# 
# to.plot$normalized <- (to.plot$V1 - min(to.plot$V1)) / (max(to.plot$V1) - min(to.plot$V1))
# 
# plot.temp <- to.plot[-c(1:2),]
# 
# ggplot(data=plot.temp, 
#        aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
#   geom_density(alpha=0.5,color=NA)+
#   geom_jitter(aes(x=normalized,y=1))+
#   scale_fill_manual(values=col.pal$Timepoint)+
#   scale_color_manual(values=col.pal$Timepoint)+
#   theme_classic()+
#   xlab('Hippo State')+
#   ylab('Mesothelial Population Density')+
#   xlim(0,1)+ 
#   geom_vline(xintercept = 0, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 1, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 0.5, linetype="dotted", 
#              color = "black", size=0.5)+
#   facet_wrap(~Timepoint)
# 
# png(filename = 'Mesothelial HIPPO Status Densities.png',width = 8,height = 6,units = 'in',res = 300)
# ggplot(data=plot.temp, 
#        aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
#   geom_density(mapping = aes(group = Sample),alpha=0.5,color=NA)+
#   geom_jitter(aes(x=normalized,y=1))+
#   scale_fill_manual(values=col.pal$Timepoint)+
#   scale_color_manual(values=col.pal$Timepoint)+
#   theme_classic()+
#   xlab('Hippo State')+
#   ylab('Mesothelial Population Density')+
#   xlim(0,1)+ 
#   geom_vline(xintercept = 0, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 1, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 0.5, linetype="dotted", 
#              color = "black", size=0.5)+
#   facet_wrap(~Timepoint) # BEAUTIFUL
# dev.off()
# 
# p1 <- ggplot(data=subset(plot.temp,Timepoint=='Day 0'), 
#        aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
#   geom_density(alpha=0.5,color=NA)+
#   geom_jitter(aes(x=normalized,y=1))+
#   scale_fill_manual(values=col.pal$Timepoint)+
#   scale_color_manual(values=col.pal$Timepoint)+
#   theme_classic()+
#   xlab('Hippo State')+
#   ylab('Population Density')+
#   xlim(0,1)+ 
#   geom_vline(xintercept = 0, linetype="dotted", 
#                         color = "black", size=0.5)+
#   geom_vline(xintercept = 1, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 0.5, linetype="dotted", 
#            color = "black", size=0.5)
# p2 <- ggplot(data=subset(plot.temp,Timepoint=='Day 3'), 
#              aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
#   geom_density(alpha=0.5,color=NA)+
#   geom_jitter(aes(x=normalized,y=1))+
#   scale_fill_manual(values=col.pal$Timepoint)+
#   scale_color_manual(values=col.pal$Timepoint)+
#   theme_classic()+
#   xlab('Hippo State')+
#   ylab('Population Density')+
#   xlim(0,1)+ 
#   geom_vline(xintercept = 0, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 1, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 0.5, linetype="dotted", 
#              color = "black", size=0.5)
# p3 <- ggplot(data=subset(plot.temp,Timepoint=='Day 7'), 
#              aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
#   geom_density(alpha=0.5,color=NA)+
#   geom_jitter(aes(x=normalized,y=1))+
#   scale_fill_manual(values=col.pal$Timepoint)+
#   scale_color_manual(values=col.pal$Timepoint)+
#   theme_classic()+
#   xlab('Hippo State')+
#   ylab('Population Density')+
#   xlim(0,1)+ 
#   geom_vline(xintercept = 0, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 1, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 0.5, linetype="dotted", 
#              color = "black", size=0.5)
# p4 <- ggplot(data=subset(plot.temp,Timepoint=='Day 14'), 
#              aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
#   geom_density(alpha=0.5,color=NA)+
#   geom_jitter(aes(x=normalized,y=1))+
#   scale_fill_manual(values=col.pal$Timepoint)+
#   scale_color_manual(values=col.pal$Timepoint)+
#   theme_classic()+
#   xlab('Hippo State')+
#   ylab('Population Density')+
#   xlim(0,1)+ 
#   geom_vline(xintercept = 0, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 1, linetype="dotted", 
#              color = "black", size=0.5)+
#   geom_vline(xintercept = 0.5, linetype="dotted", 
#              color = "black", size=0.5)
# plot_grid(p1,p2,p3,p4)
# 
# # just playin
# temp <- as.matrix(pneum.clean@assays$RNA[pneum.clean@assays$RNA@var.features,])
# d <- dist(t(temp))
# # reduce dimension
# fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
