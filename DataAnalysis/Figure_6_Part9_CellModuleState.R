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

#### What is the status of the HIPPO pathway in Mesothelial cells? ####
                                 
mes <- subset(pneum.clean,idents = 'Mesothelium')
Idents(mes) <- mes$Timepoint
output <- data.frame(AverageExpression(mes,features = 'Yap1',assays = 'RNA'))
#View(output)
                
# Define a module
names <- c('Yap1','Tead1','Tead2','Tead3','Tead4','Taok1',
           'Ccn2','Nf2','Lats2','Lats1','Igf1')
# Define ON state
hippo.on <- c(1,-1,-1,-1,-1,1,
              1,1,1,1,-1) 
# 1 means it should up UP if HIPPO in ON (yap is in cytoplasm and is not actively promoting transcription 
# -1 means it should be DOWN if HIPPO is ON (yap is in cytoplasm and is not actively promoting transcription)
names(hippo.on) <- names
# Define OFF state
hippo.off <- -hippo.on

# Gather data for this module from a target population (Mesothelium, here)
temp.data <- mes[names,]
meta.data <- temp.data@meta.data

# Compute ranges for each feature
ranges <- matrixStats::rowRanges(as.matrix(temp.data@assays$RNA@data))

# Scale the ON/OFF states to the extrema of these ranges for each features
for (i in 1:length(hippo.on)){
  if(hippo.on[i]<0){hippo.on[i]<-ranges[names(hippo.on[i]),1]}else{hippo.on[i]<-ranges[names(hippo.on[i]),2]}
}
for (i in 1:length(hippo.off)){
  if(hippo.off[i]<0){hippo.off[i]<-ranges[names(hippo.off[i]),1]}else{hippo.off[i]<-ranges[names(hippo.off[i]),2]}
}
# Bind ON and OFF states
hippo.stat <- data.frame(hippo.on,hippo.off)

# Bind ON and OFF states to cellular states
data <- cbind(hippo.stat,as.data.frame(temp.data@assays$RNA@data))
#View(data)

# ?? Scale these values by row so that each distance is computed within its own space??
#data.scaled <- t(scale(t(data)))

# Transform into a (cell-cell) distance matrix
d <- dist(t(data))
#d <- dist(t(data.scaled))

# Reduce dimension to 1 dimension (to visualize on a number line)
fit <- cmdscale(d,eig=TRUE, k=1) # k is the number of dim

# Organize data for plotting
to.plot <- as.data.frame(fit$points)
to.plot$Sample <- NA
to.plot[rownames(meta.data),]$Sample <- as.character(meta.data$orig.ident)
to.plot['hippo.on',]$Sample <- 'Hippo-ON'
to.plot['hippo.off',]$Sample <- 'Hippo-OFF'
to.plot$Timepoint <- NA
to.plot[rownames(meta.data),]$Timepoint <- as.character(meta.data$Timepoint)
to.plot['hippo.on',]$Timepoint <- 'Hippo-ON'
to.plot['hippo.off',]$Timepoint <- 'Hippo-OFF'
# Organize metadata for plotting
to.plot$Timepoint <- factor(to.plot$Timepoint,levels = c('Day 0','Day 3','Day 7','Day 14'))

# Scale
to.plot$scale <- scale(to.plot$V1,center = T)[,1]
# Normalize
to.plot$normalized <- (to.plot$V1 - min(to.plot$V1)) / (max(to.plot$V1) - min(to.plot$V1))
# Isolate columns of interest
plot.temp <- to.plot[-c(1:2),]

# Fractions on both sides of the mean
frac.on.0 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 0',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 0',]),digits = 3)),'%',sep='')
frac.off.0 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 0',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 0',]),digits = 3)),'%',sep='')
frac.on.3 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 3',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 3',]),digits = 3)),'%',sep='')
frac.off.3 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 3',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 3',]),digits = 3)),'%',sep='')
frac.on.7 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 7',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 7',]),digits = 3)),'%',sep='')
frac.off.7 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 7',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 7',]),digits = 3)),'%',sep='')
frac.on.14 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 14',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 14',]),digits = 3)),'%',sep='')
frac.off.14 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 14',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 14',]),digits = 3)),'%',sep='')
dat_text <- data.frame(
  frac.on = c(frac.on.0, frac.on.3, frac.on.7,frac.on.14),
  frac.off = c(frac.off.0, frac.off.3, frac.off.7,frac.off.14),
  Timepoint   = as.factor(c('Day 0','Day 3','Day 7','Day 14'))
)

# Make base plot
plot.total <- ggplot(data=plot.temp, 
                    aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
                    geom_density(mapping = aes(group = Sample),alpha=0.5,color=NA)+
                    geom_jitter(aes(x=normalized,y=1))+
                    scale_fill_manual(values=col.pal$Timepoint)+
                    scale_color_manual(values=col.pal$Timepoint)+
                    theme_classic()+
                    xlab('Hippo State')+
                    ylab('Mesothelial Cell Population Density')+
                    xlim(0,1)+ 
                    geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.5)+
                    geom_vline(xintercept = 1, linetype="dotted", color = "black", size=0.5)+
                    # geom_vline(xintercept = 0.5, linetype="dotted", 
                    #            color = "black", size=0.5)+
                    geom_vline(xintercept = mean(plot.temp$normalized), linetype="dotted", 
                               color = "black", size=0.5)+
                    facet_wrap(~Timepoint,ncol=1)+NoLegend()
plot.total
# Add annotations
plot.total <- plot.total + geom_text(
                    data    = dat_text,
                    mapping = aes(x = 0, y = 3, label = frac.off),
                    size=3.5,
                    hjust   = -0.1,
                    vjust   = -1
                  )+ geom_text(
                    data    = dat_text,
                    mapping = aes(x = 0.89, y = 3, label = frac.on),
                    size=3.5,
                    hjust   = -0.1,
                    vjust   = -1
                  )
plot.total
# Make final plot
png(filename = 'Mesothelial HIPPO Pathway Status Densities VERTICAL.png',width = 5.5,height = 6.5,units = 'in',res = 300)
plot.total
dev.off()

# Experiment: a plot of the actual network ON vs OFF state?
fit.on <- cmdscale(dist(data$hippo.on),eig = T,k=2)
to.plot.on<- data.frame(fit.on$points)
to.plot.on$Expression <- data$hippo.on
ggplot(to.plot.on,aes(y=to.plot.on$X1,
                       x=to.plot.on$X2,
                       color = to.plot.on$Expression,
                       label=rownames(data)))+
  geom_point()+
  ggrepel::geom_text_repel()+
  theme_classic()

fit.off <- cmdscale(dist(data$hippo.off),eig = T,k=2)
to.plot.off <- data.frame(fit.off$points)
to.plot.off$Expression <- data$hippo.off
ggplot(to.plot.off,aes(y=to.plot.off$X1,
                       x=to.plot.off$X2,
                       color = to.plot.off$Expression,
                       label=rownames(data)))+
  geom_point()+
  ggrepel::geom_text_repel()+
  theme_classic()





scale.factor = 0.7
#### Look at Tuft cell WNT pathway activation ####
tuft <- subset(pneum.clean,idents = 'Tuft')
Idents(tuft) <- tuft$Timepoint

# Define a module
names <- c('Lgr5','Rnf43','Lrp5',
           'Lrp6','Fzd6','Ctnnb1',
           'Gsk3b','Ccnd1','Axin2',
           'Myc','Lef1','Tcf7',
           'Tcf7l1','Tcf7l2','Tle1',
           'Apc','Csnk1a1','Dvl1')
# grouncho homologue https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5592244/
# Define WNT ON state (Rspo1 binding)
# If it is being 'used more' when Rspo1 is bound and WNT is activated, then we would expect upregulation (1) on the transcriptomic level.
# If it is being 'used less' when Rspo1 is bound and WNT is activated, then we would expect downregulation (-1) on the transcriptional level.
wnt.on <- c(1,1,-1,
              -1,-1,-1,
              1,1,1,
              1,1,1,
              1,1,1,
              -1,-1,-1) 
# APC should be transcriptionally down when WNT is ON (proteomically up and increased in cytoplasm, less degradation)
# DVL should be transcriptomically down when WNT is ON
names(wnt.on) <- names
# Define OFF state
wnt.off <- -wnt.on

# Gather data for this module from a target population (Mesothelium, here)
temp.data <- tuft[names,]
meta.data <- temp.data@meta.data

# Compute ranges for each feature
data.temp <- as.data.frame(temp.data@assays$RNA@data)
ranges <- matrixStats::rowRanges(as.matrix(data.temp),na.rm = F) # do consider zeros when making space ranges?
#data.temp[data.temp==0] <- NA # ignore zeros this way

# Scale the ON/OFF states to the extrema of these ranges for each features
for (i in 1:length(wnt.on)){
  if(wnt.on[i]<0){wnt.on[i]<-ranges[names(wnt.on[i]),1]}else{wnt.on[i]<-ranges[names(wnt.on[i]),2]}
}
for (i in 1:length(wnt.off)){
  if(wnt.off[i]<0){wnt.off[i]<-ranges[names(wnt.off[i]),1]}else{wnt.off[i]<-ranges[names(wnt.off[i]),2]}
}
# Bind ON and OFF states
wnt.stat <- data.frame(wnt.on,wnt.off)
# Bind ON and OFF states to cellular states

data <- cbind(wnt.stat,data.temp)
#View(data)

# Transform into a (cell-cell) distance matrix
d <- dist(t(data),method = 'manhattan') # THIS NEEDS TO BE CAREFULLY REFINED 
#View(as.matrix(d))
#d[is.na(d)] <- 0 # workaround following https://github.com/tidymodels/corrr/issues/78
# Reduce dimension to 1 dimension (to visualize on a number line)
fit <- cmdscale(d,eig=T, k=1) # k is the number of dim
# see here for more thoughts: https://cran.r-project.org/web/packages/ordr/vignettes/cmds-variables.html
# Organize data for plotting
to.plot <- as.data.frame(fit$points)
to.plot$Sample <- NA
to.plot[rownames(meta.data),]$Sample <- as.character(meta.data$orig.ident)
to.plot['wnt.on',]$Sample <- 'Wnt-ON'
to.plot['wnt.off',]$Sample <- 'Wnt-OFF'
to.plot$Timepoint <- NA
to.plot[rownames(meta.data),]$Timepoint <- as.character(meta.data$Timepoint)
to.plot['wnt.on',]$Timepoint <- 'Wnt-ON'
to.plot['wnt.off',]$Timepoint <- 'Wnt-OFF'
# Organize metadata for plotting
to.plot$Timepoint <- factor(to.plot$Timepoint,levels = c('Day 0','Day 3','Day 7','Day 14'))

# Scale
to.plot$scale <- scale(to.plot$V1,center = T)[,1]
# Normalize
to.plot$normalized <- (to.plot$V1 - min(to.plot$V1)) / (max(to.plot$V1) - min(to.plot$V1))
# Isolate columns of interest
plot.temp <- to.plot[-c(1:2),]

# Fractions on both sides of the mean
frac.on.0 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 0',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 0',]),digits = 3)),'%',sep='')
frac.off.0 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 0',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 0',]),digits = 3)),'%',sep='')
frac.on.3 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 3',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 3',]),digits = 3)),'%',sep='')
frac.off.3 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 3',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 3',]),digits = 3)),'%',sep='')
frac.on.7 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 7',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 7',]),digits = 3)),'%',sep='')
frac.off.7 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 7',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 7',]),digits = 3)),'%',sep='')
frac.on.14 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 14',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 14',]),digits = 3)),'%',sep='')
frac.off.14 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 14',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 14',]),digits = 3)),'%',sep='')
dat_text <- data.frame(
  frac.on = c(frac.on.0, frac.on.3, frac.on.7,frac.on.14),
  frac.off = c(frac.off.0, frac.off.3, frac.off.7,frac.off.14),
  Timepoint   = as.factor(c('Day 0','Day 3','Day 7','Day 14'))
)

# Make base plot
plot.total <- ggplot(data=plot.temp, 
                     aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
  geom_density(mapping = aes(group = Sample),alpha=0.5,color=NA)+
  geom_jitter(aes(x=normalized,y=1))+
  scale_fill_manual(values=col.pal$Timepoint)+
  scale_color_manual(values=col.pal$Timepoint)+
  theme_classic()+
  xlab('WNT State')+
  ylab('Tuft Population Density')+
  #xlim(0,1)+ 
  # geom_vline(xintercept = 0, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 1, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 0.5, linetype="dotted", 
  #            color = "black", size=0.5)+
  geom_vline(xintercept = mean(plot.temp$normalized), linetype="dotted", 
             color = "black", size=0.5)+
  facet_wrap(~Timepoint,ncol=1)+ NoLegend()
# Add annotations
plot.total <- plot.total + geom_text(
  data    = dat_text,
  mapping = aes(x = 0.20, y = 5, label = frac.off),
  size=3.5,
  hjust   = -0.1,
  vjust   = -1
)+ geom_text(
  data    = dat_text,
  mapping = aes(x = 0.60, y = 5, label = frac.on),
  size=3.5,
  hjust   = -0.1,
  vjust   = -1
)

plot.total

# Make final plot
png(filename = 'Tuft WNT Pathway Status Densities VERTICAL.png',width = 5*scale.factor,height = 9*scale.factor,units = 'in',res = 300)
plot.total
dev.off()

# Stash
wnt.values <- plot.temp

# Experiment: a plot of the actual network ON vs OFF state?
fit.on <- cmdscale(dist(data$wnt.on),eig = T,k=2)
to.plot.on<- data.frame(fit.on$points)
to.plot.on$Expression <- data$wnt.on
ggplot(to.plot.on,aes(y=to.plot.on$X1,
                      x=to.plot.on$X2,
                      color = to.plot.on$Expression,
                      label=rownames(data)))+
  geom_point()+
  ggrepel::geom_text_repel()+
  theme_classic()

fit.off <- cmdscale(dist(data$wnt.off),eig = T,k=2)
to.plot.off <- data.frame(fit.off$points)
to.plot.off$Expression <- data$wnt.off
ggplot(to.plot.off,aes(y=to.plot.off$X1,
                       x=to.plot.off$X2,
                       color = to.plot.off$Expression,
                       label=rownames(data)))+
  geom_point()+
  ggrepel::geom_text_repel()+
  theme_classic()






#### Look at Tuft cell BMP-Smad pathway activation ####
tuft <- subset(pneum.clean,idents = 'Tuft')
Idents(tuft) <- tuft$Timepoint

# Define a module
names <- c('Smad1','Smad2','Smad3','Smad4',
           'Smad5','Smad6','Smad7','Bambi')

# Define BMP-SMAD ON state (Bmp4 binding)
# If it is being 'used more' when Bmp4 is bound and BMP is activated, then we would expect upregulation (1) on the transcriptomic level.
# If it is being 'used less' when Bmp4 is bound and BMP is activated, then we would expect downregulation (-1) on the transcriptional level.
bmp.smad.on <- c(1,1,1,1,
                 1,1,1,1) 
# doi:10.1006/dbio.2001.0388

names(bmp.smad.on) <- names
# Define OFF state
bmp.smad.off <- -bmp.smad.on

# Gather data for this module from a target population (Mesothelium, here)
temp.data <- tuft[names,]
meta.data <- temp.data@meta.data

# Compute ranges for each feature
data.temp <- as.data.frame(temp.data@assays$RNA@data)
ranges <- matrixStats::rowRanges(as.matrix(data.temp),na.rm = F) # do consider zeros when making space ranges?
#data.temp[data.temp==0] <- NA # ignore zeros this way

# Scale the ON/OFF states to the extrema of these ranges for each features
for (i in 1:length(bmp.smad.on)){
  if(bmp.smad.on[i]<0){bmp.smad.on[i]<-ranges[names(bmp.smad.on[i]),1]}else{bmp.smad.on[i]<-ranges[names(bmp.smad.on[i]),2]}
}
for (i in 1:length(bmp.smad.off)){
  if(bmp.smad.off[i]<0){bmp.smad.off[i]<-ranges[names(bmp.smad.off[i]),1]}else{bmp.smad.off[i]<-ranges[names(bmp.smad.off[i]),2]}
}
# Bind ON and OFF states
bmp.smad.stat <- data.frame(bmp.smad.on,bmp.smad.off)
# Bind ON and OFF states to cellular states

data <- cbind(bmp.smad.stat,data.temp)
#View(data)

# Transform into a (cell-cell) distance matrix
d <- dist(t(data),method = 'manhattan') # THIS NEEDS TO BE CAREFULLY REFINED 
#View(as.matrix(d))
#d[is.na(d)] <- 0 # workaround following https://github.com/tidymodels/corrr/issues/78
# Reduce dimension to 1 dimension (to visualize on a number line)
fit <- cmdscale(d,eig=T, k=1) # k is the number of dim
# see here for more thoughts: https://cran.r-project.org/web/packages/ordr/vignettes/cmds-variables.html
# Organize data for plotting
to.plot <- as.data.frame(fit$points)
to.plot$Sample <- NA
to.plot[rownames(meta.data),]$Sample <- as.character(meta.data$orig.ident)
to.plot['bmp.smad.on',]$Sample <- 'Bmp-ON'
to.plot['bmp.smad.off',]$Sample <- 'Bmp-OFF'
to.plot$Timepoint <- NA
to.plot[rownames(meta.data),]$Timepoint <- as.character(meta.data$Timepoint)
to.plot['bmp.smad.on',]$Timepoint <- 'Bmp-ON'
to.plot['bmp.smad.off',]$Timepoint <- 'Bmp-OFF'
# Organize metadata for plotting
to.plot$Timepoint <- factor(to.plot$Timepoint,levels = c('Day 0','Day 3','Day 7','Day 14'))

# Scale
to.plot$scale <- scale(to.plot$V1,center = T)[,1]
# Normalize
to.plot$normalized <- (to.plot$V1 - min(to.plot$V1)) / (max(to.plot$V1) - min(to.plot$V1))
# Isolate columns of interest
plot.temp <- to.plot[-c(1:2),]

# Fractions on both sides of the mean
frac.on.0 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 0',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 0',]),digits = 3)),'%',sep='')
frac.off.0 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 0',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 0',]),digits = 3)),'%',sep='')
frac.on.3 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 3',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 3',]),digits = 3)),'%',sep='')
frac.off.3 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 3',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 3',]),digits = 3)),'%',sep='')
frac.on.7 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 7',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 7',]),digits = 3)),'%',sep='')
frac.off.7 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 7',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 7',]),digits = 3)),'%',sep='')
frac.on.14 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 14',]$normalized>mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 14',]),digits = 3)),'%',sep='')
frac.off.14 <- paste(as.character(100*round(sum(plot.temp[plot.temp$Timepoint=='Day 14',]$normalized<mean(plot.temp$normalized))/nrow(plot.temp[plot.temp$Timepoint=='Day 14',]),digits = 3)),'%',sep='')
dat_text <- data.frame(
  frac.on = c(frac.on.0, frac.on.3, frac.on.7,frac.on.14),
  frac.off = c(frac.off.0, frac.off.3, frac.off.7,frac.off.14),
  Timepoint   = as.factor(c('Day 0','Day 3','Day 7','Day 14'))
)

# Make base plot
plot.total <- ggplot(data=plot.temp, 
                     aes(x=normalized, group=Timepoint, fill=Timepoint,color=Timepoint)) +
  geom_density(mapping = aes(group = Sample),alpha=0.5,color=NA)+
  geom_jitter(aes(x=normalized,y=1))+
  scale_fill_manual(values=col.pal$Timepoint)+
  scale_color_manual(values=col.pal$Timepoint)+
  theme_classic()+
  xlab('BMP-SMAD State')+
  ylab('Tuft Population Density')+
  #xlim(0,1)+ 
  # geom_vline(xintercept = 0, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 1, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 0.5, linetype="dotted", 
  #            color = "black", size=0.5)+
  geom_vline(xintercept = mean(plot.temp$normalized), linetype="dotted", 
             color = "black", size=0.5)+
  facet_wrap(~Timepoint,ncol=1)+ NoLegend()
plot.total
# Add annotations
plot.total <- plot.total + geom_text(
  data    = dat_text,
  mapping = aes(x = -0.03, y = 5, label = frac.off),
  size=3.5,
  hjust   = -0.1,
  vjust   = -1
)+ geom_text(
  data    = dat_text,
  mapping = aes(x = 0.4, y = 5, label = frac.on),
  size=3.5,
  hjust   = -0.1,
  vjust   = -1
)

plot.total

# Make final plot
png(filename = 'Tuft BMP-SMAD Pathway Status Densities VERTICAL.png',width = 5*scale.factor,height = 9*scale.factor,units = 'in',res = 300)
plot.total
dev.off()

# Stash
bmp.values <- plot.temp

# Experiment: a plot of the actual network ON vs OFF state?
fit.on <- cmdscale(dist(data$bmp.smad.on),eig = T,k=2)
to.plot.on<- data.frame(fit.on$points)
to.plot.on$Expression <- data$bmp.smad.on
ggplot(to.plot.on,aes(y=to.plot.on$X1,
                      x=to.plot.on$X2,
                      color = to.plot.on$Expression,
                      label=rownames(data)))+
  geom_point()+
  ggrepel::geom_text_repel()+
  theme_classic()

fit.off <- cmdscale(dist(data$bmp.smad.off),eig = T,k=2)
to.plot.off <- data.frame(fit.off$points)
to.plot.off$Expression <- data$bmp.smad.off
ggplot(to.plot.off,aes(y=to.plot.off$X1,
                       x=to.plot.off$X2,
                       color = to.plot.off$Expression,
                       label=rownames(data)))+
  geom_point()+
  ggrepel::geom_text_repel()+
  theme_classic() # All are zero, will not work




#### Experiment: Graph WNT vs. BMP in Tuft Cells ####
bmp <- bmp.values$normalized
wnt <- wnt.values$normalized
two.paths <- data.frame(bmp=bmp,wnt=wnt)
ggplot(two.paths,
       aes(x=bmp,y=wnt))+geom_point()
# this doesn't work because of data sparsity. good try though.

# try:
bmp.means <- bmp.values %>% group_by(Sample) %>% summarise(bmp = mean(normalized),n=n())
wnt.means <- wnt.values %>% group_by(Sample) %>% summarise(wnt = mean(normalized),n=n())

bmp.means$pathway <- 'BMP-SMAD'
wnt.means$pathway <- 'WNT'

two.paths <- data.frame(Sample = bmp.means$Sample,
                        BMP = bmp.means$bmp,
                        WNT = wnt.means$wnt)
ggplot(two.paths,
       aes(x=BMP,y=WNT,color=Sample))+geom_point()+theme_classic()

# what if we bootstrap via metacells?
# stash barcode information
bmp.values$barcode <- rownames(bmp.values)
wnt.values$barcode <- rownames(wnt.values)
# First lets get randomized samples of cells
n = 40 #(cells per sample to bootstrap)
p = 500 #(number of bootstraps)
cell.list <- list()
for(i in 1:p){
  set.seed <- i
  cell.list[[i]] <- bmp.values %>% group_by(Timepoint) %>% slice_sample(n = n) 
  cell.list[[i]] <- cell.list[[i]]$barcode
  }

# Then let's get the BMP and WNT values for each of these cell lists
bmp.list <- list()
for(i in 1:length(cell.list)){
  set.seed <- i
  bmp.list[[i]] <- bmp.values[cell.list[[i]],] 
}
wnt.list <- list()
for(i in 1:length(cell.list)){
  set.seed <- i
  wnt.list[[i]] <- wnt.values[cell.list[[i]],] 
}
# and then we can calculate the mean values for each of these data samplings
for(i in 1:length(bmp.list)){
bmp.list[[i]] <- bmp.list[[i]] %>% group_by(Sample) %>% summarise(BMP = mean(normalized),n=n())
bmp.list[[i]]$Bootstrap <- i
bmp.list[[i]] <- as.data.frame(bmp.list[[i]])
}
bmp.bootstrap.data <- bind_rows(bmp.list)
for(i in 1:length(wnt.list)){
  wnt.list[[i]] <- wnt.list[[i]] %>% group_by(Sample) %>% summarise(WNT = mean(normalized),n=n())
  wnt.list[[i]]$Bootstrap <- i
  wnt.list[[i]] <- as.data.frame(wnt.list[[i]])
}
wnt.bootstrap.data <- bind_rows(wnt.list)

wnt.bmp.joint.data <- data.frame(Sample = factor(wnt.bootstrap.data$Sample,levels = c('P0-f','P0-m','P3-f','P3-m','P7-f','P7-m','P14-f','P14-m')),
                                 WNT = wnt.bootstrap.data$WNT,
                                 BMP = bmp.bootstrap.data$BMP)
col.pal$Sample2 <- col.pal$Sample
names(col.pal$Sample2) <- c('P0-f','P0-m','P3-f','P3-m','P7-f','P7-m','P14-f','P14-m')

# Set up metadata
wnt.bmp.joint.data$Sex <- NA
wnt.bmp.joint.data[wnt.bmp.joint.data$Sample %in% c('P0-f','P3-f','P7-f','P14-f'),]$Sex <- 'Female'
wnt.bmp.joint.data[wnt.bmp.joint.data$Sample %in% c('P0-m','P3-m','P7-m','P14-m'),]$Sex <- 'Male'
wnt.bmp.joint.data$Timepoint <- NA
wnt.bmp.joint.data[wnt.bmp.joint.data$Sample %in% c('P0-f','P0-m'),]$Timepoint <- 'Day 0'
wnt.bmp.joint.data[wnt.bmp.joint.data$Sample %in% c('P3-m','P3-f'),]$Timepoint <- 'Day 3'
wnt.bmp.joint.data[wnt.bmp.joint.data$Sample %in% c('P7-f','P7-m'),]$Timepoint <- 'Day 7'
wnt.bmp.joint.data[wnt.bmp.joint.data$Sample %in% c('P14-m','P14-f'),]$Timepoint <- 'Day 14'
wnt.bmp.joint.data$Timepoint <- factor(wnt.bmp.joint.data$Timepoint,
                                levels = c('Day 0','Day 3','Day 7','Day 14'))
png(filename = 'BMP vs. WNT Scatterplot.png',width=8,height=2.4,res=300,units='in')
ggplot(wnt.bmp.joint.data,
       aes(x=BMP,y=WNT,color=Timepoint,group=Sex))+
  geom_point()+
  #scale_color_manual(values=col.pal$Sample2)+
  #scale_color_manual(values=c('red','blue'))+
  scale_color_manual(values=alpha(col.pal$Timepoint,alpha = 0.5))+
  theme_classic()+NoLegend()#+xlim(0.05,0.3)
  #ylim(0.05,0.28)+#geom_smooth(method = "lm", se = FALSE)
dev.off()





