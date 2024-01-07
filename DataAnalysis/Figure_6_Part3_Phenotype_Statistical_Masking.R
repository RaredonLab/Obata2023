# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell")
# Set Seed
set.seed(123)

# Load Packages
require(Seurat)
require(NICHES)
require(ggplot2)
require(metap)
require(ggrepel)

# Load colors and order
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneumonectomy.col.pal.Robj")

# Load previous, unfiltered, statistical significance tables
load("pneum.clean.data.split.seurat.markers.unfiltered.Robj")

#### Load niches data to create vectortype indices ####
# Load  phenotype data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneum.clean.annotated.2023-06-14.Robj")

# Inspect data
table(pneum.clean$Sample)
table(pneum.clean$CellType)

# Split into two objects, one for each sex / timecourse
pneum.clean.data <- list()
pneum.clean.data[['male']] <- subset(pneum.clean,subset = Sample %in% c('P0_m','P3_m','P7_m','P14_m'))
pneum.clean.data[['female']]  <- subset(pneum.clean,subset = Sample %in% c('P0_f','P3_f','P7_f','P14_f'))

# Split each sex into 34 CellTypes
# Capture metadata
metadata.male <- pneum.clean.data[['male']]@meta.data
metadata.female <- pneum.clean.data[['female']]@meta.data
# Split metadata by CellType
metadata.male.split <- split(x = metadata.male, f = metadata.male$CellType)
metadata.female.split <- split(x = metadata.female, f = metadata.female$CellType)

# Split phenotype data by CellType
# Identify mechanisms and CellTypes to evaluate (same for both male and female)
features <- rownames(pneum.clean)
sum(!(names(metadata.male.split) == names(metadata.female.split))) # checking that it equals 0
CellTypes <- names(metadata.male.split)
CellTypes.mod <- gsub('—','.',CellTypes) # Making string separator compatible with downstream work

#### ####
# # Add useful metadata
# for (i in 1:length(CellTypes)){
#   print(i)
#   # Male
#   if(nrow(pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]])>0){
#   pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$SendingType <- unique(metadata.male.split[[CellTypes[i]]]$SendingType)
#   pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$ReceivingType <- unique(metadata.male.split[[CellTypes[i]]]$ReceivingType)
#   if(unique(metadata.male.split[[CellTypes[i]]]$SendingType)!='Cell_cycle'){
#   pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$CellClass.Sending <- unique(metadata.male.split[[CellTypes[i]]]$CellClass.Sending)
#   }else{
#     pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$CellClass.Sending <- NA
#   }
#   if(unique(metadata.male.split[[CellTypes[i]]]$ReceivingType)!='Cell_cycle'){
#   pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$CellClass.Receiving <- unique(metadata.male.split[[CellTypes[i]]]$CellClass.Receiving)
#   }else{
#     pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$CellClass.Receiving <- NA
#   }
#   }
#   # Female
#   if(nrow(pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]])>0){
#   pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$SendingType <- unique(metadata.female.split[[CellTypes[i]]]$SendingType)
#   pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$ReceivingType <- unique(metadata.female.split[[CellTypes[i]]]$ReceivingType)
#   if(unique(metadata.female.split[[CellTypes[i]]]$SendingType)!='Cell_cycle'){
#     pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$CellClass.Sending <- unique(metadata.female.split[[CellTypes[i]]]$CellClass.Sending)
#   }else{
#     pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$CellClass.Sending <- NA
#   }
#   if(unique(metadata.female.split[[CellTypes[i]]]$ReceivingType)!='Cell_cycle'){
#     pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$CellClass.Receiving <- unique(metadata.female.split[[CellTypes[i]]]$CellClass.Receiving)
#   }else{
#     pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$CellClass.Receiving <- NA
#   }
#   }
# }

# Make into one big data frame for each sex
# Male
#male.sig.findings <- do.call(rbind,pneum.clean.data.split.seurat.markers.significant[['male']])
male.findings <- do.call(rbind,pneum.clean.data.split.seurat.markers[['male']])
# Female
#female.sig.findings <- do.call(rbind,pneum.clean.data.split.seurat.markers.significant[['female']])
female.findings <- do.call(rbind,pneum.clean.data.split.seurat.markers[['female']])


# Limit to just those things meeting the earlier, very minimal, assessment thresholds in both sexes
length(male.findings$unique)
length(female.findings$unique)
unique.intersecting <- intersect(male.findings$unique,female.findings$unique)
length(unique.intersecting)
male.findings.unique <- male.findings[male.findings$unique %in% unique.intersecting,]
female.findings.unique <- female.findings[female.findings$unique %in% unique.intersecting,]
nrow(male.findings.unique)
nrow(female.findings.unique)

# Order the same
male.findings.unique <- male.findings.unique[order(male.findings.unique$unique),]
female.findings.unique <- female.findings.unique[order(female.findings.unique$unique),]
which(is.na(male.findings.unique$unique)) # =0
which(is.na(female.findings.unique$unique)) # =0
rownames(male.findings.unique[which(is.na(male.findings.unique$unique)),])==rownames(female.findings.unique[which(is.na(female.findings.unique$unique)),]) # none
sum(!(male.findings.unique$unique==female.findings.unique$unique),na.rm = T) # = 0 --> ordered the same

# Combine into a meta findings list (no real thresholding yet)
# Calculate combined p-values
combined.p.vals <- metapod::combineParallelPValues(list(male.findings.unique$p_val_adj,female.findings.unique$p_val_adj),
                                method = 'fisher',
                                log.p=FALSE)
combined.p.vals$p.value
combined.findings <- data.frame(Day=male.findings.unique$cluster,
                                CellType = male.findings.unique$CellType,
                                #ReceivingType = male.findings.unique$ReceivingType,
                                #CellClass.Sending = male.findings.unique$CellClass.Sending,
                                #CellClass.Receiving = male.findings.unique$CellClass.Receiving,
                                Gene = male.findings.unique$gene,
                                male.avg_log2FC = male.findings.unique$avg_log2FC,
                                female.avg_log2FC = female.findings.unique$avg_log2FC,
                                male.p_val_adj = male.findings.unique$p_val_adj,
                                female.p_val_adj = female.findings.unique$p_val_adj,
                                male.pct.1 = male.findings.unique$pct.1,
                                female.pct.1 = female.findings.unique$pct.1,
                                male.pct.2 = male.findings.unique$pct.2,
                                female.pct.2 = female.findings.unique$pct.2,
                                male.logFC.sign = male.findings.unique$avg_log2FC>0,
                                female.logFC.sign = female.findings.unique$avg_log2FC>0,
                                mean.pct.1=rowMeans(cbind(male.findings.unique$pct.1,female.findings.unique$pct.1)),
                                mean.pct.2=rowMeans(cbind(male.findings.unique$pct.2,female.findings.unique$pct.2)),
                                mean.log2FC=rowMeans(cbind(male.findings.unique$avg_log2FC,female.findings.unique$avg_log2FC)),
                                comb.p.val.adj=combined.p.vals$p.value,
                                #VectorType = male.findings.unique$VectorType,
                                unique = male.findings.unique$unique)
combined.findings$delta <- 'NA'
combined.findings$delta[combined.findings$mean.log2FC>0] <- 'Positive'
combined.findings$delta[combined.findings$mean.log2FC<0] <- 'Negative'

# Establish consensus between sexes
combined.findings$consensus <- FALSE
combined.findings$consensus[combined.findings$male.logFC.sign==combined.findings$female.logFC.sign] <- TRUE
View(combined.findings) #Beautiful! # 1,820,774

# Save point
save(combined.findings,file = 'combined.findings.phenotype.2023-06-28.Robj')

# Apply statistical thresholding
# P-value threshold
p.adj.thresh <- 0.05
#combined.findings.thresh <- combined.findings[combined.findings$male.p_val_adj<p.adj.thresh,] # 33,482
#combined.findings.thresh <- combined.findings.thresh[combined.findings.thresh$female.p_val_adj<p.adj.thresh,] # 7483
combined.findings.thresh <- combined.findings[combined.findings$comb.p.val.adj<p.adj.thresh,] # 89,518

# Require consensus signal change direction across sexes
combined.findings.thresh <- combined.findings.thresh[combined.findings.thresh$consensus==TRUE,] # 48,631
View(combined.findings.thresh)

# Save point
save(combined.findings.thresh,file = 'combined.findings.phenotype.thresh.2023-06-28.Robj')

# Investigate trends in LIGANDS
combined.findings.thresh.ligands <- subset(combined.findings.thresh[combined.findings.thresh$Gene %in% NICHES::ncomms8866_rat$Ligand.ApprovedSymbol,])
exclude <- grep('Scgb|Sftp|Hsp',combined.findings.thresh.ligands$Gene)
combined.findings.thresh.ligands <- combined.findings.thresh.ligands[-exclude,]
save(combined.findings.thresh.ligands,file = 'combined.findings.thresh.ligands.2023-06-28.Robj')

day0 <- subset(combined.findings.thresh.ligands,Day=='Day 0')
day3 <- subset(combined.findings.thresh.ligands,Day=='Day 3')
day7 <- subset(combined.findings.thresh.ligands,Day=='Day 7')
day14 <- subset(combined.findings.thresh.ligands,Day=='Day 14')
day0.sum <- table(day0$CellType,day0$delta)
day3.sum <- table(day3$CellType,day3$delta)
day7.sum <- table(day7$CellType,day7$delta)
day14.sum <- table(day14$CellType,day14$delta)
day0.sum <- cbind(day0.sum,rowSums(day0.sum))
day3.sum <- cbind(day3.sum,rowSums(day3.sum))
day7.sum <- cbind(day7.sum,rowSums(day7.sum))
day14.sum <- cbind(day14.sum,rowSums(day14.sum))
colnames(day0.sum)[3] <- 'Total'
colnames(day3.sum)[3] <- 'Total'
colnames(day7.sum)[3] <- 'Total'
colnames(day14.sum)[3] <- 'Total'
day0.sum <- data.frame(day0.sum)
day3.sum <- data.frame(day3.sum)
day7.sum <- data.frame(day7.sum)
day14.sum <- data.frame(day14.sum)
day0.sum$Day <- 'Day 0'
day3.sum$Day <- 'Day 3'
day7.sum$Day <- 'Day 7'
day14.sum$Day <- 'Day 14'
day0.sum$CellType <- rownames(day0.sum)
day3.sum$CellType <- rownames(day3.sum)
day7.sum$CellType <- rownames(day7.sum)
day14.sum$CellType <- rownames(day14.sum)
day0.sum$Negative <- -day0.sum$Negative
day3.sum$Negative <- -day3.sum$Negative
day7.sum$Negative <- -day7.sum$Negative
day14.sum$Negative <- -day14.sum$Negative
day0.sum <- reshape2::melt(day0.sum)
day3.sum <- reshape2::melt(day3.sum)
day7.sum <- reshape2::melt(day7.sum)
day14.sum <- reshape2::melt(day14.sum)
day0.sum$CellType <- factor(day0.sum$CellType,levels = rev(names(col.pal$CellType)))
day3.sum$CellType <- factor(day3.sum$CellType,levels = rev(names(col.pal$CellType)))
day7.sum$CellType <- factor(day7.sum$CellType,levels = rev(names(col.pal$CellType)))
day14.sum$CellType <- factor(day14.sum$CellType,levels = rev(names(col.pal$CellType)))

# Define plotting limits
total.data <- rbind(day0.sum,day3.sum,day7.sum,day14.sum)
x.max <- max(total.data[total.data$variable=='Positive','value'])
x.min <- min(total.data[total.data$variable=='Negative','value'])
if(x.max>-x.min){lim <- x.max}else{lim <- -x.min}
# Plot as a stacked bar chart
p1 <- ggplot(data = day0.sum[day0.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
               color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Ligands')+
  labs(fill = 'Delta')+
  ggtitle('Day 0 Ligands')

p2 <- ggplot(data = day3.sum[day3.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Ligands')+
  labs(fill = 'Delta')+
  ggtitle('Day 3 Ligands')

p3 <- ggplot(data = day7.sum[day7.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Ligands')+
  labs(fill = 'Delta')+
  ggtitle('Day 7 Ligands')

p4 <- ggplot(data = day14.sum[day14.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Ligands')+
  labs(fill = 'Delta')+
  ggtitle('Day 14 Ligands')

# Investigate trends in RECEPTORS
combined.findings.thresh.receptors <- subset(combined.findings.thresh[combined.findings.thresh$Gene %in% NICHES::ncomms8866_rat$Receptor.ApprovedSymbol,])
save(combined.findings.thresh.receptors,file = 'combined.findings.thresh.receptors.2023-06-28.Robj')

day0 <- subset(combined.findings.thresh.receptors,Day=='Day 0')
day3 <- subset(combined.findings.thresh.receptors,Day=='Day 3')
day7 <- subset(combined.findings.thresh.receptors,Day=='Day 7')
day14 <- subset(combined.findings.thresh.receptors,Day=='Day 14')
day0.sum <- table(day0$CellType,day0$delta)
day3.sum <- table(day3$CellType,day3$delta)
day7.sum <- table(day7$CellType,day7$delta)
day14.sum <- table(day14$CellType,day14$delta)
day0.sum <- cbind(day0.sum,rowSums(day0.sum))
day3.sum <- cbind(day3.sum,rowSums(day3.sum))
day7.sum <- cbind(day7.sum,rowSums(day7.sum))
day14.sum <- cbind(day14.sum,rowSums(day14.sum))
colnames(day0.sum)[3] <- 'Total'
colnames(day3.sum)[3] <- 'Total'
colnames(day7.sum)[3] <- 'Total'
colnames(day14.sum)[3] <- 'Total'
day0.sum <- data.frame(day0.sum)
day3.sum <- data.frame(day3.sum)
day7.sum <- data.frame(day7.sum)
day14.sum <- data.frame(day14.sum)
day0.sum$Day <- 'Day 0'
day3.sum$Day <- 'Day 3'
day7.sum$Day <- 'Day 7'
day14.sum$Day <- 'Day 14'
day0.sum$CellType <- rownames(day0.sum)
day3.sum$CellType <- rownames(day3.sum)
day7.sum$CellType <- rownames(day7.sum)
day14.sum$CellType <- rownames(day14.sum)
day0.sum$Negative <- -day0.sum$Negative
day3.sum$Negative <- -day3.sum$Negative
day7.sum$Negative <- -day7.sum$Negative
day14.sum$Negative <- -day14.sum$Negative
day0.sum <- reshape2::melt(day0.sum)
day3.sum <- reshape2::melt(day3.sum)
day7.sum <- reshape2::melt(day7.sum)
day14.sum <- reshape2::melt(day14.sum)
day0.sum$CellType <- factor(day0.sum$CellType,levels = rev(names(col.pal$CellType)))
day3.sum$CellType <- factor(day3.sum$CellType,levels = rev(names(col.pal$CellType)))
day7.sum$CellType <- factor(day7.sum$CellType,levels = rev(names(col.pal$CellType)))
day14.sum$CellType <- factor(day14.sum$CellType,levels = rev(names(col.pal$CellType)))

# Define plotting limits
total.data <- rbind(day0.sum,day3.sum,day7.sum,day14.sum)
x.max <- max(total.data[total.data$variable=='Positive','value'])
x.min <- min(total.data[total.data$variable=='Negative','value'])
if(x.max>-x.min){lim <- x.max}else{lim <- -x.min}
# Plot as a stacked bar chart
p5 <- ggplot(data = day0.sum[day0.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Receptors')+
  labs(fill = 'Delta')+
  ggtitle('Day 0 Receptors')

p6 <- ggplot(data = day3.sum[day3.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Receptors')+
  labs(fill = 'Delta')+
  ggtitle('Day 3 Receptors')

p7 <- ggplot(data = day7.sum[day7.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Receptors')+
  labs(fill = 'Delta')+
  ggtitle('Day 7 Receptors')

p8 <- ggplot(data = day14.sum[day14.sum$variable!='Total',],
       aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Receptors')+
  labs(fill = 'Delta')+
  ggtitle('Day 14 Receptors')

plot_grid(p1,p2,p3,p4,
          p5,p6,p7,p8,nrow=2,align = 'hv')


# Investigate trends in TOTAL GENES
combined.findings.thresh.to.plot <- combined.findings.thresh
exclude <- grep('Scgb|Sftp|Hsp|AABR|7SK.|LOC|Rp|AC',ignore.case = FALSE,combined.findings.thresh.to.plot$Gene)
combined.findings.thresh.to.plot <- combined.findings.thresh.to.plot[-exclude,]
save(combined.findings.thresh.to.plot,file = 'combined.findings.thresh.to.plot.2023-06-28.Robj')

day0 <- subset(combined.findings.thresh.to.plot,Day=='Day 0')
day3 <- subset(combined.findings.thresh.to.plot,Day=='Day 3')
day7 <- subset(combined.findings.thresh.to.plot,Day=='Day 7')
day14 <- subset(combined.findings.thresh.to.plot,Day=='Day 14')
day0.sum <- table(day0$CellType,day0$delta)
day3.sum <- table(day3$CellType,day3$delta)
day7.sum <- table(day7$CellType,day7$delta)
day14.sum <- table(day14$CellType,day14$delta)
day0.sum <- cbind(day0.sum,rowSums(day0.sum))
day3.sum <- cbind(day3.sum,rowSums(day3.sum))
day7.sum <- cbind(day7.sum,rowSums(day7.sum))
day14.sum <- cbind(day14.sum,rowSums(day14.sum))
colnames(day0.sum)[3] <- 'Total'
colnames(day3.sum)[3] <- 'Total'
colnames(day7.sum)[3] <- 'Total'
colnames(day14.sum)[3] <- 'Total'
day0.sum <- data.frame(day0.sum)
day3.sum <- data.frame(day3.sum)
day7.sum <- data.frame(day7.sum)
day14.sum <- data.frame(day14.sum)
day0.sum$Day <- 'Day 0'
day3.sum$Day <- 'Day 3'
day7.sum$Day <- 'Day 7'
day14.sum$Day <- 'Day 14'
day0.sum$CellType <- rownames(day0.sum)
day3.sum$CellType <- rownames(day3.sum)
day7.sum$CellType <- rownames(day7.sum)
day14.sum$CellType <- rownames(day14.sum)
day0.sum$Negative <- -day0.sum$Negative
day3.sum$Negative <- -day3.sum$Negative
day7.sum$Negative <- -day7.sum$Negative
day14.sum$Negative <- -day14.sum$Negative
day0.sum <- reshape2::melt(day0.sum)
day3.sum <- reshape2::melt(day3.sum)
day7.sum <- reshape2::melt(day7.sum)
day14.sum <- reshape2::melt(day14.sum)
day0.sum$CellType <- factor(day0.sum$CellType,levels = rev(names(col.pal$CellType)))
day3.sum$CellType <- factor(day3.sum$CellType,levels = rev(names(col.pal$CellType)))
day7.sum$CellType <- factor(day7.sum$CellType,levels = rev(names(col.pal$CellType)))
day14.sum$CellType <- factor(day14.sum$CellType,levels = rev(names(col.pal$CellType)))

# Define plotting limits
total.data <- rbind(day0.sum,day3.sum,day7.sum,day14.sum)
x.max <- max(total.data[total.data$variable=='Positive','value'])
x.min <- min(total.data[total.data$variable=='Negative','value'])
if(x.max>-x.min){lim <- x.max}else{lim <- -x.min}

# Plot as a stacked bar chart
p9 <- ggplot(data = day0.sum[day0.sum$variable!='Total',],
             aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Genes')+
  labs(fill = 'Delta')+
  ggtitle('Day 0 Genes')

p10 <- ggplot(data = day3.sum[day3.sum$variable!='Total',],
             aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Genes')+
  labs(fill = 'Delta')+
  ggtitle('Day 3 Genes')

p11 <- ggplot(data = day7.sum[day7.sum$variable!='Total',],
             aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Genes')+
  labs(fill = 'Delta')+
  ggtitle('Day 7 Genes')

p12 <- ggplot(data = day14.sum[day14.sum$variable!='Total',],
             aes(x=value,y=CellType,fill=variable))+
  geom_bar(position="stack", stat="identity")+
  xlim(-lim,lim)+
  scale_y_discrete(drop=FALSE)+
  scale_fill_manual(values=c('red3','blue3'))+
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", linewidth=0.25)+
  theme_classic()+
  xlab('Number of Genes')+
  labs(fill = 'Delta')+
  ggtitle('Day 14 Genes')

# Make a big summary plot for supplement
png(filename='differential.expression.phenotype.stacked.bar.png',width=22,height=11,units='in',res=300)
plot_grid(p1,p2,p3,p4,
          p5,p6,p7,p8,
          p9,p10,p11,p12,nrow=3,align = 'hv')
dev.off()
