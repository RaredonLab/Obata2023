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

# Load previous, unfiltered, phenotype significance findings and stash to prevent overwriting
load('combined.findings.phenotype.2023-06-28.Robj')
combined.findings.phenotype <- combined.findings # unthresholded phenotype significance testing
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/combined.findings.phenotype.thresh.2023-06-28.Robj")
combined.findings.phenotype.thresh <- combined.findings.thresh # limited phenotype significance testing

# Load previous, unfiltered, statistical significance tables
load("niches.imputed.data.split.seurat.markers.by.timepoint.unfiltered.Robj")

#### Load niches data to create vectortype indices ####
# Load NICHES data
load('cell.to.cell.imputed.by.timepoint.Robj')

# Inspect data
table(cell.to.cell.imputed.by.timepoint$Sample.Sending)
table(cell.to.cell.imputed.by.timepoint$VectorType)

# Split into two objects, one for each sex / timecourse
niches.imputed.data <- list()
niches.imputed.data[['male']] <- subset(cell.to.cell.imputed.by.timepoint,subset = Sample.Sending %in% c('P0_m','P3_m','P7_m','P14_m'))
niches.imputed.data[['female']]  <- subset(cell.to.cell.imputed.by.timepoint,subset = Sample.Sending %in% c('P0_f','P3_f','P7_f','P14_f'))

# Split each sex into 1156 vectortypes
# Capture metadata
metadata.male <- niches.imputed.data[['male']]@meta.data
metadata.female <- niches.imputed.data[['female']]@meta.data
# Split metadata by vectortype
metadata.male.split <- split(x = metadata.male, f = metadata.male$VectorType)
metadata.female.split <- split(x = metadata.female, f = metadata.female$VectorType)

# Split NICHES data by Vectortype
# Identify mechanisms and vectortypes to evaluate (same for both male and female)
mechs <- rownames(cell.to.cell.imputed.by.timepoint)
sum(!(names(metadata.male.split) == names(metadata.female.split))) # checking that it equals 0
VectorTypes <- names(metadata.male.split)
VectorTypes.mod <- gsub('—','.',VectorTypes) # Making string separator compatible with downstream work

#### ####
# Add useful metadata
for (i in 1:length(VectorTypes)){
  print(i)
  # Male
  if(nrow(niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]])>0){
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$SendingType <- unique(metadata.male.split[[VectorTypes[i]]]$SendingType)
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$ReceivingType <- unique(metadata.male.split[[VectorTypes[i]]]$ReceivingType)
  if(unique(metadata.male.split[[VectorTypes[i]]]$SendingType)!='Cell_cycle'){
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$CellClass.Sending <- unique(metadata.male.split[[VectorTypes[i]]]$CellClass.Sending)
  }else{
    niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$CellClass.Sending <- NA
  }
  if(unique(metadata.male.split[[VectorTypes[i]]]$ReceivingType)!='Cell_cycle'){
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$CellClass.Receiving <- unique(metadata.male.split[[VectorTypes[i]]]$CellClass.Receiving)
  }else{
    niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$CellClass.Receiving <- NA
  }
  }
  # Female
  if(nrow(niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]])>0){
  niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$SendingType <- unique(metadata.female.split[[VectorTypes[i]]]$SendingType)
  niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$ReceivingType <- unique(metadata.female.split[[VectorTypes[i]]]$ReceivingType)
  if(unique(metadata.female.split[[VectorTypes[i]]]$SendingType)!='Cell_cycle'){
    niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$CellClass.Sending <- unique(metadata.female.split[[VectorTypes[i]]]$CellClass.Sending)
  }else{
    niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$CellClass.Sending <- NA
  }
  if(unique(metadata.female.split[[VectorTypes[i]]]$ReceivingType)!='Cell_cycle'){
    niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$CellClass.Receiving <- unique(metadata.female.split[[VectorTypes[i]]]$CellClass.Receiving)
  }else{
    niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$CellClass.Receiving <- NA
  }
  }
}

# Make into one big data frame for each sex
# Male
#male.sig.findings <- do.call(rbind,niches.imputed.data.split.seurat.markers.significant[['male']])
male.findings <- do.call(rbind,niches.imputed.data.split.seurat.markers[['male']])
# Female
#female.sig.findings <- do.call(rbind,niches.imputed.data.split.seurat.markers.significant[['female']])
female.findings <- do.call(rbind,niches.imputed.data.split.seurat.markers[['female']])


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
                                SendingType = male.findings.unique$SendingType,
                                ReceivingType = male.findings.unique$ReceivingType,
                                CellClass.Sending = male.findings.unique$CellClass.Sending,
                                CellClass.Receiving = male.findings.unique$CellClass.Receiving,
                                Mechanism = male.findings.unique$gene,
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
                                VectorType = male.findings.unique$VectorType,
                                unique = male.findings.unique$unique)
combined.findings$delta <- 'NA'
combined.findings$delta[combined.findings$mean.log2FC>0] <- 'Positive'
combined.findings$delta[combined.findings$mean.log2FC<0] <- 'Negative'

# Establish consensus between sexes
combined.findings$consensus <- FALSE
combined.findings$consensus[combined.findings$male.logFC.sign==combined.findings$female.logFC.sign] <- TRUE
View(combined.findings) #Beautiful! # 4,100,011

# Save point
save(combined.findings,file = 'combined.findings.imputed.by.timepoint.2023-06-25.Robj') #"combined.findings"
load('combined.findings.imputed.by.timepoint.2023-06-25.Robj')


#### Apply statistical thresholding ####

# P-value threshold
p.adj.thresh <- 0.0001
combined.findings.thresh <- combined.findings[combined.findings$male.p_val_adj<p.adj.thresh,] # 324,559
combined.findings.thresh <- combined.findings.thresh[combined.findings.thresh$female.p_val_adj<p.adj.thresh,] # 136,046
#combined.findings.thresh <- combined.findings[combined.findings$comb.p.val.adj<p.adj.thresh,] # 695,428

# Require consensus signal change direction across sexes
combined.findings.thresh <- combined.findings.thresh[combined.findings.thresh$consensus==TRUE,] # 113,726

# Require that trends agree with the phenotype findings, on either the sending or receiving side
# Add Ligand and receptor data
combined.findings.thresh$Ligand <- stringr::str_split_fixed(combined.findings.thresh$Mechanism,'—',n=2)[,1]
combined.findings.thresh$Receptor <- stringr::str_split_fixed(combined.findings.thresh$Mechanism,'—',n=2)[,2]
# Build consensus handles for connectomic data
combined.findings.thresh$Day.SendingType.Ligand.global.delta <- paste(combined.findings.thresh$Day,
                                                      combined.findings.thresh$SendingType,
                                                      combined.findings.thresh$Ligand,
                                                      combined.findings.thresh$delta)
combined.findings.thresh$Day.ReceivingType.Receptor.global.delta <- paste(combined.findings.thresh$Day,
                                                         combined.findings.thresh$ReceivingType,
                                                         combined.findings.thresh$Receptor,
                                                         combined.findings.thresh$delta)
combined.findings.thresh$Day.SendingType.Ligand <- paste(combined.findings.thresh$Day,
                                                                      combined.findings.thresh$SendingType,
                                                                      combined.findings.thresh$Ligand)
combined.findings.thresh$Day.ReceivingType.Receptor <- paste(combined.findings.thresh$Day,
                                                                          combined.findings.thresh$ReceivingType,
                                                                          combined.findings.thresh$Receptor)
# Build consensus handles for phenotype data (using non-thresholded phenotype data)
combined.findings.phenotype$Day.CellType.Gene.delta <- paste(combined.findings.phenotype$Day,
                                                             combined.findings.phenotype$CellType,
                                                             combined.findings.phenotype$Gene,
                                                             combined.findings.phenotype$delta)
combined.findings.phenotype$Day.CellType.Gene <- paste(combined.findings.phenotype$Day,
                                                             combined.findings.phenotype$CellType,
                                                             combined.findings.phenotype$Gene)
# Limit findings to things with consistent trends between male and female in the RNA data
combined.findings.phenotype.limited <- combined.findings.phenotype[combined.findings.phenotype$consensus==TRUE,]

# Filter connectomic findings based on pass/nopass columns
# Pct1 must be above threshold
combined.findings.phenotype.limited <- combined.findings.phenotype.limited[combined.findings.phenotype.limited$mean.pct.1>0.04,]

combined.findings.niches.thresh <- combined.findings.thresh[combined.findings.thresh$Day.SendingType.Ligand %in%
                                                              combined.findings.phenotype.limited$Day.CellType.Gene,] #60853
combined.findings.niches.thresh <- combined.findings.niches.thresh[combined.findings.niches.thresh$Day.ReceivingType.Receptor %in%
                                                              combined.findings.phenotype.limited$Day.CellType.Gene,] #37959
# Delta must agree
combined.findings.niches.thresh <- combined.findings.niches.thresh[combined.findings.niches.thresh$Day.SendingType.Ligand.global.delta %in% combined.findings.phenotype.limited$Day.CellType.Gene.delta |
                                                                     combined.findings.niches.thresh$Day.ReceivingType.Receptor.global.delta %in% combined.findings.phenotype.limited$Day.CellType.Gene.delta,] # 37,555
View(combined.findings.niches.thresh)

# Save point
save(combined.findings.niches.thresh,file = 'combined.findings.niches.thresh.2023-06-28.Robj')
