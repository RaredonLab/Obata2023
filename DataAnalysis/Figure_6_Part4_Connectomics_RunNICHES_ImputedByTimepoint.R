# Figure 6 Part 2 Pneumonectomy Connectomics

# Set WD
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell")
# Set Seed
set.seed(123)

# Load Packages
require(Seurat)
require(NICHES)
require(ggplot2)

# Load and Inspect Data
load('pneum.clean.annotated.2023-06-14.Robj')
table(Idents(pneum.clean))
table(pneum.clean$Sample)

####  Split by Timepoint for Imputation #### 
split <- SplitObject(pneum.clean,split.by = 'Timepoint')

####  Impute selected genes ####
# Pick which genes to impute
num.cells.per.feature <- 50
genes.to.impute <- rownames(pneum.clean)[rowSums(pneum.clean@assays$RNA@counts>0)>num.cells.per.feature]
# impute using ALRA
options(warn = 1)
for(i in 1:length(split)){
  print(i)
  split[[i]] <- SeuratWrappers::RunALRA(split[[i]],
                                        genes.use = genes.to.impute)
  gc()
}
# Save imputed data by timepoint
gc()
pneum.clean.imputed.by.timepoint <- merge(split[[1]],split[2:length(split)])
save(pneum.clean.imputed.by.timepoint,file='pneum.clean.imputed.by.timepoint.Robj')

####  Run NICHES [IMPUTED] #### 
# Need to re-split by sample first
gc()
split <- SplitObject(pneum.clean.imputed.by.timepoint,split.by='Sample')
names(split)

# CTC 
scc.list <- list()
for(i in 1:length(split)){
  print(i)
  scc.list[[i]] <- RunNICHES(split[[i]],
                             LR.database="fantom5",
                             species="rat",
                             assay="alra",
                             cell_types = "CellType",
                             meta.data.to.map = names(split[[i]]@meta.data),
                             SystemToCell = F,
                             CellToCell = T,
                             CellToSystem=F) # Need to standarize the number of cells in the System measurements
}
names(scc.list) <- names(split)

# Merge outputs
temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToCell # Isolate CellToCell Signaling
  gc()
}
cell.to.cell.imputed.by.timepoint <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.cell.imputed.by.timepoint

# Format metadata for plotting
cell.to.cell.imputed.by.timepoint$Timepoint.Sending <- factor(cell.to.cell.imputed.by.timepoint$Timepoint.Sending,
                                                              levels = c('Day 0','Day 3','Day 7','Day 14'))
cell.to.cell.imputed.by.timepoint$Sample.Sending <- factor(cell.to.cell.imputed.by.timepoint$Sample.Sending,
                                                           levels = c('P0_f','P0_m','P3_f','P3_m',
                                                                      'P7_f','P7_m','P14_f','P14_m'))

save(cell.to.cell.imputed.by.timepoint,file='cell.to.cell.imputed.by.timepoint.Robj')

# STC and CTS
# Standardizing cell number for each sample
split <- SplitObject(pneum.clean.imputed.by.timepoint,split.by = 'Sample')
max.cells.per.sample <- min(table(pneum.clean.imputed.by.timepoint$Sample)) #6122

scc.list <- list()
for(i in 1:length(split)){
  print(i)
  Idents(split[[i]]) <- split[[i]]$orig.ident
  split[[i]] <- subset(split[[i]],cells = WhichCells(split[[i]],downsample = max.cells.per.sample))
  scc.list[[i]] <- RunNICHES(split[[i]],
                             LR.database="fantom5",
                             species="rat",
                             assay="alra",
                             cell_types = "CellType",
                             meta.data.to.map = names(split[[i]]@meta.data),
                             SystemToCell = T,
                             CellToCell = F,
                             CellToSystem=T) 
}
names(scc.list) <- names(split)

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$SystemToCell # Isolate SystemToCell Signaling
}
system.to.cell <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
system.to.cell
save(system.to.cell,file='system.to.cell.imputed.by.timepoint.Robj')

temp.list <- list()
for(i in 1:length(scc.list)){
  temp.list[[i]] <- scc.list[[i]]$CellToSystem # Isolate CellToSystem Signaling
}
cell.to.system <- merge(temp.list[[1]],temp.list[2:length(temp.list)])
cell.to.system
save(cell.to.system,file='cell.to.system.imputed.by.timepoint.Robj')

# # Clean NICHES Data
# VlnPlot(cell.to.cell,
#         features = 'nFeature_CellToCell',
#         group.by = 'Sample.Joint',
#         pt.size=0,log = T,raster = T)+NoLegend()
# pneum.ctc.sub <- subset(cell.to.cell,nFeature_CellToCell>100) # Choose this limit based on the above violin plots or similar
# pneum.ctc.sub
# save(pneum.ctc.sub,file = 'pneum.cell.to.cell.sub.Robj')
# gc()
# VlnPlot(system.to.cell,
#         features = 'nFeature_SystemToCell',
#         group.by = 'Sample',
#         pt.size=0,log = T)+NoLegend()
# pneum.stc.sub <- subset(system.to.cell,nFeature_SystemToCell>100) # Choose this limit based on the above violin plots or similar
# pneum.stc.sub
# save(pneum.stc.sub,file = 'pneum.system.to.cell.sub.Robj')
# gc()
# VlnPlot(cell.to.system,
#         features = 'nFeature_CellToSystem',
#         group.by = 'Sample',
#         pt.size=0,log = T)+NoLegend()
# pneum.cts.sub <- subset(cell.to.system,nFeature_CellToSystem>100) # Choose this limit based on the above violin plots or similar
# pneum.cts.sub
# save(pneum.cts.sub,file = 'pneum.cell.to.system.sub.Robj')
# gc()


