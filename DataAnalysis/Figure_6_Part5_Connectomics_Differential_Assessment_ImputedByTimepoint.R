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

# Load NICHES data
load('cell.to.cell.imputed.by.timepoint.Robj') # note that this has not been filtered or cleaned - best input for statistical assessment done below

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

# Initialize
niches.imputed.data.split <- list()

# Male
niches.imputed.data.split[['male']] <- list()
for (i in 1:length(VectorTypes)){
  print(i)
  niches.imputed.data.split[['male']][[VectorTypes.mod[i]]] <- niches.imputed.data[['male']]@assays$CellToCell@data[mechs,rownames(metadata.male.split[[VectorTypes[i]]])]
}

# Female
niches.imputed.data.split[['female']] <- list()
for (i in 1:length(VectorTypes)){
  print(i)
  niches.imputed.data.split[['female']][[VectorTypes.mod[i]]] <- niches.imputed.data[['female']]@assays$CellToCell@data[mechs,rownames(metadata.female.split[[VectorTypes[i]]])]
}

# Reduce size of metadata to add
metadata.male.limited <- metadata.male[,c('orig.ident','VectorType','Timepoint.Sending')]
metadata.female.limited <- metadata.female[,c('orig.ident','VectorType','Timepoint.Sending')]

# Convert to a seurat object and add metadata
# Initialize
niches.imputed.data.split.seurat <- list()
# Male
niches.imputed.data.split.seurat[['male']] <- list()
for (i in 1:length(VectorTypes)){
  print(i)
  niches.imputed.data.split.seurat[['male']][[VectorTypes.mod[i]]] <- CreateSeuratObject(counts = niches.imputed.data.split[['male']][[VectorTypes.mod[i]]],
                                                                                         assay='CellToCell',
                                                                                         meta.data = metadata.male.limited[colnames(niches.imputed.data.split[['male']][[VectorTypes.mod[i]]]),] )
}
# Female
niches.imputed.data.split.seurat[['female']] <- list()
for (i in 1:length(VectorTypes)){
  print(i)
  niches.imputed.data.split.seurat[['female']][[VectorTypes.mod[i]]] <- CreateSeuratObject(counts = niches.imputed.data.split[['female']][[VectorTypes.mod[i]]],
                                                                                           assay='CellToCell',
                                                                                           meta.data = metadata.female.limited[colnames(niches.imputed.data.split[['female']][[VectorTypes.mod[i]]]),] )
}

# Save point
save(niches.imputed.data.split.seurat,file = 'niches.imputed.data.split.seurat.by.timepoint.unfiltered.Robj')

# Calculate markers for timecourse
niches.imputed.data.split.seurat.markers <- list()

# Male
niches.imputed.data.split.seurat.markers[['male']] <- list()
for (i in 1:length(VectorTypes)){
  print(i)
  # Set Idents to Timepoint
  Idents(niches.imputed.data.split.seurat[['male']][[VectorTypes.mod[i]]]) <- 
    niches.imputed.data.split.seurat[['male']][[VectorTypes.mod[i]]]$Timepoint.Sending
  # Calculate all marker values with reasonable thresholds
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]] <-
    FindAllMarkers(object = niches.imputed.data.split.seurat[['male']][[VectorTypes.mod[i]]],
                   features = mechs,
                   # min.pct = 0.1,
                   # logfc.threshold = 0.1,
                   # min.cells.feature = 3,
                   # min.cells.group = 3,
                   # return.thresh = 0.01,
                   min.pct = 0,
                   logfc.threshold = 0,
                   min.cells.feature = 3,
                   min.cells.group = 3,
                   return.thresh = 1,
                   only.pos = FALSE)
  # Add additional metadata
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$Sex <- 'Male'
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$VectorType <- VectorTypes[i]
}

# Female
niches.imputed.data.split.seurat.markers[['female']] <- list()
for (i in 1:length(VectorTypes)){
  print(i)
  # Set Idents to Timepoint
  Idents(niches.imputed.data.split.seurat[['female']][[VectorTypes.mod[i]]]) <- 
    niches.imputed.data.split.seurat[['female']][[VectorTypes.mod[i]]]$Timepoint.Sending
  # Calculate all marker values with reasonable thresholds
  niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]] <-
    FindAllMarkers(object = niches.imputed.data.split.seurat[['female']][[VectorTypes.mod[i]]],
                   features = mechs,
                   # min.pct = 0.1,
                   # logfc.threshold = 0.1,
                   # min.cells.feature = 3,
                   # min.cells.group = 3,
                   # return.thresh = 0.01,
                   min.pct = 0,
                   logfc.threshold = 0,
                   min.cells.feature = 3,
                   min.cells.group = 3,
                   return.thresh = 1,
                   only.pos = FALSE)
  # Add additional metadata
  if(nrow(niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]])>0){
    niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$Sex <- 'Female'
    niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$VectorType <- VectorTypes[i]
  }
}

# Make a unique column ID to identify agreement in significant perturbation between the sexes
for (i in 1:length(VectorTypes)){
  print(i)
  niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$unique <- 
    paste(niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$cluster,
          niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$VectorType,
          niches.imputed.data.split.seurat.markers[['male']][[VectorTypes.mod[i]]]$gene)
  niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$unique <- 
    paste(niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$cluster,
          niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$VectorType,
          niches.imputed.data.split.seurat.markers[['female']][[VectorTypes.mod[i]]]$gene)
}

# Save point
#save(niches.imputed.data.split.seurat.markers,file = 'niches.imputed.data.split.seurat.markers.Robj')
save(niches.imputed.data.split.seurat.markers,file = 'niches.imputed.data.split.seurat.markers.by.timepoint.unfiltered.Robj')



