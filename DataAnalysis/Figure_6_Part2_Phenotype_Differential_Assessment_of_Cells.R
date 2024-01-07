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

# Load  phenotype data
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Tuft_Sox9_Pneumonectomy_Project/Pneumonectomy_Single_Cell/pneum.clean.annotated.2023-06-14.Robj")

# Format metadata for plotting
pneum.clean$Timepoint <- factor(pneum.clean$Timepoint,
                                levels = c('Day 0','Day 3','Day 7','Day 14'))
pneum.clean$Sample <- factor(pneum.clean$Sample,
                             levels = c('P0_f','P0_m','P3_f','P3_m',
                                        'P7_f','P7_m','P14_f','P14_m'))
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

# Initialize
pneum.clean.data.split <- list()

# Male
pneum.clean.data.split[['male']] <- list()
for (i in 1:length(CellTypes)){
  print(i)
  pneum.clean.data.split[['male']][[CellTypes.mod[i]]] <- pneum.clean.data[['male']]@assays$RNA@data[features,rownames(metadata.male.split[[CellTypes[i]]])]
}

# Female
pneum.clean.data.split[['female']] <- list()
for (i in 1:length(CellTypes)){
  print(i)
  pneum.clean.data.split[['female']][[CellTypes.mod[i]]] <- pneum.clean.data[['female']]@assays$RNA@data[features,rownames(metadata.female.split[[CellTypes[i]]])]
}

# Reduce size of metadata to add
metadata.male.limited <- metadata.male#[,c('orig.ident','CellType','Timepoint')]
metadata.female.limited <- metadata.female#[,c('orig.ident','CellType','Timepoint')]

# Convert to a seurat object and add metadata
# Initialize
pneum.clean.data.split.seurat <- list()
# Male
pneum.clean.data.split.seurat[['male']] <- list()
for (i in 1:length(CellTypes)){
  print(i)
  pneum.clean.data.split.seurat[['male']][[CellTypes.mod[i]]] <- CreateSeuratObject(counts = pneum.clean.data.split[['male']][[CellTypes.mod[i]]],
                                                                                         assay='RNA',
                                                                                         meta.data = metadata.male.limited[colnames(pneum.clean.data.split[['male']][[CellTypes.mod[i]]]),] )
}
# Female
pneum.clean.data.split.seurat[['female']] <- list()
for (i in 1:length(CellTypes)){
  print(i)
  pneum.clean.data.split.seurat[['female']][[CellTypes.mod[i]]] <- CreateSeuratObject(counts = pneum.clean.data.split[['female']][[CellTypes.mod[i]]],
                                                                                           assay='RNA',
                                                                                           meta.data = metadata.female.limited[colnames(pneum.clean.data.split[['female']][[CellTypes.mod[i]]]),] )
}

# Save point
#save(pneum.clean.data.split.seurat,file = 'pneum.clean.data.split.seurat.by.timepoint.unfiltered.Robj')

# Calculate markers for timecourse
plan("multisession")
pneum.clean.data.split.seurat.markers <- list()

# Male
pneum.clean.data.split.seurat.markers[['male']] <- list()
for (i in 1:length(CellTypes)){
  print(i)
  # Set Idents to Timepoint
  Idents(pneum.clean.data.split.seurat[['male']][[CellTypes.mod[i]]]) <- 
    pneum.clean.data.split.seurat[['male']][[CellTypes.mod[i]]]$Timepoint
  # Calculate all marker values with reasonable thresholds
  pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]] <-
    FindAllMarkers(object = pneum.clean.data.split.seurat[['male']][[CellTypes.mod[i]]],
                   features = features,
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
  pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$Sex <- 'Male'
  pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$CellType <- CellTypes[i]
}

# Female
pneum.clean.data.split.seurat.markers[['female']] <- list()
for (i in 1:length(CellTypes)){
  print(i)
  # Set Idents to Timepoint
  Idents(pneum.clean.data.split.seurat[['female']][[CellTypes.mod[i]]]) <- 
    pneum.clean.data.split.seurat[['female']][[CellTypes.mod[i]]]$Timepoint
  # Calculate all marker values with reasonable thresholds
  pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]] <-
    FindAllMarkers(object = pneum.clean.data.split.seurat[['female']][[CellTypes.mod[i]]],
                   features = features,
                   # min.pct = 0.1,
                   # logfc.threshold = 0.1,
                   # min.cells.feature = 3,
                   # min.cells.group = 3,
                   # return.thresh = 0.01,
                   min.pct = 0,
                   logfc.threshold = 0,
                   min.cells.feature = 3,
                   min.cells.group = 3, # Need to adjust to "1" for i=29 only (Rag1+ cells) to get to run
                   return.thresh = 1,
                   only.pos = FALSE)
  # Add additional metadata
  if(nrow(pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]])>0){
    pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$Sex <- 'Female'
    pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$CellType <- CellTypes[i]
  }
}

# Make a unique column ID to identify agreement in significant perturbation between the sexes
for (i in 1:length(CellTypes)){
  print(i)
  pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$unique <- 
    paste(pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$cluster,
          pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$CellType,
          pneum.clean.data.split.seurat.markers[['male']][[CellTypes.mod[i]]]$gene)
  pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$unique <- 
    paste(pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$cluster,
          pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$CellType,
          pneum.clean.data.split.seurat.markers[['female']][[CellTypes.mod[i]]]$gene)
}

# Save point
#save(pneum.clean.data.split.seurat.markers,file = 'pneum.clean.data.split.seurat.markers.Robj')
save(pneum.clean.data.split.seurat.markers,file = 'pneum.clean.data.split.seurat.markers.unfiltered.Robj')





