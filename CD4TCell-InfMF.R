################################################################################
#Project full name: Human memory CD4+ T-cells recognize Mycobacterium tuberculosis-infected macrophages amid broader pathogen-specific responses
#Project abbreviated name: CD4TCell-InfMF
################################################################################

###The following R script consists of the followings parts:

# I. Preintegration processing
# 1. Preparing the session
# 2. Loading data sets and transforming them into Seurat objects
# 3. Adjusting the meta.data structure
# 4. Adding cell barcodes to TCR read tables
# 5. Quantify percentages of mitochondrial, ribosomal, TCR genes
# 6. Quality control: removal low quality cells
# 7. Cell cycle scoring
# 8. Combining scRNA- and scTCR-seq data data
# 9. Labeling of clonally expanded cells
# 10. Merging samples
# 11. Variance stabilization

# II. Sample integration/Batch correction

# III. Postintegration analysis
# 1. Cell clustering
# 2. Quality control: removal of monocytic contamination
# 3. Gene expression analysis
# 4. Assessment of clonal expansion
# 5. Single-cell trajectory analysis
# 6. Mapping CDR3 motifs on the UMAP
# 7. Differential gene expression analysis
# 8. Biological theme comparison/Reactome pathway overrepresentation analysis
# 9. Cell-cell communication analysis

################################################################################
# I. Preintegration processing #################################################
################################################################################

#=== 1. Preparing the session ==================================================

librarian::lib_startup(librarian, global = T)
librarian::shelf(Seurat, BiocManager, devtools, scRepertoire, BiocFileCache, speedglm, BiocGenerics, 
                 tidyverse, sctransform, gridExtra, patchwork, SingleCellExperiment, cowplot, monocle3, 
                 SeuratObject, DelayedArray, DelayedMatrixStats, limma, lme4, S4Vectors, batchelor, 
                 HDF5Array, terra, ggrastr, SummarizedExperiment, Trex, S4Vectors, monocle3, Matrix,
                 multtest, metap, fields, KernSmooth, ROCR, parallel, SingleR, SoupX, reticulate,
                 dplyr, ggplot2, tidyverse, hdf5r, patchwork, BPCells, fs, glmGamPoi, DropletUtils,
                 rstan, SeuratWrappers, DESeq2, MAST, EnhancedVolcano, clusterProfiler, qs, scater,
                 AnnotationDbi, readxl, clustree, viridis, scCustomize, tidyseurat, circlize, scales,
                 celldex, pheatmap, scDblFinder, irlba, data.table, SCpubr, ComplexHeatmap, Nebulosa,
                 speckle, statmod, edgeR, CellBench, RColorBrewer, msigdbr, SPIA, fgsea, svglite, 
                 org.Hs.eg.db, enrichplot, ggpubr, ggupset, ReactomePA, reactome.db, future, dittoSeq,
                 nichenetr, promises, AUCell, doMC, doRNG, doSNOW, mixtools, DT, plotly, NMF, d3heatmap,
                 dynamicTreeCut, R2HTML, Rtsne, zoo, magic, update_all = T, quiet = T)

set.seed(123)

#=== 2. Loading data sets and transforming them into Seurat objects ============

# Setting up the directory folder
dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/data/CL_Lys")

# Loading the original combined GEX data for Inf vs Inf+Lysate experiment 
CL003_Inf_LYSATErun.data <- Read10X_h5(filename = 'CL003_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL003_InfLys.data <- Read10X_h5(filename = 'CL003_InfLys_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL007_Inf_LYSATErun.data <- Read10X_h5(filename = 'CL007_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL007_InfLys.data <- Read10X_h5(filename = 'CL007_InfLys_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL011_Inf_LYSATErun.data <- Read10X_h5(filename = 'CL011_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL011_InfLys.data <- Read10X_h5(filename = 'CL011_InfLys_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL012_Inf_LYSATErun.data <- Read10X_h5(filename = 'CL012_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL012_InfLys.data <- Read10X_h5(filename = 'CL012_InfLys_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL013_Inf_LYSATErun.data <- Read10X_h5(filename = 'CL013_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL013_InfLys.data <- Read10X_h5(filename = 'CL013_InfLys_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL018_Inf_LYSATErun.data <- Read10X_h5(filename = 'CL018_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
P079_Inf_LYSATErun.data <- Read10X_h5(filename = 'P079_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)

# Setting up the directory folder for non-LTBI controls

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/data/Non-LTBI-InfOnly")

# Loading the original combined GEX data from non-LTBI controls with InfOnly

CL002_Inf.data <- Read10X_h5(filename = 'CL002ba_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
CL010_Inf.data <- Read10X_h5(filename = 'CL010_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
Lkpk2_Inf.data <- Read10X_h5(filename = 'Lkpk2_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
Lkpk3_Inf.data <- Read10X_h5(filename = 'Lkpk3_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
Lkpk4_Inf.data <- Read10X_h5(filename = 'Lkpk4_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)
Lkpk5_Inf.data <- Read10X_h5(filename = 'Lkpk5_Inf_sample_filtered_feature_bc_matrix.h5', use.names = T, unique.features = T)

# Setting up the directory folder for Lysate run

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/data/CL_Lys")

# Read in TCR clonotype read counts from Cellranger VDJ output

CL003_Inf_LYSATErun_clones <- read.csv('CL003_Inf_filtered_contig_annotations.csv')
CL003_InfLys_clones <- read.csv('CL003_InfLys_filtered_contig_annotations.csv')
CL007_Inf_LYSATErun_clones <- read.csv('CL007_Inf_filtered_contig_annotations.csv')
CL007_InfLys_clones <- read.csv('CL007_InfLys_filtered_contig_annotations.csv')
CL011_Inf_LYSATErun_clones <- read.csv('CL011_Inf_filtered_contig_annotations.csv')
CL011_InfLys_clones <- read.csv('CL011_InfLys_filtered_contig_annotations.csv')
CL012_Inf_LYSATErun_clones <- read.csv('CL012_Inf_filtered_contig_annotations.csv')
CL012_InfLys_clones <- read.csv('CL012_InfLys_filtered_contig_annotations.csv')
CL013_Inf_LYSATErun_clones <- read.csv('CL013_Inf_filtered_contig_annotations.csv')
CL013_InfLys_clones <- read.csv('CL013_InfLys_filtered_contig_annotations.csv')
CL018_Inf_LYSATErun_clones <- read.csv('CL018_Inf_filtered_contig_annotations.csv')
P079_Inf_LYSATErun_clones <- read.csv('P079_Inf_filtered_contig_annotations.csv')

# Setting up the directory folder for non-LTBI controls

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/data/Non-LTBI-InfOnly")

# Read in TCR clonotype read counts from Cellranger VDJ output for non-LTBI controls

CL002_Inf_clones <- read.csv('CL002ba_Inf_filtered_contig_annotations.csv')
CL010_Inf_clones <- read.csv('CL010_Inf_filtered_contig_annotations.csv')
Lkpk2_Inf_clones <- read.csv('Lkpk2_Inf_filtered_contig_annotations.csv')
Lkpk3_Inf_clones <- read.csv('Lkpk3_Inf_filtered_contig_annotations.csv')
Lkpk4_Inf_clones <- read.csv('Lkpk4_Inf_filtered_contig_annotations.csv')
Lkpk5_Inf_clones <- read.csv('Lkpk5_Inf_filtered_contig_annotations.csv')

# Setting up the directory folder

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")

# Creating Seurat objects from GEX data 

CL003_Inf_LYSATErun <- CreateSeuratObject(CL003_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL003_InfLys <- CreateSeuratObject(CL003_InfLys.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL007_Inf_LYSATErun <- CreateSeuratObject(CL007_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL007_InfLys <- CreateSeuratObject(CL007_InfLys.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL011_Inf_LYSATErun <- CreateSeuratObject(CL011_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL011_InfLys <- CreateSeuratObject(CL011_InfLys.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL012_Inf_LYSATErun <- CreateSeuratObject(CL012_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL012_InfLys <- CreateSeuratObject(CL012_InfLys.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL013_Inf_LYSATErun <- CreateSeuratObject(CL013_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL013_InfLys <- CreateSeuratObject(CL013_InfLys.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL018_Inf_LYSATErun <- CreateSeuratObject(CL018_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
P079_Inf_LYSATErun <- CreateSeuratObject(P079_Inf_LYSATErun.data, project = 'LTBI', min.cells = 3, min.features = 200)
CL002_Inf <- CreateSeuratObject(CL002_Inf.data, project = 'nonLTBI', min.cells = 3, min.features = 200)
CL010_Inf <- CreateSeuratObject(CL010_Inf.data, project = 'nonLTBI', min.cells = 3, min.features = 200)
Lkpk2_Inf <- CreateSeuratObject(Lkpk2_Inf.data, project = 'nonLTBI', min.cells = 3, min.features = 200)
Lkpk3_Inf <- CreateSeuratObject(Lkpk3_Inf.data, project = 'nonLTBI', min.cells = 3, min.features = 200)
Lkpk4_Inf <- CreateSeuratObject(Lkpk4_Inf.data, project = 'nonLTBI', min.cells = 3, min.features = 200)
Lkpk5_Inf <- CreateSeuratObject(Lkpk5_Inf.data, project = 'nonLTBI', min.cells = 3, min.features = 200)

#--- Delete files that are no longer needed

rm(CL003_Inf_LYSATErun.data, CL003_InfLys.data,
   CL007_Inf_LYSATErun.data, CL007_InfLys.data,
   CL011_Inf_LYSATErun.data, CL011_InfLys.data,
   CL012_Inf_LYSATErun.data, CL012_InfLys.data,
   CL013_Inf_LYSATErun.data, CL013_InfLys.data,
   CL018_Inf_LYSATErun.data, P079_Inf_LYSATErun.data,
   CL002_Inf.data, CL010_Inf.data, 
   Lkpk2_Inf.data, Lkpk3_Inf.data,
   Lkpk4_Inf.data, Lkpk5_Inf.data)

#=== 3. Adjusting the meta.data structure ======================================

# Adding prefix to cell barcodes corresponding to sample of origin to prepare for scRepertoire compatibility

CL003_Inf_LYSATErun <- RenameCells(CL003_Inf_LYSATErun, add.cell.id = 'CL003_LTBI:Inf')
CL003_InfLys <- RenameCells(CL003_InfLys, add.cell.id = 'CL003_LTBI:InfLys')
CL007_Inf_LYSATErun <- RenameCells(CL007_Inf_LYSATErun, add.cell.id = 'CL007_LTBI:Inf')
CL007_InfLys <- RenameCells(CL007_InfLys, add.cell.id = 'CL007_LTBI:InfLys')
CL011_Inf_LYSATErun <- RenameCells(CL011_Inf_LYSATErun, add.cell.id = 'CL011_LTBI:Inf')
CL011_InfLys <- RenameCells(CL011_InfLys, add.cell.id = 'CL011_LTBI:InfLys')
CL012_Inf_LYSATErun <- RenameCells(CL012_Inf_LYSATErun, add.cell.id = 'CL012_LTBI:Inf')
CL012_InfLys <- RenameCells(CL012_InfLys, add.cell.id = 'CL012_LTBI:InfLys')
CL013_Inf_LYSATErun <- RenameCells(CL013_Inf_LYSATErun, add.cell.id = 'CL013_LTBI:Inf')
CL013_InfLys <- RenameCells(CL013_InfLys, add.cell.id = 'CL013_LTBI:InfLys')
CL018_Inf_LYSATErun <- RenameCells(CL018_Inf_LYSATErun, add.cell.id = 'CL018_LTBI:Inf')
P079_Inf_LYSATErun <- RenameCells(P079_Inf_LYSATErun, add.cell.id = 'P079_LTBI:Inf')
CL002_Inf <- RenameCells(CL002_Inf, add.cell.id = 'CL002_nonLTBI:Inf')
CL010_Inf <- RenameCells(CL010_Inf, add.cell.id = 'CL010_nonLTBI:Inf')
Lkpk2_Inf <- RenameCells(Lkpk2_Inf, add.cell.id = 'Lkpk2_nonLTBI:Inf')
Lkpk3_Inf <- RenameCells(Lkpk3_Inf, add.cell.id = 'Lkpk3_nonLTBI:Inf')
Lkpk4_Inf <- RenameCells(Lkpk4_Inf, add.cell.id = 'Lkpk4_nonLTBI:Inf')
Lkpk5_Inf <- RenameCells(Lkpk5_Inf, add.cell.id = 'Lkpk5_nonLTBI:Inf')

# Reorganizing the meta.data in more intuitive way

CL003_Inf_LYSATErun$sample <- rownames(CL003_Inf_LYSATErun@meta.data)
CL003_InfLys$sample <- rownames(CL003_InfLys@meta.data)
CL007_Inf_LYSATErun$sample <- rownames(CL007_Inf_LYSATErun@meta.data)
CL007_InfLys$sample <- rownames(CL007_InfLys@meta.data)
CL011_Inf_LYSATErun$sample <- rownames(CL011_Inf_LYSATErun@meta.data)
CL011_InfLys$sample <- rownames(CL011_InfLys@meta.data)
CL012_Inf_LYSATErun$sample <- rownames(CL012_Inf_LYSATErun@meta.data)
CL012_InfLys$sample <- rownames(CL012_InfLys@meta.data)
CL013_Inf_LYSATErun$sample <- rownames(CL013_Inf_LYSATErun@meta.data)
CL013_InfLys$sample <- rownames(CL013_InfLys@meta.data)
CL018_Inf_LYSATErun$sample <- rownames(CL018_Inf_LYSATErun@meta.data)
P079_Inf_LYSATErun$sample <- rownames(P079_Inf_LYSATErun@meta.data)
CL002_Inf$sample <- rownames(CL002_Inf@meta.data)
CL010_Inf$sample <- rownames(CL010_Inf@meta.data)
Lkpk2_Inf$sample <- rownames(Lkpk2_Inf@meta.data)
Lkpk3_Inf$sample <- rownames(Lkpk3_Inf@meta.data)
Lkpk4_Inf$sample <- rownames(Lkpk4_Inf@meta.data)
Lkpk5_Inf$sample <- rownames(Lkpk5_Inf@meta.data)

# Modifying the meta.data structure for the analysis

CL003_Inf_LYSATErun@meta.data <- tidyr::separate(CL003_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL003_Inf_LYSATErun@meta.data <- unite(CL003_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL003_Inf_LYSATErun@meta.data$Barcode <- NULL
CL003_Inf_LYSATErun@meta.data$barcode <- NULL

CL003_InfLys@meta.data <- tidyr::separate(CL003_InfLys@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL003_InfLys@meta.data <- unite(CL003_InfLys@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL003_InfLys@meta.data$Barcode <- NULL
CL003_InfLys@meta.data$barcode <- NULL

CL007_Inf_LYSATErun@meta.data <- tidyr::separate(CL007_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL007_Inf_LYSATErun@meta.data <- unite(CL007_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL007_Inf_LYSATErun@meta.data$Barcode <- NULL
CL007_Inf_LYSATErun@meta.data$barcode <- NULL

CL007_InfLys@meta.data <- tidyr::separate(CL007_InfLys@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL007_InfLys@meta.data <- unite(CL007_InfLys@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL007_InfLys@meta.data$Barcode <- NULL
CL007_InfLys@meta.data$barcode <- NULL

CL011_Inf_LYSATErun@meta.data <- tidyr::separate(CL011_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL011_Inf_LYSATErun@meta.data <- unite(CL011_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL011_Inf_LYSATErun@meta.data$Barcode <- NULL
CL011_Inf_LYSATErun@meta.data$barcode <- NULL

CL011_InfLys@meta.data <- tidyr::separate(CL011_InfLys@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL011_InfLys@meta.data <- unite(CL011_InfLys@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL011_InfLys@meta.data$Barcode <- NULL
CL011_InfLys@meta.data$barcode <- NULL

CL012_Inf_LYSATErun@meta.data <- tidyr::separate(CL012_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL012_Inf_LYSATErun@meta.data <- unite(CL012_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL012_Inf_LYSATErun@meta.data$Barcode <- NULL
CL012_Inf_LYSATErun@meta.data$barcode <- NULL

CL012_InfLys@meta.data <- tidyr::separate(CL012_InfLys@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL012_InfLys@meta.data <- unite(CL012_InfLys@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL012_InfLys@meta.data$Barcode <- NULL
CL012_InfLys@meta.data$barcode <- NULL

CL013_Inf_LYSATErun@meta.data <- tidyr::separate(CL013_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL013_Inf_LYSATErun@meta.data <- unite(CL013_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL013_Inf_LYSATErun@meta.data$Barcode <- NULL
CL013_Inf_LYSATErun@meta.data$barcode <- NULL

CL013_InfLys@meta.data <- tidyr::separate(CL013_InfLys@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL013_InfLys@meta.data <- unite(CL013_InfLys@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL013_InfLys@meta.data$Barcode <- NULL
CL013_InfLys@meta.data$barcode <- NULL

CL018_Inf_LYSATErun@meta.data <- tidyr::separate(CL018_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL018_Inf_LYSATErun@meta.data <- unite(CL018_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL018_Inf_LYSATErun@meta.data$Barcode <- NULL
CL018_Inf_LYSATErun@meta.data$barcode <- NULL

P079_Inf_LYSATErun@meta.data <- tidyr::separate(P079_Inf_LYSATErun@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
P079_Inf_LYSATErun@meta.data <- unite(P079_Inf_LYSATErun@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
P079_Inf_LYSATErun@meta.data$Barcode <- NULL
P079_Inf_LYSATErun@meta.data$barcode <- NULL

CL002_Inf@meta.data <- tidyr::separate(CL002_Inf@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL002_Inf@meta.data <- unite(CL002_Inf@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL002_Inf@meta.data$Barcode <- NULL
CL002_Inf@meta.data$barcode <- NULL

CL010_Inf@meta.data <- tidyr::separate(CL010_Inf@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
CL010_Inf@meta.data <- unite(CL010_Inf@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
CL010_Inf@meta.data$Barcode <- NULL
CL010_Inf@meta.data$barcode <- NULL

Lkpk2_Inf@meta.data <- tidyr::separate(Lkpk2_Inf@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
Lkpk2_Inf@meta.data <- unite(Lkpk2_Inf@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
Lkpk2_Inf@meta.data$Barcode <- NULL
Lkpk2_Inf@meta.data$barcode <- NULL

Lkpk3_Inf@meta.data <- tidyr::separate(Lkpk3_Inf@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
Lkpk3_Inf@meta.data <- unite(Lkpk3_Inf@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
Lkpk3_Inf@meta.data$Barcode <- NULL
Lkpk3_Inf@meta.data$barcode <- NULL

Lkpk4_Inf@meta.data <- tidyr::separate(Lkpk4_Inf@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
Lkpk4_Inf@meta.data <- unite(Lkpk4_Inf@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
Lkpk4_Inf@meta.data$Barcode <- NULL
Lkpk4_Inf@meta.data$barcode <- NULL

Lkpk5_Inf@meta.data <- tidyr::separate(Lkpk5_Inf@meta.data, col = 'sample', into = c('Donor', 'Condition', 'Barcode'), sep = '_', remove = TRUE)
Lkpk5_Inf@meta.data <- unite(Lkpk5_Inf@meta.data, col = Sample, c('Donor', 'Condition'), remove = FALSE, sep = '_', na.rm = FALSE)
Lkpk5_Inf@meta.data$Barcode <- NULL
Lkpk5_Inf@meta.data$barcode <- NULL

#=== 4. Adding cell barcodes to TCR read tables ================================

CL003_Inf_tcr_LYSATErun <- combineTCR(CL003_Inf_LYSATErun_clones, samples = c('CL003'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
CL003_InfLys_tcr <- combineTCR(CL003_InfLys_clones, samples = c('CL003'), ID = c('LTBI:InfLys'), removeNA = T, filterMulti = T)
CL007_Inf_tcr_LYSATErun <- combineTCR(CL007_Inf_LYSATErun_clones, samples = c('CL007'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
CL007_InfLys_tcr <- combineTCR(CL007_InfLys_clones, samples = c('CL007'), ID = c('LTBI:InfLys'), removeNA = T, filterMulti = T)
CL011_Inf_tcr_LYSATErun <- combineTCR(CL011_Inf_LYSATErun_clones, samples = c('CL011'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
CL011_InfLys_tcr <- combineTCR(CL011_InfLys_clones, samples = c('CL011'), ID = c('LTBI:InfLys'), removeNA = T, filterMulti = T)
CL012_Inf_tcr_LYSATErun <- combineTCR(CL012_Inf_LYSATErun_clones, samples = c('CL012'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
CL012_InfLys_tcr <- combineTCR(CL012_InfLys_clones, samples = c('CL012'), ID = c('LTBI:InfLys'), removeNA = T, filterMulti = T)
CL013_Inf_tcr_LYSATErun <- combineTCR(CL013_Inf_LYSATErun_clones, samples = c('CL013'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
CL013_InfLys_tcr <- combineTCR(CL013_InfLys_clones, samples = c('CL013'), ID = c('LTBI:InfLys'), removeNA = T, filterMulti = T)
CL018_Inf_tcr_LYSATErun <- combineTCR(CL018_Inf_LYSATErun_clones, samples = c('CL018'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
P079_Inf_tcr_LYSATErun <- combineTCR(P079_Inf_LYSATErun_clones, samples = c('P079'), ID = c('LTBI:Inf'), removeNA = T, filterMulti = T)
CL002_Inf_tcr <- combineTCR(CL002_Inf_clones, samples = c('CL002'), ID = c('nonLTBI:Inf'), removeNA = T, filterMulti = T)
CL010_Inf_tcr <- combineTCR(CL010_Inf_clones, samples = c('CL010'), ID = c('nonLTBI:Inf'), removeNA = T, filterMulti = T)
Lkpk2_Inf_tcr <- combineTCR(Lkpk2_Inf_clones, samples = c('Lkpk2'), ID = c('nonLTBI:Inf'), removeNA = T, filterMulti = T)
Lkpk3_Inf_tcr <- combineTCR(Lkpk3_Inf_clones, samples = c('Lkpk3'), ID = c('nonLTBI:Inf'), removeNA = T, filterMulti = T)
Lkpk4_Inf_tcr <- combineTCR(Lkpk4_Inf_clones, samples = c('Lkpk4'), ID = c('nonLTBI:Inf'), removeNA = T, filterMulti = T)
Lkpk5_Inf_tcr <- combineTCR(Lkpk5_Inf_clones, samples = c('Lkpk5'), ID = c('nonLTBI:Inf'), removeNA = T, filterMulti = T)

#--- Delete files that are no longer needed

rm(CL003_Inf_LYSATErun_clones, CL003_InfLys_clones,
   CL007_Inf_LYSATErun_clones, CL007_InfLys_clones,
   CL011_Inf_LYSATErun_clones, CL011_InfLys_clones,
   CL012_Inf_LYSATErun_clones, CL012_InfLys_clones,
   CL013_Inf_LYSATErun_clones, CL013_InfLys_clones,
   CL018_Inf_LYSATErun_clones, P079_Inf_LYSATErun_clones,
   CL002_Inf_clones, CL010_Inf_clones,
   Lkpk2_Inf_clones, Lkpk3_Inf_clones,
   Lkpk4_Inf_clones, Lkpk5_Inf_clones)

#=== 5. Quantify percentages of mitochondrial, ribosomal, TCR genes ============

#--- % mitochondrial genes reads

CL003_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(CL003_Inf_LYSATErun, pattern = "^MT-")
CL003_InfLys[['percent.mt']] <- PercentageFeatureSet(CL003_InfLys, pattern = "^MT-")
CL007_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(CL007_Inf_LYSATErun, pattern = "^MT-")
CL007_InfLys[['percent.mt']] <- PercentageFeatureSet(CL007_InfLys, pattern = "^MT-")
CL011_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(CL011_Inf_LYSATErun, pattern = "^MT-")
CL011_InfLys[['percent.mt']] <- PercentageFeatureSet(CL011_InfLys, pattern = "^MT-")
CL012_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(CL012_Inf_LYSATErun, pattern = "^MT-")
CL012_InfLys[['percent.mt']] <- PercentageFeatureSet(CL012_InfLys, pattern = "^MT-")
CL013_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(CL013_Inf_LYSATErun, pattern = "^MT-")
CL013_InfLys[['percent.mt']] <- PercentageFeatureSet(CL013_InfLys, pattern = "^MT-")
CL018_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(CL018_Inf_LYSATErun, pattern = "^MT-")
P079_Inf_LYSATErun[['percent.mt']] <- PercentageFeatureSet(P079_Inf_LYSATErun, pattern = "^MT-")
CL002_Inf[['percent.mt']] <- PercentageFeatureSet(CL002_Inf, pattern = "^MT-")
CL010_Inf[['percent.mt']] <- PercentageFeatureSet(CL010_Inf, pattern = "^MT-")
Lkpk2_Inf[['percent.mt']] <- PercentageFeatureSet(Lkpk2_Inf, pattern = "^MT-")
Lkpk3_Inf[['percent.mt']] <- PercentageFeatureSet(Lkpk3_Inf, pattern = "^MT-")
Lkpk4_Inf[['percent.mt']] <- PercentageFeatureSet(Lkpk4_Inf, pattern = "^MT-")
Lkpk5_Inf[['percent.mt']] <- PercentageFeatureSet(Lkpk5_Inf, pattern = "^MT-")

#--- % ribosomal genes reads

CL003_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(CL003_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL003_InfLys[['percent.rb']] <- PercentageFeatureSet(CL003_InfLys, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL007_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(CL007_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL007_InfLys[['percent.rb']] <- PercentageFeatureSet(CL007_InfLys, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL011_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(CL011_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL011_InfLys[['percent.rb']] <- PercentageFeatureSet(CL011_InfLys, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL012_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(CL012_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL012_InfLys[['percent.rb']] <- PercentageFeatureSet(CL012_InfLys, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL013_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(CL013_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL013_InfLys[['percent.rb']] <- PercentageFeatureSet(CL013_InfLys, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL018_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(CL018_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
P079_Inf_LYSATErun[['percent.rb']] <- PercentageFeatureSet(P079_Inf_LYSATErun, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL002_Inf[['percent.rb']] <- PercentageFeatureSet(CL002_Inf, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
CL010_Inf[['percent.rb']] <- PercentageFeatureSet(CL010_Inf, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
Lkpk2_Inf[['percent.rb']] <- PercentageFeatureSet(Lkpk2_Inf, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
Lkpk3_Inf[['percent.rb']] <- PercentageFeatureSet(Lkpk3_Inf, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
Lkpk4_Inf[['percent.rb']] <- PercentageFeatureSet(Lkpk4_Inf, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
Lkpk5_Inf[['percent.rb']] <- PercentageFeatureSet(Lkpk5_Inf, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

#--- % TCR genes reads

CL003_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(CL003_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL003_InfLys[['percent.tcr']] <- PercentageFeatureSet(CL003_InfLys, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL007_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(CL007_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL007_InfLys[['percent.tcr']] <- PercentageFeatureSet(CL007_InfLys, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*]")
CL011_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(CL011_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*]")
CL011_InfLys[['percent.tcr']] <- PercentageFeatureSet(CL011_InfLys, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL012_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(CL012_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL012_InfLys[['percent.tcr']] <- PercentageFeatureSet(CL012_InfLys, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL013_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(CL013_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL013_InfLys[['percent.tcr']] <- PercentageFeatureSet(CL013_InfLys, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL018_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(CL018_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
P079_Inf_LYSATErun[['percent.tcr']] <- PercentageFeatureSet(P079_Inf_LYSATErun, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL002_Inf[['percent.tcr']] <- PercentageFeatureSet(CL002_Inf, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
CL010_Inf[['percent.tcr']] <- PercentageFeatureSet(CL010_Inf, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
Lkpk2_Inf[['percent.tcr']] <- PercentageFeatureSet(Lkpk2_Inf, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
Lkpk3_Inf[['percent.tcr']] <- PercentageFeatureSet(Lkpk3_Inf, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
Lkpk4_Inf[['percent.tcr']] <- PercentageFeatureSet(Lkpk4_Inf, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")
Lkpk5_Inf[['percent.tcr']] <- PercentageFeatureSet(Lkpk5_Inf, pattern = "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")

#=== 6. Quality control: removal low quality cells =============================

CL003_Inf_LYSATErun <- subset(CL003_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL003_InfLys <- subset(CL003_InfLys, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL007_Inf_LYSATErun <- subset(CL007_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL007_InfLys <- subset(CL007_InfLys, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL011_Inf_LYSATErun <- subset(CL011_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL011_InfLys <- subset(CL011_InfLys, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL012_Inf_LYSATErun <- subset(CL012_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL012_InfLys <- subset(CL012_InfLys, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL013_Inf_LYSATErun <- subset(CL013_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL013_InfLys <- subset(CL013_InfLys, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL018_Inf_LYSATErun <- subset(CL018_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
P079_Inf_LYSATErun <- subset(P079_Inf_LYSATErun, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL002_Inf <- subset(CL002_Inf, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
CL010_Inf <- subset(CL010_Inf, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
Lkpk2_Inf <- subset(Lkpk2_Inf, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
Lkpk3_Inf <- subset(Lkpk3_Inf, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
Lkpk4_Inf <- subset(Lkpk4_Inf, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
Lkpk5_Inf <- subset(Lkpk5_Inf, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

#=== 7. Cell-Cycle Scoring =====================================================

# Normalize the data and perform cell-cycle scoring

CL003_Inf_LYSATErun <- NormalizeData(CL003_Inf_LYSATErun)
CL003_Inf_LYSATErun <- CellCycleScoring(object = CL003_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
CL003_InfLys <- NormalizeData(CL003_InfLys)
CL003_InfLys <- CellCycleScoring(object = CL003_InfLys, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

CL007_Inf_LYSATErun <- NormalizeData(CL007_Inf_LYSATErun)
CL007_Inf_LYSATErun <- CellCycleScoring(object = CL007_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
CL007_InfLys <- NormalizeData(CL007_InfLys)
CL007_InfLys <- CellCycleScoring(object = CL007_InfLys, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

CL011_Inf_LYSATErun <- NormalizeData(CL011_Inf_LYSATErun)
CL011_Inf_LYSATErun <- CellCycleScoring(object = CL011_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
CL011_InfLys <- NormalizeData(CL011_InfLys)
CL011_InfLys <- CellCycleScoring(object = CL011_InfLys, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

CL012_Inf_LYSATErun <- NormalizeData(CL012_Inf_LYSATErun)
CL012_Inf_LYSATErun <- CellCycleScoring(object = CL012_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
CL012_InfLys <- NormalizeData(CL012_InfLys)
CL012_InfLys <- CellCycleScoring(object = CL012_InfLys, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

CL013_Inf_LYSATErun <- NormalizeData(CL013_Inf_LYSATErun)
CL013_Inf_LYSATErun <- CellCycleScoring(object = CL013_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
CL013_InfLys <- NormalizeData(CL013_InfLys)
CL013_InfLys <- CellCycleScoring(object = CL013_InfLys, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

CL018_Inf_LYSATErun <- NormalizeData(CL018_Inf_LYSATErun)
CL018_Inf_LYSATErun <- CellCycleScoring(object = CL018_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

P079_Inf_LYSATErun <- NormalizeData(P079_Inf_LYSATErun)
P079_Inf_LYSATErun <- CellCycleScoring(object = P079_Inf_LYSATErun, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

CL002_Inf <- NormalizeData(CL002_Inf)
CL002_Inf <- CellCycleScoring(object = CL002_Inf, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
CL010_Inf <- NormalizeData(CL010_Inf)
CL010_Inf <- CellCycleScoring(object = CL010_Inf, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
Lkpk2_Inf <- NormalizeData(Lkpk2_Inf)
Lkpk2_Inf <- CellCycleScoring(object = Lkpk2_Inf, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
Lkpk3_Inf <- NormalizeData(Lkpk3_Inf)
Lkpk3_Inf <- CellCycleScoring(object = Lkpk3_Inf, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
Lkpk4_Inf <- NormalizeData(Lkpk4_Inf)
Lkpk4_Inf <- CellCycleScoring(object = Lkpk4_Inf, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
Lkpk5_Inf <- NormalizeData(Lkpk5_Inf)
Lkpk5_Inf <- CellCycleScoring(object = Lkpk5_Inf, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)

#=== 8. Combining scRNA- and scTCR-seq data data ===============================

CL003_Inf_LYSATErun <- combineExpression(CL003_Inf_tcr_LYSATErun, CL003_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL003_InfLys <- combineExpression(CL003_InfLys_tcr, CL003_InfLys, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL007_Inf_LYSATErun <- combineExpression(CL007_Inf_tcr_LYSATErun, CL007_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL007_InfLys <- combineExpression(CL007_InfLys_tcr, CL007_InfLys, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL011_Inf_LYSATErun <- combineExpression(CL011_Inf_tcr_LYSATErun, CL011_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL011_InfLys <- combineExpression(CL011_InfLys_tcr, CL011_InfLys, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL012_Inf_LYSATErun <- combineExpression(CL012_Inf_tcr_LYSATErun, CL012_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL012_InfLys <- combineExpression(CL012_InfLys_tcr, CL012_InfLys, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL013_Inf_LYSATErun <- combineExpression(CL013_Inf_tcr_LYSATErun, CL013_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL013_InfLys <- combineExpression(CL013_InfLys_tcr, CL013_InfLys, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL018_Inf_LYSATErun <- combineExpression(CL018_Inf_tcr_LYSATErun, CL018_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
P079_Inf_LYSATErun <- combineExpression(P079_Inf_tcr_LYSATErun, P079_Inf_LYSATErun, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL002_Inf <- combineExpression(CL002_Inf_tcr, CL002_Inf, cloneCall = 'gene', chain = "both", proportion = TRUE)
CL010_Inf <- combineExpression(CL010_Inf_tcr, CL010_Inf, cloneCall = 'gene', chain = "both", proportion = TRUE)
Lkpk2_Inf <- combineExpression(Lkpk2_Inf_tcr, Lkpk2_Inf, cloneCall = 'gene', chain = "both", proportion = TRUE)
Lkpk3_Inf <- combineExpression(Lkpk3_Inf_tcr, Lkpk3_Inf, cloneCall = 'gene', chain = "both", proportion = TRUE)
Lkpk4_Inf <- combineExpression(Lkpk4_Inf_tcr, Lkpk4_Inf, cloneCall = 'gene', chain = "both", proportion = TRUE)
Lkpk5_Inf <- combineExpression(Lkpk5_Inf_tcr, Lkpk5_Inf, cloneCall = 'gene', chain = "both", proportion = TRUE)

#=== 9. Labeling of clonally expanded cells ====================================

# Labeling strategy: version 1

CL003_Inf_LYSATErun@meta.data <- CL003_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL003_InfLys@meta.data <- CL003_InfLys@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL007_Inf_LYSATErun@meta.data <- CL007_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL007_InfLys@meta.data <- CL007_InfLys@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL011_Inf_LYSATErun@meta.data <- CL011_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL011_InfLys@meta.data <- CL011_InfLys@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL012_Inf_LYSATErun@meta.data <- CL012_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL012_InfLys@meta.data <- CL012_InfLys@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL013_Inf_LYSATErun@meta.data <- CL013_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL013_InfLys@meta.data <- CL013_InfLys@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL018_Inf_LYSATErun@meta.data <- CL018_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
P079_Inf_LYSATErun@meta.data <- P079_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL002_Inf@meta.data <- CL002_Inf@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
CL010_Inf@meta.data <- CL010_Inf@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
Lkpk2_Inf@meta.data <- Lkpk2_Inf@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
Lkpk3_Inf@meta.data <- Lkpk3_Inf@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
Lkpk4_Inf@meta.data <- Lkpk4_Inf@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 
Lkpk5_Inf@meta.data <- Lkpk5_Inf@meta.data %>% dplyr::mutate(cloneType_v1 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 3 ~ 'Nonexpanded', clonalFrequency >= 3 ~ 'Expanded')) 

# Labeling strategy: version 2

CL003_Inf_LYSATErun@meta.data <- CL003_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL003_InfLys@meta.data <- CL003_InfLys@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL007_Inf_LYSATErun@meta.data <- CL007_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL007_InfLys@meta.data <- CL007_InfLys@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL011_Inf_LYSATErun@meta.data <- CL011_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL011_InfLys@meta.data <- CL011_InfLys@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL012_Inf_LYSATErun@meta.data <- CL012_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL012_InfLys@meta.data <- CL012_InfLys@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL013_Inf_LYSATErun@meta.data <- CL013_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL013_InfLys@meta.data <- CL013_InfLys@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL018_Inf_LYSATErun@meta.data <- CL018_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
P079_Inf_LYSATErun@meta.data <- P079_Inf_LYSATErun@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL002_Inf@meta.data <- CL002_Inf@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
CL010_Inf@meta.data <- CL010_Inf@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
Lkpk2_Inf@meta.data <- Lkpk2_Inf@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
Lkpk3_Inf@meta.data <- Lkpk3_Inf@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
Lkpk4_Inf@meta.data <- Lkpk4_Inf@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 
Lkpk5_Inf@meta.data <- Lkpk5_Inf@meta.data %>% dplyr::mutate(cloneType_v2 = case_when(clonalFrequency == 0 ~ 'Absent', clonalFrequency < 2 ~ 'Nonexpanded', clonalFrequency >= 2 ~ 'Expanded')) 

#=== 10. Merging samples ========================================================

CL_merged_all <- merge(CL002_Inf, y = list(CL010_Inf, Lkpk2_Inf, Lkpk3_Inf, Lkpk4_Inf, Lkpk5_Inf,
                                           CL003_Inf_LYSATErun, CL003_InfLys, CL007_Inf_LYSATErun, CL007_InfLys,
                                           CL011_Inf_LYSATErun, CL011_InfLys, CL012_Inf_LYSATErun, CL012_InfLys, 
                                           CL013_Inf_LYSATErun, CL013_InfLys,
                                           CL018_Inf_LYSATErun, P079_Inf_LYSATErun))

#--- Delete files that are no longer needed
rm(CL002_Inf, CL010_Inf, Lkpk2_Inf, Lkpk3_Inf, Lkpk4_Inf, Lkpk5_Inf,
   CL003_Inf_LYSATErun, CL003_InfLys, CL007_Inf_LYSATErun, CL007_InfLys,
   CL011_Inf_LYSATErun, CL011_InfLys, CL012_Inf_LYSATErun, CL012_InfLys, 
   CL013_Inf_LYSATErun, CL013_InfLys, CL018_Inf_LYSATErun, P079_Inf_LYSATErun)

rm(CL003_Inf_tcr_LYSATErun, CL003_InfLys_tcr, 
   CL007_Inf_tcr_LYSATErun, CL007_InfLys_tcr,
   CL011_Inf_tcr_LYSATErun, CL011_InfLys_tcr, 
   CL012_Inf_tcr_LYSATErun, CL012_InfLys_tcr, 
   CL013_Inf_tcr_LYSATErun, CL013_InfLys_tcr,
   CL018_Inf_tcr_LYSATErun, P079_Inf_tcr_LYSATErun,
   CL002_Inf_tcr, CL010_Inf_tcr, 
   Lkpk2_Inf_tcr, Lkpk3_Inf_tcr, 
   Lkpk4_Inf_tcr, Lkpk5_Inf_tcr)

#=== 11. Variance stabilization =================================================

CL_merged_all[["RNA"]] <- as(CL_merged_all[["RNA"]], Class = "Assay5")
CL_merged_all <- SCTransform(CL_merged_all, vars.to.regress = c('percent.mt', 'percent.rb', 'percent.tcr','S.Score', 'G2M.Score'), method = 'glmGamPoi', vst.flavor = "v2", verbose = T, variable.features.n = 2000)
CL_merged_all <- quietTCRgenes(CL_merged_all)
CL_merged_all <- Seurat::RunPCA(CL_merged_all, verbose = TRUE)

################################################################################
# II. Sample integration/Batch correction ######################################
################################################################################

# Perform integration/batch correction using Harmony
CL_harmony_all <- IntegrateLayers(object = CL_merged_all, method = HarmonyIntegration, normalization.method = "SCT", new.reduction = "harmony", dims = 1:20, verbose = T)

#--- Delete files that are no longer needed
rm(CL_merged_all)

################################################################################
# III. Postintegration analysis ################################################
################################################################################

#=== 1. Cell clustering ====

#--- Perform cell clustering ----

CL_harmony_all[["RNA"]] <- JoinLayers(CL_harmony_all[["RNA"]])
CL_harmony_all <- Seurat::RunUMAP(CL_harmony_all, reduction = "harmony", dims = 1:30, n.epochs = NULL, n.neighbors = 50)
CL_harmony_all <- Seurat::FindNeighbors(CL_harmony_all, reduction = "harmony", dims = 1:30)
CL_harmony_all <- Seurat::FindClusters(CL_harmony_all, resolution = 0.3, algorithm = 2)

#--- Visualize cell clusters ----

umap_general <- DimPlot_scCustom(CL_harmony_all, group.by = 'seurat_clusters', reduction = "umap", label = F, label.box = F, label.size = 5, repel = TRUE, pt.size = 0.2, raster = F, order = TRUE)+ 
  ggtitle(NULL)+
  guides(color = guide_legend(override.aes = list(size=7), ncol=1))+ 
  theme(legend.title = element_text(size = 25,face="bold"), 
        legend.text = element_text(size = 25, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=0, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 0, color="black"),
        axis.ticks.length=unit(0, "cm"))

# Save Supplemental Figure 4A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
umap_general
dev.off()

#--- Quantify the number of cells per cluster per condition -----

Idents(CL_harmony_all) <- CL_harmony_all$seurat_clusters
cellProps <- propeller(clusters = CL_harmony_all$seurat_clusters, sample = CL_harmony_all$Sample, group = CL_harmony_all$Condition)
cellProps$ident <- rownames(cellProps)
CL_harmony_all$Condition <- factor(CL_harmony_all$Condition, levels = c("nonLTBI:Inf", "LTBI:Inf", "LTBI:InfLys"))
props <- getTransformedProps(CL_harmony_all$Condition, CL_harmony_all$seurat_clusters, transform=NULL)
df_spekle <- as.data.frame(props$Counts)
cell.dyn <- ggplot(df_spekle, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_colour_manual(values = rainbow(3)) + 
  coord_cartesian(ylim = c(0, 42000)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 14, angle = 0, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 14, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save data for Supplemental Figure 4B

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
write_excel_csv(df_spekle,"Spekle_byClusters.csv", col_names = TRUE)
saveRDS(df_spekle, "Spekle_byClusters.rds")

#=== 2. Quality control: removal of monocytic contamination ====

#--- Visualize gene expression density of monocytic markers within the dataset ----

macrof <- Nebulosa::plot_density(CL_harmony_all, c('ITGAX', 'MRC1'), joint = TRUE, combine = FALSE)

# Save Supplemental Figure 4C

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
macrof
dev.off()

#--- Visualize gene expression levels of monocytic markers within the dataset ----

macrof_genes <- c('ITGAX', 'MRC1')
dotPlot_macrof <- SCpubr::do_DotPlot(sample = CL_harmony_all, features = macrof_genes, cluster = TRUE, dot.scale = 11, flip = FALSE)

# Save Supplemental Figure 4D

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*100, res = 300, pointsize = 1)     
dotPlot_macrof
dev.off()

#--- Find differentially expressed conserved markers for all cell clusters----

CL_harmony_copy <- CL_harmony_all
CL_harmony_copy$Condition[which(CL_harmony_copy$Condition == "nonLTBI:Inf")] <- "nonLTBI_Inf"
CL_harmony_copy$Condition[which(CL_harmony_copy$Condition == "LTBI:Inf")] <- "LTBI_Inf"
CL_harmony_copy$Condition[which(CL_harmony_copy$Condition == "LTBI:InfLys")] <- "LTBI_InfLys"
DefaultAssay(CL_harmony_copy) <- "RNA"

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(CL_harmony_copy,
                       ident.1 = cluster,
                       grouping.var = "Condition",
                       only.pos = TRUE,
                       min.pct = 0.25, 
                       logfc.threshold = 0.25) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers_v4 <- map_dfr(0:15, get_conserved)

#--- Delete files that are no longer needed

rm(CL_harmony_copy)

# Save the data

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
write_excel_csv(conserved_markers_v4,"conserved_markers_19clusters.csv", col_names = TRUE)
saveRDS(conserved_markers_v4, "conserved_markers_19clusters.rds")

# Extract top 5 markers per cluster
topn <- conserved_markers_v4 %>% 
  mutate(avg_fc = (nonLTBI_Inf_avg_log2FC + LTBI_Inf_avg_log2FC + LTBI_InfLys_avg_log2FC)/3) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 100, wt = avg_fc)

top_markers <- Extract_Top_Markers(marker_dataframe = topn, num_genes = 5, named_vector = FALSE,
                                   make_unique = TRUE,  rank_by = "avg_fc", group_by = 'cluster_id')

view(top_markers)

Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters
Idents(object = CL_harmony_all)

Unique_heatmap <- DoHeatmap(CL_harmony_all, features = top_markers)

# Save Supplemental Figure 4E

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
Unique_heatmap
dev.off()

#--- Remove cell clusters contaminated with monocytes ----
CL_harmony_all <- subset(CL_harmony_all, seurat_clusters %in% c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'))

#--- Visualize cell clusters ----

umap_general_subset <- DimPlot_scCustom(CL_harmony_all, group.by = 'seurat_clusters', reduction = "umap", label = F, label.box = F, label.size = 5, repel = TRUE, pt.size = 0.2, raster = F, order = TRUE)+ 
  ggtitle(NULL)+
  guides(color = guide_legend(override.aes = list(size=7), ncol=1))+ 
  theme(legend.title = element_text(size = 25,face="bold"), 
        legend.text = element_text(size = 25, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=0, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 0, color="black"),
        axis.ticks.length=unit(0, "cm"))

# Save Figure 5A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
umap_general_subset
dev.off()

#=== 3. Gene expression analysis ====

#--- Visualize gene expression density for CD4 T cell canonical markers ----

nebulosa_ifng <- Nebulosa::plot_density(CL_harmony_all, c('IFNG')) & NoLegend() & NoAxes()
nebulosa_tnf <- Nebulosa::plot_density(CL_harmony_all, c('TNF')) & NoLegend() & NoAxes()
nebulosa_il2 <- Nebulosa::plot_density(CL_harmony_all, c('IL2')) & NoLegend() & NoAxes()
nebulosa_csf2 <- Nebulosa::plot_density(CL_harmony_all, c('CSF2')) & NoLegend() & NoAxes()
nebulosa_gnly <- Nebulosa::plot_density(CL_harmony_all, c('GNLY')) & NoLegend() & NoAxes()
nebulosa_gzmb <- Nebulosa::plot_density(CL_harmony_all, c('GZMB')) & NoLegend() & NoAxes()
nebulosa_il17a <- Nebulosa::plot_density(CL_harmony_all, c('IL17A')) & NoLegend() & NoAxes()
nebulosa_il17f <- Nebulosa::plot_density(CL_harmony_all, c('IL17F')) & NoLegend() & NoAxes()
nebulosa_il4 <- Nebulosa::plot_density(CL_harmony_all, c('IL4')) & NoLegend() & NoAxes()
nebulosa_il13 <- Nebulosa::plot_density(CL_harmony_all, c('IL13')) & NoLegend() & NoAxes()
nebulosa_tgfb1 <- Nebulosa::plot_density(CL_harmony_all, c('TGFB1')) & NoLegend() & NoAxes()
nebulosa_il10 <- Nebulosa::plot_density(CL_harmony_all, c('IL10')) & NoLegend() & NoAxes()
nebulosa_tbx21 <- Nebulosa::plot_density(CL_harmony_all, c('TBX21')) & NoLegend() & NoAxes()
nebulosa_rorc <- Nebulosa::plot_density(CL_harmony_all, c('RORC')) & NoLegend() & NoAxes()
nebulosa_gata3 <- Nebulosa::plot_density(CL_harmony_all, c('GATA3')) & NoLegend() & NoAxes()
nebulosa_foxp3 <- Nebulosa::plot_density(CL_harmony_all, c('FOXP3')) & NoLegend() & NoAxes()

# arranging in a grid

nebulosa_canon <- grid.arrange(nebulosa_ifng,
             nebulosa_tnf,
             nebulosa_il2,
             nebulosa_csf2,
             nebulosa_gnly,
             nebulosa_gzmb,
             nebulosa_il17a,
             nebulosa_il17f,
             nebulosa_il4,
             nebulosa_il13,
             nebulosa_il10,
             nebulosa_tgfb1,
             nebulosa_tbx21,
             nebulosa_rorc,
             nebulosa_gata3,
             nebulosa_foxp3,
             nrow = 4)

# Save Figure 5B

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
nebulosa_canon
dev.off()

#--- Gene expression density plot for naive/Tcm CD4T cell canonical markers ----

naive_nebulosa <- Nebulosa::plot_density(CL_harmony_all, c('CCR7', 'SELL'), joint = TRUE, combine = FALSE)

# Save Figure 5G

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
naive_nebulosa
dev.off()

#--- Find differentially expressed conserved markers for all cell clusters after QC ----

CL_harmony_copy <- CL_harmony_all
CL_harmony_copy$Condition[which(CL_harmony_copy$Condition == "nonLTBI:Inf")] <- "nonLTBI_Inf"
CL_harmony_copy$Condition[which(CL_harmony_copy$Condition == "LTBI:Inf")] <- "LTBI_Inf"
CL_harmony_copy$Condition[which(CL_harmony_copy$Condition == "LTBI:InfLys")] <- "LTBI_InfLys"
DefaultAssay(CL_harmony_copy) <- "RNA"
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(CL_harmony_copy,
                       ident.1 = cluster,
                       grouping.var = "Condition",
                       only.pos = TRUE,
                       min.pct = 0.25, 
                       logfc.threshold = 0.25) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
conserved_markers_v3 <- map_dfr(0:15, get_conserved)

# Delete files that are no longer needed

rm(CL_harmony_copy)

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
write_excel_csv(conserved_markers_v3,"conserved_markers_v3.csv", col_names = TRUE)
saveRDS(conserved_markers_v3, "conserved_markers_v3.rds")

# Extract top 20 markers per cluster

topn <- conserved_markers %>% 
  mutate(avg_fc = (nonLTBI_Inf_avg_log2FC + LTBI_Inf_avg_log2FC + LTBI_InfLys_avg_log2FC)/3) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 100, wt = avg_fc)

top_markers <- Extract_Top_Markers(marker_dataframe = topn, num_genes = 5, named_vector = FALSE,
                                   make_unique = TRUE,  rank_by = "avg_fc", group_by = 'cluster_id')

view(top_markers)

Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters
hierarchical_cl <- Clustered_DotPlot(seurat_object = CL_harmony_all, features = top_markers, k = 0, flip = F)

# Save Figure 5C

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*1000, res = 300, pointsize = 1)     
hierarchical_cl
dev.off()

#=== 4. Assessment of clonal expansion ====

#--- Visualize clonally expanded vs nonexpanded CD4 T cells ----

CL_harmony_all$Condition <- factor(CL_harmony_all$Condition, levels = c("nonLTBI:Inf", "LTBI:Inf", "LTBI:InfLys"))
clone.size <- DimPlot(CL_harmony_all, group.by = 'cloneType_v2', split.by = 'Condition', label = F, repel = TRUE, order = TRUE, pt.size = 0.2, raster = FALSE, ncol = 3) + ggtitle(NULL) +
  ggtitle(NULL)+
  guides(color = guide_legend(override.aes = list(size=7), ncol=1))+ 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=0, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 0, color="black"),
        axis.ticks.length=unit(0, "cm"),
        strip.text = element_blank())

# Save Figure 5D

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*1000, res = 300, pointsize = 1)     
clone.size
dev.off()

#--- Quantify percentages of expanded/nonexpanded of the total number of cells in each cluster ----

harmony_subset <- subset(CL_harmony_all, subset = orig.ident %in% 'LTBI')
props <- getTransformedProps(harmony_subset$cloneType_v2, harmony_subset$seurat_clusters, transform=NULL)
df_spekle_counts <- as.data.frame(props$Counts)
# Delete files that are no longer needed
rm(harmony_subset)

# Save data for Figure 5E

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
write_excel_csv(df_spekle_counts,"ExpVSnExp.clone.counts.csv", col_names = TRUE)
saveRDS(df_spekle_counts, "ExpVSnExp.clone.counts.rds")

#--- Quantify percentages of expanded/nonexpanded of the total number of TCR reads in each cluster ----

harmony_subset <- subset(CL_harmony_all, subset = orig.ident %in% 'LTBI')
props <- getTransformedProps(harmony_subset$cloneType_v2, harmony_subset$seurat_clusters, transform=NULL)
df_spekle_props <- as.data.frame(props$Proportions)
# Delete files that are no longer needed
rm(harmony_subset)

# Save data for Figure 5F

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
write_excel_csv(df_spekle_props,"ExpVSnExp.clone.props.csv", col_names = TRUE)
saveRDS(df_spekle_props, "ExpVSnExp.clone.props.rds")

#=== 5. Single-cell trajectory analysis ====

cds <- as.cell_data_set(CL_harmony_all, group.by = 'seurat_clusters')

# plot
colData(cds)

cds <- cluster_cells(cds, resolution = 1e-5)

# Trajectory analysis
cds <- learn_graph(cds, use_partition = TRUE)

# plot
plot_cells(cds, color_cells_by = 'seurat_clusters', label_groups_by_cluster = FALSE, group_label_size = 5) + theme(legend.position = 'right')

# Create gene annotation file
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds)

# Order the cells by pseudotime
cds <- order_cells(cds)

# Plot cells in pseudotime
plot_cells(cds, color_cells_by = 'pseudotime')

# Save Figure 5H

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cds
dev.off()

#=== 6. Mapping CDR3 motifs on the UMAP  =====

#--- Mapping the CDR3 motifs known to be Mtb specific ----

CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(Known_ags = case_when(CTaa == 'CAGGWGGNNRKLIW_CASSLRSRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'CAVEGSQGNLIF_CASSLRSRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'CAVSLNDYKLSF_CASSLRSRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'CAVTLNDYKLSF_CASSLRSRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'CAVTPNDYKLSF_CASSLRSRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'CAVTRSNYKLTF_CASSLRTRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'NA_CASSLRSRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'NA_CASSLRTRSYEQYF' ~ 'SLR%RSYE (mIHF)',
                                                                                             CTaa == 'CALMNSGGYQKVTF_CASSPGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALSYSGGYQKVTF_CASSPGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALTYSGGYQKVTF_CASSPGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'NA_CASSPGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'NA_CASSPGTESNSPLHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALILSGGYQKVTF_CASSFGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALLNSGGYQKVTF_CASSTGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALSASGGYQKVTF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALSNSGGYQKVTF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALSTSGGYQKVTF_CASSVGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALTNSGGYQKVTF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CASPSNSGYALNF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'NA_CASSFGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'NA_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'NA_CASSSGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALLNSGGYQKVTF_CASSLGTESNSPLHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALSASGGYQKVTF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALSNSGGYQKVTF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CALTNSGGYQKVTF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CASPSNSGYALNF_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'NA_CASSLGTESNQPQHF' ~ 'TESN (EspA)',
                                                                                             CTaa == 'CAVRDPENTDKLIF_CSARSSGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPGNTDKLIF_CSARASGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPGNTDKLIF_CSARTSGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPLNTDKLIF_CSARSSGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPLNTDKLIF_CSARTSGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPYNTDKLIF_CSARTSGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'NA_CSARSSGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPGNTDKLIF_CSARASGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPINTDKLIF_CSARAGGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVRDPLNTDKLIF_CSARAGGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAVSRAGAGSYQLTF_CSARAGGGEAKNIQYF' ~ 'GEAK (CFP10)',
                                                                                             CTaa == 'CAAAGSSNTGKLIF_CASSRSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CAACGSSNTGKLIF_CASSSSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CAASGSSNTGKLIF_CASSSSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CAASGSSNTGKLIF_CASSTSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CAAVGSSNTGKLIF_CASSRSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CAASPGDSGGYNKLIF_CASSTSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CAASSPGGQKLLF_CASSKSGTKYNEQFF' ~ 'S%SGTKYNE (EccE3)',
                                                                                             CTaa == 'CATAKTGANNLFF_CASSSPGQGGLNYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATARTGANNLFF_CASSSPGQGGNNYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATARTGANNLFF_CASSSPGQGGSNYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATARTGANNLFF_CASSSPGQGGANYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATPNSGNTPLVF_CASSSPGQGGANYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATASTGANNLFF_CASSSPGQGGVNYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATPNSGNTPLVF_CASSSPGQGGANYGYTF' ~ 'SSPGQGG%NYG (DacB1)',
                                                                                             CTaa == 'CATSRTGANNLFF_CASSSPGQGGANYGYTF' ~ 'SSPGQGG%NYG (DacB1)'))

# Visualizing the data

Known_ags <- DimPlot(CL_harmony_all, split.by = 'orig.ident', group.by = 'Known_ags', reduction = "umap", label = F, repel = TRUE, pt.size = 3, raster = F, order = TRUE, na.value = 'lightgray') + 
  ggtitle('Known to be M.tb specific') + labs(color = "Motifs:") +
  guides(color = guide_legend(override.aes = list(size=5), ncol=1)) + 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_blank())

# Save Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
Known_ags
dev.off()

#--- Quantify GLIPH2 motifs per cluster ----

props_Known_ags <- getTransformedProps(CL_harmony_all$Known_ags, CL_harmony_all$seurat_clusters, transform=NULL)
df_known <- as.data.frame(props_Known_ags$Counts)
bars_df_knowng <- ggplot(df_known, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") + 
  coord_cartesian(ylim = c(0, 20)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 12, angle = 45, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 12, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save Figure 6B

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
bars_df_knowng
dev.off()

# Save data for Figure 6I, Figure 7A, and Figure 7B, Supplemental Figure 5E, Supplemental Figure 5G, and Supplemental Figure 5H

write_excel_csv(df_known,"df_known.csv", col_names = TRUE)

#--- Mapping consistent CDR3 motifs likely to be Mtb specific ----

CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(Resp_to_inf_mf_strong = case_when(CTaa == 'CAASGSARQLTF_CASSFVDSSYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSFVDSSYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSFVDSTYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CAPSGSARQLTF_CASSFVDSSYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CAPSGSARQLTF_CASSFVDSTYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CAVSGSARQLTF_CASSFVDSSYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CVXETSGSRLTF_CASSFVDSSYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'NA_CASSFVDSSYEQYF' ~ 'SFVDS%YE',
                                                                                                         CTaa == 'CAVYTSGTYKYIF_CSANMPEAFF' ~ '%MPE',
                                                                                                         CTaa == 'CAVYTSGTYKYIF_CSASMPEAFF' ~ '%MPE',
                                                                                                         CTaa == 'CIVRDGSSNTGKLIF_CSANMPEAFF' ~ '%MPE',
                                                                                                         CTaa == 'NA_CSASMPEAFF' ~ '%MPE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSLVAGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSLVGGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSLVQGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSLVQGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'CALSGSARQLTF_CASSLVSGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'CAVSGSARQLTF_CASSLVEGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'NA_CASSLVEGPYEQYF' ~ 'SLV%GPYE',
                                                                                                         CTaa == 'CAGHKAGGTSYGKLTF_CASSQGTGGKYEQYF' ~ 'S%GTGGKYE',
                                                                                                         CTaa == 'CAMSSAGGTSYGKLTF_CASSHGTGGKYEQYF' ~ 'S%GTGGKYE',
                                                                                                         CTaa == 'CAMSSSGGTSYGKLTF_CASSRGTGGKYEQYF' ~ 'S%GTGGKYE',
                                                                                                         CTaa == 'CAVSIAGGTSYGKLTF_CASSLGTGGKYEQYF' ~ 'S%GTGGKYE',
                                                                                                         CTaa == 'NA_CASSMGTGGKYEQYV' ~ 'S%GTGGKYE',
                                                                                                         CTaa == 'CAASGSARQLTF_CASSLVTSGTYEQYF' ~ 'S%VTSGTYE',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSFVTSGTYEQYF' ~ 'S%VTSGTYE',
                                                                                                         CTaa == 'CALSGSARQLTF_CASSLVTSGTYEQYF' ~ 'S%VTSGTYE',
                                                                                                         CTaa == 'CAVTGGGSQGNLIF_CASSESGGSNQPQHF' ~ 'SESGG%NQP',
                                                                                                         CTaa == 'CAVTNGGSQGNLIF_CASSESGGSNQPQHF' ~ 'SESGG%NQP',
                                                                                                         CTaa == 'CAVTPGGSQGNLIF_CASSESGGSNQPQHF' ~ 'SESGG%NQP',
                                                                                                         CTaa == 'CAFMKHHRRRGGTSYGKLTF_CASSESGGSNQPQHF' ~ 'SESGG%NQP',
                                                                                                         CTaa == 'CAVHPSGGYNKLIF_CASSESGGTNQPQHF' ~ 'SESGG%NQP',
                                                                                                         CTaa == 'CAVRSNDYKLSF_CASSESGGSNQPQHF' ~ 'SESGG%NQP',
                                                                                                         CTaa == 'CAGSGSARQLTF_CASSFADSNQPQHF' ~ 'S%ADSNQP',
                                                                                                         CTaa == 'CAISGSARQLTF_CASSFADSNQPQHF' ~ 'S%ADSNQP',
                                                                                                         CTaa == 'CAGAGNNRKLIW_CASSFADSNQPQHF' ~ 'S%ADSNQP',
                                                                                                         CTaa == 'CAGLYYKLSF_CASSLADSNQPQHF' ~ 'S%ADSNQP',
                                                                                                         CTaa == 'NA_CASSFADSNQPQHF' ~ 'S%ADSNQP'))

# Visualizing the data

Resp_to_inf_mf_strong <- DimPlot(CL_harmony_all, split.by = 'orig.ident', group.by = 'Resp_to_inf_mf_strong', reduction = "umap", label = F, repel = TRUE, pt.size = 3, raster = F, order = TRUE, na.value = 'lightgray') + 
  ggtitle('Likely to be M.tb specific | Strong CDR3alpha') + labs(color = "Motifs:") +
  guides(color = guide_legend(override.aes = list(size=5), ncol=1)) + 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_blank())

# Save Figure 6C

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
Resp_to_inf_mf_strong
dev.off()

#--- Quantify GLIPH2 motifs per cluster ----

props_Resp_to_inf_mf_strong <- getTransformedProps(CL_harmony_all$Resp_to_inf_mf_strong, CL_harmony_all$seurat_clusters, transform=NULL)
df_likely_strong <- as.data.frame(props_Resp_to_inf_mf_strong$Counts)
bars_df_strong <- ggplot(df_likely_strong, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") + 
  coord_cartesian(ylim = c(0, 17)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 12, angle = 45, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 12, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save Figure 6D

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
bars_df_strong
dev.off()

# Save data for Figure 6I, Figure 7A, and Figure 7B, Supplemental Figure 5E, Supplemental Figure 5G, and Supplemental Figure 5H

write_excel_csv(df_likely_strong,"df_likely_strong.csv", col_names = TRUE)

#--- Mapping inconsistent CDR3 motifs likely to be Mtb specific ----
CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(Resp_to_inf_mf_weak = case_when(CTaa == 'CAASQSYGGATNKLIF_CASSLWEGEQYF' ~ 'SLW%GE',
                                                                                                       CTaa == 'CAASRETSGSRLTF_CASSLWRGEQYF' ~ 'SLW%GE',
                                                                                                       CTaa == 'CAVKDSGGYQKVTF_CASSLWQGEQFF' ~ 'SLW%GE',
                                                                                                       CTaa == 'CAPIGSSNTGKLIF_CASSQDPQRAEQYF' ~ 'QDPQ',
                                                                                                       CTaa == 'CAVREPTGGFKTIF_CASSQDPQGGSNQPQHF' ~ 'QDPQ',
                                                                                                       CTaa == 'CIVRLHNYGQNFVF_CASSQDPQGAMNTEAFF' ~ 'QDPQ',
                                                                                                       CTaa == 'CVVKRDDKIIF_CASSQDPQGIGKLFF' ~ 'QDPQ',
                                                                                                       CTaa == 'NA_CASSQDPQGDTQYF' ~ 'QDPQ',
                                                                                                       CTaa == 'CAAAPITQGGSEKLVF_CSVSVTNTEAFF' ~ 'S%TNTE',
                                                                                                       CTaa == 'CIVRSQTGTGGFKTIF_CSVSVTNTEAFF' ~ 'S%TNTE',
                                                                                                       CTaa == 'CAFMKVNRDDKIIF_CASSTTNTEAFF' ~ 'S%TNTE',
                                                                                                       CTaa == 'CALRPSNTGKLIF_CSVSTTNTEAFF' ~ 'S%TNTE',
                                                                                                       CTaa == 'CAPSRGNQGGKLIF_CATSRDDQPQHF' ~ 'S%DDQP',
                                                                                                       CTaa == 'CAVGAKEYGNKLVF_CATSQDDQPQHF' ~ 'S%DDQP',
                                                                                                       CTaa == 'NA_CATSRDDQPQHF' ~ 'S%DDQP',
                                                                                                       CTaa == 'CALSYNSGGYQKVTF_CASRGDGYGYTF' ~ 'RGD%YG',
                                                                                                       CTaa == 'CAVKEVGGSEKLVF_CATRGDNYGYTF' ~ 'RGD%YG',
                                                                                                       CTaa == 'CAVVSTGGFKTIF_CASRGDNYGYTF' ~ 'RGD%YG',
                                                                                                       CTaa == 'CAATPDRGSTLGRLYF_CASREGGNTIYF' ~ 'REGG%T',
                                                                                                       CTaa == 'CALTLQTSGSRLTF_CAWREGGETQYF' ~ 'REGG%T',
                                                                                                       CTaa == 'CATDHARLMF_CSAREGGETQYF' ~ 'REGG%T',
                                                                                                       CTaa == 'NA_CAWREGGETQYF' ~ 'REGG%T'))

# Visualizing the data

Resp_to_inf_mf_weak <- DimPlot(CL_harmony_all, split.by = 'orig.ident', group.by = 'Resp_to_inf_mf_weak', reduction = "umap", label = F, repel = TRUE, pt.size = 3, raster = F, order = TRUE, na.value = 'lightgray') + 
  ggtitle('Likely to be M.tb specific | Weak CDR3alpha') + labs(color = "Motifs:") +
  guides(color = guide_legend(override.aes = list(size=5), ncol=1)) + 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_blank())

# Save Supplemental Figure 5A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
Resp_to_inf_mf_weak
dev.off()

#--- Quantify GLIPH2 motifs per cluster ----

props_Resp_to_inf_mf_weak <- getTransformedProps(CL_harmony_all$Resp_to_inf_mf_weak, CL_harmony_all$seurat_clusters, transform=NULL)
df_likely_weak <- as.data.frame(props_Resp_to_inf_mf_weak$Counts)
bars_df_weak <- ggplot(df_likely_weak, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  coord_cartesian(ylim = c(0, 30)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 12, angle = 45, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 12, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save Supplemental Figure 5B

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
bars_df_weak
dev.off()

# Save data for Supplemental Figure 5E, Supplemental Figure 5G, and Supplemental Figure 5H

write_excel_csv(df_likely_weak,"df_likely_weak.csv", col_names = TRUE)

#--- Mapping CDR3 motifs found in response to MTB300/Lysate stimulation ----

CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(Resp_to_stim = case_when(CTaa == 'CAASARYGGATNKLIF_CASSALFGLETQYF' ~ 'ALFG (MTB300)',
                                                                                                CTaa == 'CAGAPNYGGATNKLIF_CASSLALFGGYTF' ~ 'ALFG (MTB300)',
                                                                                                CTaa == 'CAGPSGTSYGKLTF_CASSVALFGEGYTF' ~ 'ALFG (MTB300)',
                                                                                                CTaa == 'NA_CASSVALFGETQYF' ~ 'ALFG (MTB300)',
                                                                                                CTaa == 'NA_CASSVALFGNTIYF' ~ 'ALFG (MTB300)',
                                                                                                CTaa == 'CAVHGSSNTGKLIF_CASSAHRDQPQHF' ~ 'SA%RDQP (MTB300)',
                                                                                                CTaa == 'NA_CASSANRDQPQHF' ~ 'SA%RDQP (MTB300)',
                                                                                                CTaa == 'CAASYGNNRLAF_CASSRGQGNSPLHF' ~ 'SRGQGN%P (MTB300)',
                                                                                                CTaa == 'CAVEKDSNYQLIW_CASSRGQGNQPQHF' ~ 'SRGQGN%P (MTB300)',
                                                                                                CTaa == 'CVVSLGGGFKTIF_CASSRGQGNEPQHF' ~ 'SRGQGN%P (MTB300)',
                                                                                                CTaa == 'CAASARAGGTSYGKLTF_CSARQDSYNEQFF' ~ 'R%DSYNE (Lysate)',
                                                                                                CTaa == 'CAGLVNRDDKIIF_CASRRDSYNEQFF' ~ 'R%DSYNE (Lysate)',
                                                                                                CTaa == 'CAYRSASDRTSGSRLTF_CSARRDSYNEQFF' ~ 'R%DSYNE (Lysate)',
                                                                                                CTaa == 'CAIREVETSGSRLTF_CASSGTGTGADAFF' ~ 'SGTGTGA% (Lysate)',
                                                                                                CTaa == 'CAVRSIETSGSRLTF_CASSGTGTGAEQYF' ~ 'SGTGTGA% (Lysate)',
                                                                                                CTaa == 'CALPQGGSEKLVF_CASSLRTRETQYF' ~ 'SLRT%ET (Lysate)',
                                                                                                CTaa == 'CALSRSNYQLIW_CASSLRTKETQYF' ~ 'SLRT%ET (Lysate)',
                                                                                                CTaa == 'CAVSRSNYQLIW_CASSLRTKETQYF' ~ 'SLRT%ET (Lysate)',
                                                                                                CTaa == 'CAVTLNDYKLSF_CASSLRTRETQYF' ~ 'SLRT%ET (Lysate)',
                                                                                                CTaa == 'CAAPPETGANNLFF_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'CAVGFMEYGNKLVF_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'CAGLEKNDMRF_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'CAVGANEYGNKLVF_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'CAVGFMEYGNKLVF_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'CAVGPMEYGNKLVF_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'CAVRQGGFGNVLHC_CASSRDNSPLHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'NA_CATSRDNQPQHF' ~ 'SRDN%P (Lysate)',
                                                                                                CTaa == 'NA_CATSRDNSPLHF' ~ 'SRDN%P (Lysate)'))

# Visualizing the data

Resp_to_stim <- DimPlot(CL_harmony_all, split.by = 'orig.ident', group.by = 'Resp_to_stim', reduction = "umap", label = F, repel = TRUE, pt.size = 3, raster = F, order = TRUE, na.value = 'lightgray') + 
  ggtitle('Specific for MTB300/Lysate') + labs(color = "Motifs:") +
  guides(color = guide_legend(override.aes = list(size=5), ncol=1)) + 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_blank())

# Save Figure 6E

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
Resp_to_stim
dev.off()

#--- Quantify GLIPH2 motifs per cluster ----

props_Resp_to_stim <- getTransformedProps(CL_harmony_all$Resp_to_stim, CL_harmony_all$seurat_clusters, transform=NULL)
df_stim <- as.data.frame(props_Resp_to_stim$Counts)
bars_stim <- ggplot(df_stim, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  coord_cartesian(ylim = c(0, 5)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 12, angle = 45, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 12, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save Figure 6F

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
bars_stim
dev.off()

#--- Mapping consistent CDR3 motifs likely to be  Ctrl/Viral peptide specific ----

CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(Ctrl_ags_strong = case_when(CTaa == 'CAARGAQKLVF_CASSPTSGSTDTQYF' ~ 'SPTS%STDT (SARS/EBV)',
                                                                                                   CTaa == 'CAVRGAQKLVF_CASSPTSGSTDTQYF' ~ 'SPTS%STDT (SARS/EBV)',
                                                                                                   CTaa == 'CAVSSLTQGKLIF_CASSPTSSSTDTQYF' ~ 'SPTS%STDT (SARS/EBV)',
                                                                                                   CTaa == 'NA_CASSPTSGSTDTQYF' ~ 'SPTS%STDT (SARS/EBV)',
                                                                                                   CTaa == 'NA_CASSPTSTSTDTQYF' ~ 'SPTS%STDT (SARS/EBV)',
                                                                                                   CTaa == 'CAEDDGSARQLTF_CSARADSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAFDRGSARQLTF_CSARGDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAMGPPSGYQKVTF_CASRMDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAMRSGNARLMF_CSARTDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAVAFHNAGNMLTF_CASRIDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAVGLGGGNKLTF_CASRIDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAVLTNTNAGKSTF_CASRVDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAVSLSGYSTLTF_CASRVDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CAWRTGANNLFF_CSARSDSNQPQHF' ~ 'R%DSNQP (InfluenzaA)',
                                                                                                   CTaa == 'CALRPNSGYSTLTF_CASSLSDNEQFF' ~ 'SL%DNE (SARS)',
                                                                                                   CTaa == 'CALRPNSGYSTLTF_CASSLTDNEQFF' ~ 'SL%DNE (SARS)',
                                                                                                   CTaa == 'CALSGKSGSNYKLTF_CASSLADNEQFF' ~ 'SL%DNE (SARS)',
                                                                                                   CTaa == 'NA_CASSLTDNEQFF' ~ 'SL%DNE (SARS)',
                                                                                                   CTaa == 'CAALISGSARQLTF_CASSQDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CAAMVSGSARQLTF_CASSQDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CAFMKLFTGTASKLTF_CASSHDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CAGARYTGFQKLVF_CASSQDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CAVMDSNYQLIW_CASSYDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CAVMDSNYQLIW_CATSRDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CAVRDSNYQLIW_CATSRDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'NA_CASSQDSSGANVLTF' ~ 'S%DSSGANV (SARS)',
                                                                                                   CTaa == 'CALSPSGNQFYF_CASSLGSSYNEQFF' ~ 'SLG%SYNE (SARS)',
                                                                                                   CTaa == 'CALSPSGNQFYF_CASSLGSSYNEQFF' ~ 'SLG%SYNE (SARS)',
                                                                                                   CTaa == 'CASVSGGYNKLIF_CASSLGSSYNEQFF' ~ 'SLG%SYNE (SARS)',
                                                                                                   CTaa == 'CAVKGGGYQKVTF_CASSLGASYNEQFF' ~ 'SLG%SYNE (SARS)',
                                                                                                   CTaa == 'NA_CASSLGASYNEQFF' ~ 'SLG%SYNE (SARS)',
                                                                                                   CTaa == 'NA_CASSLGDSYNEQFF' ~ 'SLG%SYNE (SARS)',
                                                                                                   CTaa == 'CAALRQGGKLIF_CASSLGGTSTDTQYF' ~ 'S%GGTSTDT (SARS)',
                                                                                                   CTaa == 'CAMSAGGAQKLVF_CASSAGGTSTDTQYF' ~ 'S%GGTSTDT (SARS)',
                                                                                                   CTaa == 'CAMSAGGAQKLVF_CASSGGGTSTDTQYF' ~ 'S%GGTSTDT (SARS)',
                                                                                                   CTaa == 'CAVKMGDDKIIF_CASSVGGTSTDTQYF' ~ 'S%GGTSTDT (SARS)',
                                                                                                   CTaa == 'NA_CAWSVGGTSTDTQYF' ~ 'S%GGTSTDT (SARS)',
                                                                                                   CTaa == 'CAASWNNQGGKLIF_CASSRGGNYEQYF' ~ 'SRGG%YE (SARS)',
                                                                                                   CTaa == 'CAVIFYNQGGKLIF_CASSRGGNYEQYF' ~ 'SRGG%YE (SARS)',
                                                                                                   CTaa == 'CAVSPIQAGWGKLQF_CASSRGGTYEQYF' ~ 'SRGG%YE (SARS)',
                                                                                                   CTaa == 'CGAADTPGGTSYGKLTF_CASSRGGGYEQYF' ~ 'SRGG%YE (SARS)',
                                                                                                   CTaa == 'CLVGGPNTGFQKLVF_CASSRGGSYEQYF' ~ 'SRGG%YE (SARS)',
                                                                                                   CTaa == 'CAARIQTGANNLFF_CASSGGNYGYTF' ~ '%GGNYG (SARS)',
                                                                                                   CTaa == 'CAFLSGNTPLVF_CSVSGGNYGYTF' ~ '%GGNYG (SARS)',
                                                                                                   CTaa == 'CAGLSGNTPLVF_CSVSGGNYGYTF' ~ '%GGNYG (SARS)',
                                                                                                   CTaa == 'NA_CSGTGGNYGYTF' ~ '%GGNYG (SARS)',
                                                                                                   CTaa == 'CALNSGAGSYQLTF_CASNPGQGYEQYF' ~ '%PGQGYE (InfluenzaA/SARS/YFV)',
                                                                                                   CTaa == 'CALNSGAGSYQLTF_CASSPGQGYEQYF' ~ '%PGQGYE (InfluenzaA/SARS/YFV)',
                                                                                                   CTaa == 'CVVWGGFKTIF_CASSPGQGYEQYF' ~ '%PGQGYE (InfluenzaA/SARS/YFV)',
                                                                                                   CTaa == 'CAGIAGGTSYGKLTF_CASSAGTSSYNEQFF' ~ 'S%GTSSYNE (SARS/CMV)',
                                                                                                   CTaa == 'CAGISGGTSYGKLTF_CASSSGTSSYNEQFF' ~ 'S%GTSSYNE (SARS/CMV)',
                                                                                                   CTaa == 'NA_CASSDGTSSYNEQFF' ~ 'S%GTSSYNE (SARS/CMV)',
                                                                                                   CTaa == 'CAMSHHQGAQKLVF_CASSLVDSTYEQYF' ~ 'SLVDS%YE (SARS)',
                                                                                                   CTaa == 'CAPSGSARQLTF_CASSLVDSSYEQYF' ~ 'SLVDS%YE (SARS)',
                                                                                                   CTaa == 'CAVSGSARQLTF_CASSLVDSAYEQYF' ~ 'SLVDS%YE (SARS)',
                                                                                                   CTaa == 'CVVRTGAGNMLTF_CASSLVDSSYEQYF' ~ 'SLVDS%YE (SARS)',
                                                                                                   CTaa == 'NA_CASSLVDSAYEQYF' ~ 'SLVDS%YE (SARS)',
                                                                                                   CTaa == 'CAGVDSNYQLIW_CASSDTSGGADTQYF' ~ 'SD%SGGADT (SARS/InfluenzaA)',
                                                                                                   CTaa == 'CAVMDSSYKLIF_CASSDSSGGADTQYF' ~ 'SD%SGGADT (SARS/InfluenzaA)',
                                                                                                   CTaa == 'CAVRDRNYQLIW_CASSDSSGGADTQYF' ~ 'SD%SGGADT (SARS/InfluenzaA)',
                                                                                                   CTaa == 'NA_CASSDSSGGADTQYF' ~ 'SD%SGGADT (SARS/InfluenzaA)',
                                                                                                   CTaa == 'CAASYRDDKIIF_CASSVAGTSTDTQYF' ~ 'SV%GTSTDT (SARS)',
                                                                                                   CTaa == 'CAVKMGDDKIIF_CASSVGGTSTDTQYF' ~ 'SV%GTSTDT (SARS)',
                                                                                                   CTaa == 'NA_CAWSVGGTSTDTQYF' ~ 'SV%GTSTDT (SARS)',
                                                                                                   CTaa == 'CAIRSGGYNKLIF_CASSEDGGYGYTF' ~ 'SED%GYG (SARS)',
                                                                                                   CTaa == 'CATRSGGYNKLIF_CASSEDSGYGYTF' ~ 'SED%GYG (SARS)',
                                                                                                   CTaa == 'CAAISNDYKLSF_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CAASASGGSQGNLIF_CASSLGPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CAFFSGGGADGLTF_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CAFYSGGGADGLTF_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CATPSYGGSQGNLIF_CASSLEPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CAVAPDTGRRALTF_CASSLDPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CAVRWANSKLTF_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CDYKLSF_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'CIVRAYSGNTGKLIF_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)',
                                                                                                   CTaa == 'NA_CASSLSPGNYGYTF' ~ 'SL%PGNYG (InfluenzaA)'))

# Visualizing the data

Resp_to_ctrl_strong <- DimPlot(CL_harmony_all, split.by = 'orig.ident', group.by = 'Ctrl_ags_strong', reduction = "umap", label = F, repel = TRUE, pt.size = 3, raster = F, order = TRUE, na.value = 'lightgray') + 
  ggtitle('Known to be for Ctrl/Viral peptide-specific') + labs(color = "Motifs:") +
  guides(color = guide_legend(override.aes = list(size=5), ncol=1)) + 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_blank())

# Save Figure 6G

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
Resp_to_ctrl_strong
dev.off()

#--- Quantify GLIPH2 motifs per cluster ----

props_Ctrl_ags_strong <- getTransformedProps(CL_harmony_all$Ctrl_ags_strong, CL_harmony_all$seurat_clusters, transform=NULL)
df_ctrl_strong <- as.data.frame(props_Ctrl_ags_strong$Counts)
bars_ctrl_strong <- ggplot(df_ctrl_strong, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  coord_cartesian(ylim = c(0, 30)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 12, angle = 45, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 12, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save Figure 6H

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
bars_ctrl_strong
dev.off()

# Save data for Figure 6I, Figure 7A, and Figure 7B, Supplemental Figure 5F, Supplemental Figure 5G, and Supplemental Figure 5H

write_excel_csv(df_ctrl_strong,"df_ctrl_strong.csv", col_names = TRUE)

#--- Mapping inconsistent CDR3 motifs likely to be  Ctrl/Viral peptide specific ----

CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(Ctrl_ags_weak = case_when(CTaa == 'CAEDNTGTASKLTF_CASSLAQGEQYF' ~ 'SL%QGE (InfluenzaA)',
                                                                                                 CTaa == 'CALSDYNQGGKLIF_CASSLAQGEQYF' ~ 'SL%QGE (InfluenzaA)',
                                                                                                 CTaa == 'CAMSSNSNSGYALNF_CSASLTQGEQFF' ~ 'SL%QGE (InfluenzaA)',
                                                                                                 CTaa == 'CATDGRYGNKLVF_CASSLGQGEQFF' ~ 'SL%QGE (InfluenzaA)',
                                                                                                 CTaa == 'CAFMKPWDASGSRLTF_CASSVSGNTEAFF' ~ 'SV%GNTE (SARS/Mtb)',
                                                                                                 CTaa == 'CALSSRSGAGSYQLTF_CAISVQGNTEAFF' ~ 'SV%GNTE (SARS/Mtb)',
                                                                                                 CTaa == 'CVRGPSSNTGKLIF_CASSVGGNTEAFF' ~ 'SV%GNTE (SARS/Mtb)',
                                                                                                 CTaa == 'CAASRTGGTSYGKLTF_CASSYSGDTQYF' ~ 'SYSG%T (SARS)',
                                                                                                 CTaa == 'CAATLSGNTPLVF_CASSYSGETQYF' ~ 'SYSG%T (SARS)',
                                                                                                 CTaa == 'CALSDLISGGGADGLTF_CASSYSGETQYF' ~ 'SYSG%T (SARS)',
                                                                                                 CTaa == 'CAVRDDSNYQLIW_CASSYSGETQYF' ~ 'SYSG%T (SARS)',
                                                                                                 CTaa == 'NA_CASSYSGETQYF' ~ 'SYSG%T (SARS)',
                                                                                                 CTaa == 'CAVARSGQPAGGTSYGKLTF_CASSYGPDTQYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'CAVKVNTNAGKSTF_CASSLGPDTIYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'CAVNLGGSTLGRLYF_CASSLGPDTQYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'CAVRPGGATNKLIF_CASSYGPDTQYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'CVVNAAGTASKLTF_CASSFGPDTQYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'NA_CASSLGPDTIYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'NA_CASSWGPDTQYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'NA_CASSYGPDTQYF' ~ 'S%GPDT (SARS)',
                                                                                                 CTaa == 'CALVGARLMF_CASSLEGGSYNEQFF' ~ 'SLEG%SYNE (SARS/HBV)',
                                                                                                 CTaa == 'CAVTASQGSARQLTF_CASSLEGSSYNEQFF' ~ 'SLEG%SYNE (SARS/HBV)',
                                                                                                 CTaa == 'CAVVGAPLVF_CASSLEGGSYNEQFF' ~ 'SLEG%SYNE (SARS/HBV)',
                                                                                                 CTaa == 'NA_CASSLEGSSYNEQFF' ~ 'SLEG%SYNE (SARS/HBV)',
                                                                                                 CTaa == 'NA_CASSLEGSSYNEQSF' ~ 'SLEG%SYNE (SARS/HBV)',
                                                                                                 CTaa == 'CAVAAGVGGSQGNLIF_CASSELENTEAFF' ~ 'S%LENTE (SARS)',
                                                                                                 CTaa == 'CAVGAGTGGYNKLIF_CASSELENTEAFF' ~ 'S%LENTE (SARS)',
                                                                                                 CTaa == 'CAVGVREGANNLFF_CASSELENTEAFF' ~ 'S%LENTE (SARS)',
                                                                                                 CTaa == 'CAVRGTNAGKSTF_CASSELENTEAFF' ~ 'S%LENTE (SARS)',
                                                                                                 CTaa == 'NA_CASSQLENTEAFF' ~ 'S%LENTE (SARS)',
                                                                                                 CTaa == 'CAAIVRGSNYKLTF_CASSNSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAALSGSARQLTF_CASSKSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAAPLYNQGGKLIF_CAISNSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAFMKHVLTSGTYKYIF_CAISDSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAGHYGGSQGNLIF_CASSYSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CALGMDSSYKLIF_CASSHSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAMSPDYGQNFVF_CASSGSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAVEEGXDTGRRALTF_CASSRSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CAVSYASGGSYIPTF_CASSRSNQPQHF' ~ 'S%SNQP (SARS/InfluenzaA)',
                                                                                                 CTaa == 'CALSPTGRRALTF_CASSQDGNYGYTF' ~ 'SQ%GNYG (CMV)',
                                                                                                 CTaa == 'CVVIRDARLMF_CASSQEGNYGYTF' ~ 'SQ%GNYG (CMV)',
                                                                                                 CTaa == 'NA_CASSQEGNYGYTF' ~ 'SQ%GNYG (CMV)',
                                                                                                 CTaa == 'CAVEPGAGSYQLTF_CSASRTGTYEQYF' ~ 'SRTG%YE (InfluenzaA)',
                                                                                                 CTaa == 'CAVKVDSSYKLIF_CASSRTGTYEQYF' ~ 'SRTG%YE (InfluenzaA)',
                                                                                                 CTaa == 'CILTSYDKVIF_CSASRTGVYEQYF' ~ 'SRTG%YE (InfluenzaA)',
                                                                                                 CTaa == 'CAVAAGVGGSQGNLIF_CASSELENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CAVADSSYKLIF_CASSEAENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CAVGAGTGGYNKLIF_CASSELENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CAVGAREYGNKLVF_CASSEIENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CAVGVREGANNLFF_CASSELENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CAVRGTNAGKSTF_CASSELENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CVVNINDPGGGADGLTF_CASSEVENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'NA_CASSEIENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'NA_CASSEVENTEAFF' ~ 'SE%ENTE (SARS)',
                                                                                                 CTaa == 'CAFMKHVQGAQKLVF_CSATSYNEQFF' ~ 'T%YNE (SARS/CMV)',
                                                                                                 CTaa == 'CAGDWDSSYKLIF_CASTDYNEQFF' ~ 'T%YNE (SARS/CMV)',
                                                                                                 CTaa == 'CATDAVPFLGAQKLVF_CASTGYNEQFF' ~ 'T%YNE (SARS/CMV)',
                                                                                                 CTaa == 'NA_CSATSYNEQFF' ~ 'T%YNE (SARS/CMV)',
                                                                                                 CTaa == 'CAESLDQAGTALIF_CASSFSETQYF' ~ 'SF%ET (SARS/CMV)',
                                                                                                 CTaa == 'CAVGANTGNQFYF_CASSFHETQYF' ~ 'SF%ET (SARS/CMV)',
                                                                                                 CTaa == 'CIVKTNSGGSNYKLTF_CASSFEETQYF' ~ 'SF%ET (SARS/CMV)',
                                                                                                 CTaa == 'CVVSDQAGTASKLTF_CASSFGETQYF' ~ 'SF%ET (SARS/CMV)',
                                                                                                 CTaa == 'NA_CASSFEETQYF' ~ 'SF%ET (SARS/CMV)'))

# Visualizing the data

Resp_to_ctrl_weak <- DimPlot(CL_harmony_all, split.by = 'orig.ident', group.by = 'Ctrl_ags_weak', reduction = "umap", label = F, repel = TRUE, pt.size = 3, raster = F, order = TRUE, na.value = 'lightgray') + 
  ggtitle('Known to be for Ctrl/Viral peptide-specific') + labs(color = "Motifs:") +
  guides(color = guide_legend(override.aes = list(size=5), ncol=1)) + 
  theme(legend.title = element_text(size = 14,face="bold"), 
        legend.text = element_text(size = 14, face="bold"),
        axis.title = element_text(size = 14,face= "bold"),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.line.x = element_line(linewidth=1, color="black"),
        axis.line.y = element_line(linewidth=1, color="black"),
        axis.ticks = element_line(linewidth = 1, color="black"),
        axis.ticks.length=unit(.25, "cm"),
        strip.text = element_blank())

# Save Supplemental Figure 5C

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
Resp_to_ctrl_weak
dev.off()

#--- Quantify GLIPH2 motifs per cluster ----

props_Ctrl_ags_weak <- getTransformedProps(CL_harmony_all$Ctrl_ags_weak, CL_harmony_all$seurat_clusters, transform=NULL)
df_ctrl_weak <- as.data.frame(props_Ctrl_ags_weak$Counts)
bars_ctrl_weak <- ggplot(df_ctrl_weak, aes(fill=clusters, y=Freq, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  coord_cartesian(ylim = c(0, 20)) +
  theme(axis.text.x = element_text(face = "bold",color = "black",size = 12, angle = 45, , hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(face = "bold", color = "black",size = 12, angle = 0),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Save Supplemental Figure 5D

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1000, height = 5*500, res = 300, pointsize = 1)     
df_ctrl_weak
dev.off()

# Save data for Supplemental Figure 5F, Supplemental Figure 5G, and Supplemental Figure 5H

write_excel_csv(df_ctrl_weak,"df_ctrl_weak.csv", col_names = TRUE)

#=== 7. Differential gene expression analysis ====

#--- Perform DEG analysis for each cluster of the integrated dataset ----

Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters
CL_h.markers <- FindAllMarkers(CL_harmony_all, only.pos = F, min.pct = 0.25, logfc.threshold = 0, test.use = 'MAST', recorrect_umi=FALSE)

#--- Add percentage difference to DEG results
CL_h.markers <- Add_Pct_Diff(CL_h.markers, pct.1_name = "pct.1", pct.2_name = "pct.2", overwrite = FALSE)

# Save data for Figure 7C, Figure 7D, Figure 7E, Figure 7F, Supplemental Figure 6A, and Supplemental Figure 6B

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
write_excel_csv(CL_h.markers,"FindAllMarkers_16_clusters.csv", col_names = TRUE)
saveRDS(CL_h.markers, "FindAllMarkers_16_clusters.rds")

#--- Visualization and hierarchichal clustering

CL_h.markers %>% dplyr::filter(!str_detect(row.names(CL_h.markers), "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA | ^MT- | ^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*")) %>% 
  group_by(cluster)

Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters
CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(pathogen_spec = case_when(seurat_clusters == '6' ~ 'nonspecific',
                                                                                                 seurat_clusters == '12' ~ 'nonspecific',
                                                                                                 seurat_clusters == '14' ~ 'nonspecific',
                                                                                                 seurat_clusters == '4' ~ 'Mtb-specific',
                                                                                                 seurat_clusters == '11' ~ 'Mtb-specific',
                                                                                                 seurat_clusters == '13' ~ 'Mtb-specific',
                                                                                                 seurat_clusters == '15' ~ 'Mtb-specific'))

CL_harmony_all@meta.data <- CL_harmony_all@meta.data %>% dplyr::mutate(cluster_order = case_when(seurat_clusters == '6' ~ '1',
                                                                                                 seurat_clusters == '12' ~ '2',
                                                                                                 seurat_clusters == '14' ~ '3',
                                                                                                 seurat_clusters == '4' ~ '4',
                                                                                                 seurat_clusters == '11' ~ '5',
                                                                                                 seurat_clusters == '13' ~ '6',
                                                                                                 seurat_clusters == '15' ~ '7'))


CL_harmony_sub <- subset(CL_harmony_all, seurat_clusters %in% c('6', '12', '14', '11', '13', '15'))

DEG_markers_sub <- subset(CL_h.markers, cluster %in% c('6', '12', '14', '4', '11', '13', '15'))
top_markers_sub <- Extract_Top_Markers(marker_dataframe = DEG_markers_sub, num_genes = 10, named_vector = FALSE,
                                       make_unique = TRUE,  rank_by = "avg_log2FC")


Idents(CL_harmony_sub) <- CL_harmony_sub@meta.data$seurat_clusters
CL_dittoHeatmap_sub <- dittoHeatmap(CL_harmony_sub, top_markers_sub, annot.by = c('seurat_clusters', 'pathogen_spec'), 
                                    order.by = "cluster_order",
                                    scaled.to.max = T,
                                    complex = TRUE,
                                    use_raster = TRUE, row_km = 0, 
                                    annot.colors = c('white','white', 'white', 'white', 'palegreen3', 'white', 'orange', 
                                                     'white', 'white', 'white', 'white', 'plum3', 'mediumpurple4', 
                                                     'rosybrown2', 'royalblue4', 'brown4', 'white', 'white', 'white', 
                                                     'blue3', 'firebrick3'))


# Save Figure 7C

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*1000, res = 300, pointsize = 1)     
CL_dittoHeatmap_sub
dev.off()

#=== 8. Biological theme comparison/Reactome pathway overrepresentation analysis ====

#--- cluster 6 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c6 <- subset(CL_h.markers, cluster %in% '6')
deg_test_c6 <- subset(deg_test_c6, subset = avg_log2FC > 0)
view(deg_test_c6)
# Map the genes and convert to Entrez IDs
mapped_genes_c6 <- data.frame(GeneName = deg_test_c6$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c6$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c6$entrezID <- mapped_genes_c6$entrezID 

#
geneList_c6 <- deg_test_c6[,2]
names(geneList_c6) = as.character(deg_test_c6[,9])
view(geneList_c6)
geneList_c6 = sort(geneList_c6, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c6 <- na.omit(geneList_c6)

de_c6 <- names(geneList_c6)[abs(geneList_c6) > 1.0]

ora_c6 <- enrichPathway(gene = de_c6, pvalueCutoff = 0.05, readable=TRUE)
#ora_c6@result <- ora_c6[ora_c6@result$Count > 1]
view(ora_c6)

ora_c6 <- setReadable(ora_c6, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c6, showCategory=10)
cnetplot_c6 <- cnetplot(ora_c6, categorySize="pvalue", color.params = list(foldChange = geneList_c6, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c6
dev.off()

#--- cluster 12 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c12 <- subset(CL_h.markers, cluster %in% '12')
deg_test_c12 <- subset(deg_test_c12, subset = avg_log2FC > 0)
view(deg_test_c12)
# Map the genes and convert to Entrez IDs
mapped_genes_c12 <- data.frame(GeneName = deg_test_c12$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c12$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c12$entrezID <- mapped_genes_c12$entrezID 

#
geneList_c12 <- deg_test_c12[,2]
names(geneList_c12) = as.character(deg_test_c12[,9])
view(geneList_c12)
geneList_c12 = sort(geneList_c12, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c12 <- na.omit(geneList_c12)

de_c12 <- names(geneList_c12)[abs(geneList_c12) > 1.0]

ora_c12 <- enrichPathway(gene = de_c12, pvalueCutoff = 0.05, readable=TRUE)
view(ora_c12)

ora_c12 <- setReadable(ora_c12, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c12, showCategory=10)
cnetplot_c12 <- cnetplot(ora_c12, categorySize="pvalue", color.params = list(foldChange = geneList_c12, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c12
dev.off()

#--- cluster 14 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c14 <- subset(CL_h.markers, cluster %in% '14')
deg_test_c14 <- subset(deg_test_c14, subset = avg_log2FC > 0)
view(deg_test_c14)
# Map the genes and convert to Entrez IDs
mapped_genes_c14 <- data.frame(GeneName = deg_test_c14$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c14$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c14$entrezID <- mapped_genes_c14$entrezID 

#
geneList_c14 <- deg_test_c14[,2]
names(geneList_c14) = as.character(deg_test_c14[,9])
view(geneList_c14)
geneList_c14 = sort(geneList_c14, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c14 <- na.omit(geneList_c14)

de_c14 <- names(geneList_c14)[abs(geneList_c14) > 1.0]

ora_c14 <- enrichPathway(gene = de_c14, pvalueCutoff = 0.05, readable=TRUE)
view(ora_c14)

ora_c14 <- setReadable(ora_c14, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c14, showCategory=10)
cnetplot_c14 <- cnetplot(ora_c14, categorySize="pvalue", color.params = list(foldChange = geneList_c14, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c14
dev.off()

#--- cluster 4 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c4 <- subset(CL_h.markers, cluster %in% '4')
deg_test_c4 <- subset(deg_test_c4, subset = avg_log2FC > 0)
view(deg_test_c6)
# Map the genes and convert to Entrez IDs
mapped_genes_c4 <- data.frame(GeneName = deg_test_c4$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c4$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c4$entrezID <- mapped_genes_c4$entrezID 

#
geneList_c4 <- deg_test_c4[,2]
names(geneList_c4) = as.character(deg_test_c4[,9])
view(geneList_c4)
geneList_c4 = sort(geneList_c4, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c4 <- na.omit(geneList_c4)

de_c4 <- names(geneList_c4)[abs(geneList_c4) > 1.0]

ora_c4 <- enrichPathway(gene = de_c4, pvalueCutoff = 0.05, readable=TRUE)
#ora_c4@result <- ora_c4[ora_c4@result$Count > 1]
view(ora_c4)

ora_c4 <- setReadable(ora_c4, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c4, showCategory=10)
cnetplot_c4 <- cnetplot(ora_c4, categorySize="pvalue", color.params = list(foldChange = geneList_c4, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c4
dev.off()

#--- cluster 11 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c11 <- subset(CL_h.markers, cluster %in% '11')
deg_test_c11 <- subset(deg_test_c11, subset = avg_log2FC > 0)
view(deg_test_c11)
# Map the genes and convert to Entrez IDs
mapped_genes_c11 <- data.frame(GeneName = deg_test_c11$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c11$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c11$entrezID <- mapped_genes_c11$entrezID 

#
geneList_c11 <- deg_test_c11[,2]
names(geneList_c11) = as.character(deg_test_c11[,9])
view(geneList_c11)
geneList_c11 = sort(geneList_c11, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c11 <- na.omit(geneList_c11)

de_c11 <- names(geneList_c11)[abs(geneList_c11) > 1.0]

ora_c11 <- enrichPathway(gene = de_c11, pvalueCutoff = 0.05, readable=TRUE)
view(ora_c11)

ora_c11 <- setReadable(ora_c11, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c11, showCategory=10)
cnetplot_c11 <- cnetplot(ora_c11, categorySize="pvalue", color.params = list(foldChange = geneList_c11, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c11
dev.off()

#--- cluster 13 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c13 <- subset(CL_h.markers, cluster %in% '13')
deg_test_c13 <- subset(deg_test_c13, subset = avg_log2FC > 0)
view(deg_test_c13)
# Map the genes and convert to Entrez IDs
mapped_genes_c13 <- data.frame(GeneName = deg_test_c13$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c13$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c13$entrezID <- mapped_genes_c13$entrezID 

#
geneList_c13 <- deg_test_c13[,2]
names(geneList_c13) = as.character(deg_test_c13[,9])
view(geneList_c13)
geneList_c13 = sort(geneList_c13, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c13 <- na.omit(geneList_c13)

de_c13 <- names(geneList_c13)[abs(geneList_c13) > 1.0]

ora_c13 <- enrichPathway(gene = de_c13, pvalueCutoff = 0.05, readable=TRUE)
view(ora_c13)

ora_c13 <- setReadable(ora_c13, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c13, showCategory=10)
cnetplot_c13 <- cnetplot(ora_c13, categorySize="pvalue", color.params = list(foldChange = geneList_c13, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c13
dev.off()

#--- cluster 15 -----------------------------------------------------------------
# Load the DEG list to analyze
deg_test_c15 <- subset(CL_h.markers, cluster %in% '15')
deg_test_c15 <- subset(deg_test_c15, subset = avg_log2FC > 0)
view(deg_test_c15)
# Map the genes and convert to Entrez IDs
mapped_genes_c15 <- data.frame(GeneName = deg_test_c15$gene, entrezID = mapIds(org.Hs.eg.db, keys = deg_test_c15$gene, keytype = "SYMBOL", column="ENTREZID"))
deg_test_c15$entrezID <- mapped_genes_c15$entrezID 

#
geneList_c15 <- deg_test_c15[,2]
names(geneList_c15) = as.character(deg_test_c15[,9])
view(geneList_c15)
geneList_c15 = sort(geneList_c15, decreasing = TRUE)
#

# Ommit any NA values 
geneList_c15 <- na.omit(geneList_c15)

de_c15 <- names(geneList_c15)[abs(geneList_c15) > 1.0]

ora_c15 <- enrichPathway(gene = de_c15, pvalueCutoff = 0.05, readable=TRUE)
view(ora_c15)

ora_c15 <- setReadable(ora_c15, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
require(DOSE)
dotplot(ora_c15, showCategory=10)
cnetplot_c15 <- cnetplot(ora_c15, categorySize="pvalue", color.params = list(foldChange = geneList_c15, category = 'black', edge = F), showCategory = 7, cex.params = list(category_label = 1.2, gene_label = 1), circular = F, highlight = 'all') + scale_color_gradient2(name='Fold change', low='blue', high='firebrick')

# Save Supplemental Figure 6A

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*500, height = 5*500, res = 300, pointsize = 1)     
cnetplot_c15
dev.off()

#--- Biological theme comparison ----

de_theme <- list(de_c6, de_c12, de_c14, de_c4, de_c11, de_c13, de_c15)
str(de_theme)
names(de_theme) <- c('C6', 'C12', 'C14', 'C4', 'C11', 'C13', 'C15')
str(de_theme)

ck <- compareCluster(geneCluster = de_theme, fun = enrichPathway)
dp <- dotplot(ck, showCategory=5) + scale_y_discrete(labels=function(x) str_wrap(x, width=60))

# Save Figure 7D

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
ggsave(dp, filename="dotplot.pdf", height = 290, width = 180, units = "mm")

# Save data for Figure 7D

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
ck <- as.data.frame(ck)
write_excel_csv(ck,"BioThemeComparison.csv", col_names = TRUE)

#=== 9. Cell-cell communication analysis ====

#--- Read in the nichenetr V2 ligand-target matrices ----
lr_network <- readRDS("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/NichenetrV2_files/lr_network_human_21122021.rds")
lr_network <- lr_network %>% distinct(from, to)
ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

weighted_networks <- readRDS("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/NichenetrV2_files/weighted_networks_nsga2r_final.rds")
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from", "to"))

ligand_target_matrix <- readRDS("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/NichenetrV2_files/ligand_target_matrix_nsga2r_final.rds")

#--- Define gene sets of receiver clusters ----

deg_test_c6 <- subset(CL_h.markers, cluster %in% '6')
deg_test_c6 <- subset(deg_test_c6, subset = avg_log2FC > 0 & p_val_adj < 0.05) %>% pull('gene') %>% .[. %in% rownames(ligand_target_matrix)]

deg_test_c12 <- subset(CL_h.markers, cluster %in% '12')
deg_test_c12 <- subset(deg_test_c12, subset = avg_log2FC > 0 & p_val_adj < 0.05) %>% pull('gene') %>% .[. %in% rownames(ligand_target_matrix)]

deg_test_c14 <- subset(CL_h.markers, cluster %in% '14')
deg_test_c14 <- subset(deg_test_c14, subset = avg_log2FC > 0 & p_val_adj < 0.05) %>% pull('gene') %>% .[. %in% rownames(ligand_target_matrix)]

#--- CCC between clusters 4, 11, 13, 15 vs 12 ----
DefaultAssay(CL_harmony_all) <- 'RNA'
Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters

# a. Defining the sender and the receiver population
sender <- c('4', '11', '13', '15')
receiver <- c('12')

# b. Get expressed genes from sender and receiver
list_gex_sender <- sender %>% unique() %>% lapply(get_expressed_genes, CL_harmony_all, pct = 0.25)

gex_sender <- list_gex_sender %>% unlist() %>% unique()

list_gex_receiver <- receiver %>% unique() %>% lapply(get_expressed_genes, CL_harmony_all, pct = 0.25)

gex_receiver <- list_gex_receiver %>% unlist() %>% unique()

background_gex <- gex_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# c. Define a set of potential expressed ligands
expressed_ligands <- intersect(ligands, gex_sender)
expressed_receptors <- intersect(receptors, gex_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands

# d. Nichenetr ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = deg_test_c12, 
                                               background_expressed_genes = background_gex, 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               potential_ligands = potential_ligands)

ligand_activities <- arrange(ligand_activities, -aupr_corrected)

best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

DotPlot(CL_harmony_all, features = best_upstream_ligands, cols = 'RdYlBu') + RotatedAxis()

# e. Infer receptors and top-predicted target genes of ligands that are top ranked in the ligand-activity analysis
#--- Active target gene inference
active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = deg_test_c12, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                  ligand_target_matrix = ligand_target_matrix,
                                                                  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()

rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

p_ligand_target_network_12 <- vis_ligand_target %>% make_heatmap_ggplot('Prioritized ligands', 'Predicted target genes',color = 'purple', legend_position = 'top',legend_title = 'Regulatory potential') + 
  theme(axis.text.x = element_text(face = 'italic')) + scale_fill_gradient2(low = 'whitesmoke', high = 'purple')

p_ligand_target_network_12

# Save Supplemental Figure 6B

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_target_network_12
dev.off()

#--- Receptors of top-ranked ligands
lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from, to)

best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df <- lr_network_top_df_large %>% spread('from', 'weight', fill = 0)

lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors <- dist(lr_network_top_matrix, method = 'binary')
hclust_receptors <- hclust(dist_receptors, method = 'ward.D2')
order_receptors <- hclust_receptors$labels[hclust_receptors$order]

dist_ligands <- dist(lr_network_top_matrix %>% t(), method = 'binary')
hclust_ligands <- hclust(dist_ligands, method = 'ward.D2')
order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]

rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()

p_ligand_receptor_network_12 = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot('Ligands', 'Receptors', color = 'mediumvioletred', x_axis_position = 'top', legend_title = 'Prior interaction potential')

p_ligand_receptor_network_12

# Save Figure 7E

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_receptor_network_12
dev.off()

#--- Sender-focused approach
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands
potential_ligands_focused <- intersect(potential_ligands, gex_sender) 

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr_12 <- make_heatmap_ggplot(vis_ligand_aupr,
                                        "Prioritized ligands", "Ligand activity", 
                                        legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr_12

# Save Supplemental Figure 6B

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_aupr_12

#--- CCC between clusters 4, 11, 13, 15 vs 6 ----
DefaultAssay(CL_harmony_all) <- 'RNA'
Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters

# a. Defining the sender and the receiver population
sender <- c( '4', '11', '13', '15')
receiver <- c('6')

# b. Get expressed genes from sender and receiver
list_gex_sender <- sender %>% unique() %>% lapply(get_expressed_genes, CL_harmony_all, pct = 0.25)

gex_sender <- list_gex_sender %>% unlist() %>% unique()

list_gex_receiver <- receiver %>% unique() %>% lapply(get_expressed_genes, CL_harmony_all, pct = 0.25)

gex_receiver <- list_gex_receiver %>% unlist() %>% unique()

background_gex <- gex_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# c. Define a set of potential expressed ligands
expressed_ligands <- intersect(ligands, gex_sender)
expressed_receptors <- intersect(receptors, gex_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands

# d. Nichenetr ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = deg_test_c6, 
                                               background_expressed_genes = background_gex, 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               potential_ligands = potential_ligands)

ligand_activities <- arrange(ligand_activities, -aupr_corrected)

best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

DotPlot(CL_harmony_all, features = best_upstream_ligands, cols = 'RdYlBu') + RotatedAxis()

# e. Infer receptors and top-predicted target genes of ligands that are top ranked in the ligand-activity analysis
#--- Active target gene inference
active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = deg_test_c6, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                  ligand_target_matrix = ligand_target_matrix,
                                                                  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()

rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

p_ligand_target_network_6 <- vis_ligand_target %>% make_heatmap_ggplot('Prioritized ligands', 'Predicted target genes',color = 'purple', legend_position = 'top',legend_title = 'Regulatory potential') + 
  theme(axis.text.x = element_text(face = 'italic')) + scale_fill_gradient2(low = 'whitesmoke', high = 'purple')

p_ligand_target_network_6

# Save Supplemental Figure 6B

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_target_network_6
dev.off()

#--- Receptors of top-ranked ligands
lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from, to)

best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df <- lr_network_top_df_large %>% spread('from', 'weight', fill = 0)

lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors <- dist(lr_network_top_matrix, method = 'binary')
hclust_receptors <- hclust(dist_receptors, method = 'ward.D2')
order_receptors <- hclust_receptors$labels[hclust_receptors$order]

dist_ligands <- dist(lr_network_top_matrix %>% t(), method = 'binary')
hclust_ligands <- hclust(dist_ligands, method = 'ward.D2')
order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]

rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()

p_ligand_receptor_network_6 = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot('Ligands', 'Receptors', color = 'mediumvioletred', x_axis_position = 'top', legend_title = 'Prior interaction potential')

p_ligand_receptor_network_6

# Save Figure 7E

dirs <- setwd("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis")
tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_receptor_network_6
dev.off()

#--- Sender-focused approach
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands
potential_ligands_focused <- intersect(potential_ligands, gex_sender) 

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr_6 <- make_heatmap_ggplot(vis_ligand_aupr,
                                       "Prioritized ligands", "Ligand activity", 
                                       legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr_6

# Save Supplemental Figure 6B

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_aupr_6
dev.off()

#--- CCC between clusters 4, 11, 13, 15 vs 14 ----
DefaultAssay(CL_harmony_all) <- 'RNA'
Idents(CL_harmony_all) <- CL_harmony_all@meta.data$seurat_clusters

# a. Defining the sender and the receiver population
sender <- c('4', '11', '13', '15')
receiver <- c('14')

# b. Get expressed genes from sender and receiver
list_gex_sender <- sender %>% unique() %>% lapply(get_expressed_genes, CL_harmony_all, pct = 0.25)

gex_sender <- list_gex_sender %>% unlist() %>% unique()

list_gex_receiver <- receiver %>% unique() %>% lapply(get_expressed_genes, CL_harmony_all, pct = 0.25)

gex_receiver <- list_gex_receiver %>% unlist() %>% unique()

background_gex <- gex_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# c. Define a set of potential expressed ligands
expressed_ligands <- intersect(ligands, gex_sender)
expressed_receptors <- intersect(receptors, gex_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands

# d. Nichenetr ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = deg_test_c14, 
                                               background_expressed_genes = background_gex, 
                                               ligand_target_matrix = ligand_target_matrix, 
                                               potential_ligands = potential_ligands)

ligand_activities <- arrange(ligand_activities, -aupr_corrected)

best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

DotPlot(CL_harmony_all, features = best_upstream_ligands, cols = 'RdYlBu') + RotatedAxis()

# e. Infer receptors and top-predicted target genes of ligands that are top ranked in the ligand-activity analysis
#--- Active target gene inference
active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = deg_test_c14, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                  ligand_target_matrix = ligand_target_matrix,
                                                                  cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()

rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

p_ligand_target_network_14 <- vis_ligand_target %>% make_heatmap_ggplot('Prioritized ligands', 'Predicted target genes',color = 'purple', legend_position = 'top',legend_title = 'Regulatory potential') + 
  theme(axis.text.x = element_text(face = 'italic')) + scale_fill_gradient2(low = 'whitesmoke', high = 'purple')

p_ligand_target_network_14

# Save Figure 7F

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_target_network_14
dev.off()

#--- Receptors of top-ranked ligands
lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from, to)

best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df <- lr_network_top_df_large %>% spread('from', 'weight', fill = 0)

lr_network_top_matrix <- lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors <- dist(lr_network_top_matrix, method = 'binary')
hclust_receptors <- hclust(dist_receptors, method = 'ward.D2')
order_receptors <- hclust_receptors$labels[hclust_receptors$order]

dist_ligands <- dist(lr_network_top_matrix %>% t(), method = 'binary')
hclust_ligands <- hclust(dist_ligands, method = 'ward.D2')
order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]

rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()

p_ligand_receptor_network_14 = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot('Ligands', 'Receptors', color = 'mediumvioletred', x_axis_position = 'top', legend_title = 'Prior interaction potential')

p_ligand_receptor_network_14

# Save Figure 7E

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_receptor_network_14
dev.off()

#--- Sender-focused approach
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands
potential_ligands_focused <- intersect(potential_ligands, gex_sender) 

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr_14 <- make_heatmap_ggplot(vis_ligand_aupr,
                                        "Prioritized ligands", "Ligand activity", 
                                        legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr_14

# Save Figure 7F

tiff("C:/Users/Vlad/Documents/Carpenter_lab/TCR_LTBI_project/analysis", width = 5*1500, height = 5*500, res = 300, pointsize = 1)     
p_ligand_aupr_14
dev.off()



































































































































































































































































































































































































































