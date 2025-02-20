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
