# ##################################################
# The aim of this script is to find the most
# variable genes among the cell population
# ##################################################


## FIND THE MOST VARIABLE GENES USING THE SEURAT PACKAGE METHOD

## @knitr findVariableGenes_seuratMethod

# Find most variable features
sc10x = FindVariableFeatures( object = sc10x, 
                              selection.method = "vst",
                              loess.span = 0.3,
                              clip.max = "auto",
                              mean.function = ExpMean, 
                              dispersion.function = LogVMR,
                              nfeatures = VARIABLE_FEATURES_MAXNB,
                              verbose = TRUE);

write.csv(VariableFeatures(sc10x), file =file.path(PATH_ANALYSIS_OUTPUT, "Variable_Gene.csv"))