# #########################################
# This script reads and filters sc10x  data
# #########################################


## LOAD AND PREPARE THE DATA

## @knitr loadData

# Create the seurat object (RNA only)
#h5seurat_file = file.path(PATH_ANNDATA_H5AD_FILE)
#Convert(h5seurat_file, dest = "h5seurat", overwrite = TRUE)
#sc10x = LoadH5Seurat(PATH_ANNDATA_SEURAT_FILE);

library(Seurat)
Raw_data = Read10X(data.dir=file.path(PATH_ANNDATA_SEURAT_FILE))
metadata = read.csv(file.path(PATH_ANNDATA_SEURAT_FILE, 'metadata.csv'))

sc10x = CreateSeuratObject(counts = Raw_data, meta.data = metadata, min.cells = 3, min.features = 200)


# Attribute a numeric ID to each cell (easier than barcode)
sc10x[["numID"]] = 1:length(Cells(sc10x));

# ...................
# Heat Shock data
# ...................

# Load the list of 512 genes associated with heat shock and stress in article from Campbell et al.
hs_stress_genes_df = read.table( file = PATH_HS_STRESS_MARKER_GENES_TABLE_FILE, header = TRUE, sep=",", quote = NULL)

# Convert list of genes from human symbol to mouse symbol
converted_genes_df = convertHumanGeneList( hs_stress_genes_df$gene_symbol)
hs_stress_genes_df = merge( hs_stress_genes_df, converted_genes_df, by.x = "gene_symbol", by.y="human_symbol")
hs_stress_genes_df = hs_stress_genes_df[ order( hs_stress_genes_df$PValue, hs_stress_genes_df$logFC, decreasing = c(FALSE,TRUE)), ]

DT::datatable( hs_stress_genes_df, caption = "List of HS and Stress associated genes in human")

MODULES_GENES = list()
MODULES_GENES[[ "HeatShock"]] = unique( hs_stress_genes_df[ , "symbol"])
MODULES_GENES[[ "HeatShock_Top40"]] = unique( hs_stress_genes_df[ 1:40, "symbol"])


# FILTER DATA
# -----------

## @knitr filterData_selection


if( QC_EXPLORATION_MODE == TRUE){
  cat("<p style='color:red;'><b>WARNING: QC EXPLORATION MODE is active : Cells will not be filtered but only marked as different</b></p>")
}

### Identify mitocondrial genes in matrix
mito.genes = grep( pattern = "^mt-", x = rownames(GetAssayData(object = sc10x, slot = "counts")), value = TRUE, ignore.case = TRUE)
if(length(mito.genes)==0) 
{
  warning( "No mitochondrial genes could be identified in this dataset.");
} else 
{
  # Compute percentage of mitochondrial transcripts in each cell
  percent.mito <- PercentageFeatureSet(sc10x, features=mito.genes)
  # Add the mitocondrial gene percentage as meta information in the Seurat object
  sc10x[["percent.mito"]] <- percent.mito
}

### Identify ribosomal genes in matrix
ribo.genes = grep(pattern = "^rps|^rpl",  x = rownames(GetAssayData(object = sc10x, slot = "counts")), value = TRUE, ignore.case=TRUE)
if(length(ribo.genes)==0) 
{
  warning( "No ribosomal genes could be identified in this dataset.");
} else 
{
  # Compute percentage of ribosomal transcripts in each cell
  percent.ribo <- PercentageFeatureSet(sc10x, features=ribo.genes)
  # Add the ribosomal gene percentage as meta information in the Seurat object
  sc10x[["percent.ribo"]] <- percent.ribo
}


### Identify cells that will be rejected based on specified thresholds

# Reject cells based on UMI numbers
nUMI.drop = logical( length(Cells(sc10x)));
if( ! is.null( FILTER_UMI_MIN))
{
  nUMI.drop = nUMI.drop | (sc10x[["nCount_RNA", drop=TRUE]] < FILTER_UMI_MIN);
}
if( ! is.null( FILTER_UMI_MAX))
{
  nUMI.drop = nUMI.drop | (sc10x[["nCount_RNA", drop=TRUE]] > FILTER_UMI_MAX);
}
sc10x = AddMetaData( sc10x, metadata = nUMI.drop[ Cells( sc10x)], col.name = "outlier.nCount_RNA")

# Reject cells based on number of expressed genes
nGene.drop = logical( length(Cells(sc10x)));
if( ! is.null( FILTER_FEATURE_MIN))
{
  nGene.drop = nGene.drop | (sc10x[["nFeature_RNA", drop=TRUE]] < FILTER_FEATURE_MIN);
}
if( ! is.null( FILTER_FEATURE_MAX))
{
  nGene.drop = nGene.drop | (sc10x[["nFeature_RNA", drop=TRUE]] > FILTER_FEATURE_MAX);
}
sc10x = AddMetaData( sc10x, metadata = nGene.drop[ Cells( sc10x)], col.name = "outlier.nFeature_RNA")

# Identify cells with high percentage of mitocondrial genes
mito.drop = logical( length(Cells(sc10x)));
if( length(mito.genes) && (! is.null(FILTER_MITOPCT_MAX)))
{
  mito.drop = (sc10x[["percent.mito", drop=TRUE]] > FILTER_MITOPCT_MAX);
}
sc10x = AddMetaData( sc10x, metadata = mito.drop[ Cells( sc10x)], col.name = "outlier.percent.mito")

# Identify cells with low percentage of ribosomal genes
ribo.drop = logical( length(Cells(sc10x)));
if( length(ribo.genes) && (! is.null(FILTER_RIBOPCT_MIN)))
{
  ribo.drop = (sc10x[["percent.ribo", drop=TRUE]] < FILTER_RIBOPCT_MIN);
}
sc10x = AddMetaData( sc10x, metadata = ribo.drop[ Cells( sc10x)], col.name = "outlier.percent.ribo")


### Plot distributions of #UMIs, #Genes, %Mito, and %Ribo among cells

# The metrics used for filtering will be plotted using ggplot rasterized graphs when having to show a lot of points (instead of plotly interactive graphs)
if(length( Cells( sc10x)) < PLOT_RASTER_NBCELLS_THRESHOLD)
{
  # Do interactive plot using plotly (can cause overhead when viewing many point)
  
  # Gather data to be visualized together (cell name + numID + metrics)
  cellsData = cbind( "Cell" = colnames( sc10x), # cell names from rownames conflict with search field in datatable
                     sc10x[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                     "percent.mito" = if(length( mito.genes)) as.numeric(format(sc10x[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                     "percent.ribo" = if(length( ribo.genes)) as.numeric(format(sc10x[["percent.ribo", drop = TRUE]], digits = 5)) else NULL);
  
  # Create text to show under cursor for each cell
  hoverText = do.call(paste, c(Map( paste, 
                                    c( "", 
                                       "Cell ID: ", 
                                       "# UMIs: ", 
                                       "# Genes: ", 
                                       if(length( mito.genes)) "% Mito: ", 
                                       if(length( ribo.genes)) "% Ribo: "), 
                                    cellsData, 
                                    sep=""), 
                               sep="\n"));
  
  panelWidth =  190; # 800 when using subplot
  panelHeight = 400;
  # Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
  lypanel_umis  = plotViolinJitter( cbind( cellsData, outliers = nUMI.drop), 
                                    xAxisFormula = ~as.numeric(1), 
                                    yAxisFormula = ~nCount_RNA, 
                                    pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                    hoverText = hoverText, 
                                    xTicklabels = "# UMIs", 
                                    thresholdHigh = FILTER_UMI_MAX, 
                                    thresholdLow = FILTER_UMI_MIN, 
                                    panelWidth = panelWidth, 
                                    panelHeight = panelHeight);
  
  lypanel_genes = plotViolinJitter( cbind( cellsData, outliers = nGene.drop), 
                                    xAxisFormula = ~as.numeric(1), 
                                    yAxisFormula = ~nFeature_RNA, 
                                    pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                    hoverText = hoverText, 
                                    xTicklabels = "# Genes", 
                                    thresholdHigh = FILTER_FEATURE_MAX, 
                                    thresholdLow = FILTER_FEATURE_MIN, 
                                    panelWidth = panelWidth, 
                                    panelHeight = panelHeight);
  
  lypanel_mitos = if(length(mito.genes)) plotViolinJitter( cbind( cellsData, outliers = mito.drop), 
                                                           xAxisFormula = ~as.numeric(1), 
                                                           yAxisFormula = ~percent.mito, 
                                                           pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                                           hoverText = hoverText, 
                                                           xTicklabels = "% Mito", 
                                                           thresholdHigh = FILTER_MITOPCT_MAX, 
                                                           panelWidth = panelWidth, 
                                                           panelHeight = panelHeight) else NULL;
  
  
  lypanel_ribos = if(length(ribo.genes)) plotViolinJitter( cbind( cellsData, outliers = ribo.drop), 
                                                           xAxisFormula = ~as.numeric(1), 
                                                           yAxisFormula = ~percent.ribo, 
                                                           pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), 
                                                           hoverText = hoverText, 
                                                           xTicklabels = "% Ribo", 
                                                           thresholdLow = FILTER_RIBOPCT_MIN, 
                                                           panelWidth = panelWidth, 
                                                           panelHeight = panelHeight) else NULL;
  
  # Set panels as a list and define plotly config
  panelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
  panelsList = lapply(panelsList, config, displaylogo = FALSE, 
                      toImageButtonOptions = list(format='svg'), 
                      modeBarButtons = list(list('toImage'), 
                                            list('zoom2d', 'pan2d', 'resetScale2d')));
  
  # # Group plotly violin/jitter panels (for sizing, it uses layout of one of the plot)
  # plotPanels = layout( subplot( lypanel_umis,
  #                               lypanel_genes,
  #                               lypanel_mitos),
  #                      showlegend = FALSE) # Remove eventual legends (does not mix well with subplot)
  
  # Control layout using flex because subplot is limited regarding plot sizing and alignment
  # 'browsable' required in console, not in script/document
  browsable(
    div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
         div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
              lapply(panelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
  )
  
} else
{
  # Do the plots as simple images with ggplot when having a lot of points
  
  # Generate plotly violin/jitter panels for #umis, #genes, %mitochondrial, and %ribosomal stats
  ggpanel_umis  = ggplot( cbind(sc10x[["nCount_RNA"]], outliers = nUMI.drop), aes( y = nCount_RNA)) + 
    geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
    geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
    geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
                color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
                alpha = 0.5,
                size = 1) +
    labs(x = "# UMIs", y = "") +
    theme( axis.text.x = element_blank(), 
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));
  
  ggpanel_genes = ggplot( cbind(sc10x[["nFeature_RNA"]], outliers = nGene.drop), aes( y = nFeature_RNA)) + 
    geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
    geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
    geom_hline( yintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
                color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
                alpha = 0.5,
                size = 1) +
    labs(x = "# Genes", y = "") +
    theme( axis.text.x = element_blank(), 
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444"));
  
  ggpanel_mitos = if(length(mito.genes)) ggplot( cbind(sc10x[["percent.mito"]], outliers = mito.drop), aes( y = percent.mito)) + 
    geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
    geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
    geom_hline( yintercept = FILTER_MITOPCT_MAX, 
                color = "red", 
                alpha = 0.5,
                size = 1) +
    labs(x = "% Mito", y = "") +
    theme( axis.text.x = element_blank(), 
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;
  
  
  ggpanel_ribos = if(length(ribo.genes)) ggplot( cbind(sc10x[["percent.ribo"]], outliers = ribo.drop), aes( y = percent.ribo)) + 
    geom_jitter( aes(x = 1.5, col = outliers), width = 0.45, height = 0, size = 0.5, alpha = 0.3) +
    geom_violin( aes(x = 0.5), width = 0.9, fill = "darkblue", col = "#AAAAAA", alpha = 0.3) + 
    geom_hline( yintercept = FILTER_RIBOPCT_MIN, 
                color = "blue", 
                alpha = 0.5,
                size = 1) +
    labs(x = "% Ribo", y = "") +
    theme( axis.text.x = element_blank(), 
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none") +
    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#444444")) else NULL;
  
  # Use patchwork library to assemble panels
  print( ggpanel_umis + ggpanel_genes + ggpanel_mitos + ggpanel_ribos + plot_layout( nrow = 1))
}


cat( "<br>Number of cells removed based on number of UMIs:", sum( nUMI.drop));
cat( "<br>Number of cells removed based on number of genes:", sum( nGene.drop));
if(exists( "mito.drop")) cat( "<br>Number of cells removed based on high percentage of mitochondrial transcripts:", sum( mito.drop));
if(exists( "ribo.drop")) cat( "<br>Number of cells removed based on low percentage of ribosomal transcripts:", sum( ribo.drop));
cat( "\n<br>\n");

# Identify cells to exclude as union of cells with low nb UMI, low nb expressed genes, high pct mito genes, low pct ribo genes
sc10x[["outlier"]] = nUMI.drop  | 
  nGene.drop | 
  ( if(exists( "mito.drop")) mito.drop else FALSE ) | 
  ( if(exists( "ribo.drop")) ribo.drop else FALSE );

cat("<br><br>Removed cells after filters:", sum( unlist(sc10x[["outlier"]] )));
cat("<br>Remaining cells after filters:", sum( ! unlist(sc10x[["outlier"]] )));
cat("\n<br>\n");

### Record which cells got rejected

# Export the excluded cells to file
write.table( data.frame( cellid = names(which(nUMI.drop))), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, 
                              paste0( outputFilesPrefix, "excluded_cells_nbUMI.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

write.table( data.frame( cellid = names(which(nGene.drop))), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbGene.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

if(exists( "mito.drop"))
{
  write.table( data.frame( cellid = names(which(mito.drop))), 
               file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctMito.txt")), 
               quote = FALSE, 
               row.names = FALSE, 
               col.names = TRUE, 
               sep="\t");
}

if(exists( "ribo.drop"))
{
  write.table( data.frame( cellid = names(which(ribo.drop))), 
               file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctRibo.txt")), 
               quote = FALSE, 
               row.names = FALSE, 
               col.names = TRUE, 
               sep="\t");
}




## @knitr filterData_summaryPlot

### Plot dispersion of excluded and non-excluded cells

# number of genes and number of UMIs
ggplot( sc10x[[]][order( sc10x[["outlier"]]),], # Plot FALSE first and TRUE after
        aes( x = nFeature_RNA, 
             y = nCount_RNA, 
             color = outlier)) + 
  geom_point( size = 0.5) +
  geom_vline( xintercept = c( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_FEATURE_MIN, FILTER_FEATURE_MAX), is.null)], 
              alpha = 0.5) +
  geom_hline( yintercept = c( FILTER_UMI_MIN, FILTER_UMI_MAX), 
              linetype = 2, 
              color = c( "blue", "red")[!sapply( list( FILTER_UMI_MIN, FILTER_UMI_MAX), is.null)], 
              alpha = 0.5) +
  scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
  labs( x = "# Genes", y = "# UMIs") +
  theme( legend.position = "none")

# Mitochondrial vs ribosomal distributions
if(exists( "mito.drop") && exists( "ribo.drop"))
{
  ggplot( sc10x[[]][order( sc10x[["outlier"]]),], # Plot FALSE first and TRUE after
          aes( x = percent.ribo, 
               y = percent.mito, 
               color = outlier)) + 
    geom_point( size = 0.5) +
    geom_vline( xintercept = FILTER_RIBOPCT_MIN, linetype = 2, color = "blue", alpha = 0.5) +
    geom_hline( yintercept = FILTER_MITOPCT_MAX, linetype = 2, color = "red", alpha = 0.5) +
    scale_color_manual( values = c( "TRUE" = "#FF0000AA", "FALSE" = "#44444444")) + 
    labs( x = "% Ribosomal genes", y = "% Mitochondrial genes") +
    theme( legend.position = "none")
}


# FILTER DATA
# --------------

## @knitr filterData_filterObject

# Filter the excluded cells in the Seurat object
sc10x_nonFiltered = sc10x;
if( QC_EXPLORATION_MODE == FALSE){
  sc10x = sc10x[ , ! sc10x[[ "outlier", drop=TRUE ]] ];
}
# Save the list of remaining cells after selection during loading, and #Genes, #UMIs, pctMito, pctRibo
write.table( data.frame( cellid = Cells(sc10x)), 
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "selected_cells.txt")), 
             quote = FALSE, 
             row.names = FALSE, 
             col.names = TRUE, 
             sep="\t");

# 
# # DETECT DOUBLET
# # --------------
# 
# ## @knitr detect_doublet
# 
# library(scDblFinder)
# 
# ### Read in data as an sce object ###
# sce <- SingleCellExperiment( list( counts=Read10X( PATH_10XFILES)))
# 
# ## Calculate doublet ratio ###
# doublet_ratio <- ncol(sce)/1000*0.008
# 
# ### Calculate Singlets and Doublets ###
# sce <- scDblFinder(sce, dbr=doublet_ratio)
# 
# cell_class = sce$scDblFinder.class
# names( cell_class) = colnames( sce)
# 
# sc10x = AddMetaData( sc10x, metadata = cell_class[ Cells( sc10x)], col.name = "scDblFinder")
# 
# table( sc10x$outlier, sc10x$scDblFinder)


# NORMALIZE DATA
# --------------

## @knitr normalizeData

sc10x = NormalizeData( object = sc10x,
                       normalization.method = DATA_NORM_METHOD,
                       scale.factor = DATA_NORM_SCALEFACTOR,
                       verbose = TRUE)

#sc10x = ScaleData( object    = sc10x,
#                   do.center = DATA_CENTER,
#                   do.scale  = DATA_SCALE,
#                   vars.to.regress = DATA_VARS_REGRESS,
#                   verbose = TRUE)


