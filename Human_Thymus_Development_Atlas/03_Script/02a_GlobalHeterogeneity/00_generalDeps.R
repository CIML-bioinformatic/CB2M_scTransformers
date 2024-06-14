# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries

library( digest)# digest (hashing)
library( DT)            # datatable
library( forcats)       # fct_inorder (factors)
library( fs)            # path_sanitize
library( future)        # plan (multicore)
library( ggplot2)
library( ggpubr)
library( ggrepel)
library( gridExtra)
library( htmltools)     # browsable
library( htmlwidgets)   # JS (for datatable)
library( iheatmapr)     # iheatmap
library( kableExtra)
library( knitr)
library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
library( dplyr)
library( rmarkdown)
library( scales)        # hue_pal
library( RColorBrewer)

# Single-cell technology
library( Seurat)
library( SeuratData)
library( SeuratDisk)
library( Matrix) # Force reload Matrix 1.5.3 instead of 1.5.1 that cause a bug in FindNeighbors function

# Functional Enrichment analysis
library( biomaRt)
library( babelgene)
library( clusterProfiler)
library( org.Mm.eg.db)
library( org.Hs.eg.db)
library( rrvgo)
library( rstatix)

