# SingleCellScripts

Collection of scripts for processing and analysis of single cell RNA seq generated with 10x genomics systems or similar technologies.

## Overview

The R folder contains the following scripts:

- 1_seurat_pre_processing.Rmd
- 2_seurat_dim_reduction.Rmd
- 3_seurat_gsea_and_diff_expression.Rmd
- 4_pseudotime_analysis.Rmd
- 5_pseudotime_plots.Rmd
- 6_cytoscape_networks.R

## Detailed analysis

### (1) Preprocessing data

The preprocessing script first loads the Seurat and tidyverse packages and sets a few basic parameters for the analysis.

```r
# run sctransform, this replaces the NormalizeData, ScaleData and FindVariableFeatures
sctransform = FALSE

# Should cell cycle filtering be applied or not?
filter_cell_cycle = FALSE

# define the working directory
data_directory <- "~/GitHub/SingleCellScripts/10x/"

# define output file
seurat_object_save <- "~/GitHub/SingleCellScripts/data/Rdata/seurat_object_preprocessed.Rdata"

suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
```

Next we load the count matrices for each individual sample and then merge them into an single objecy for combined analysis.

```r
#---------------- scaffold cultures -------------------
# Analyze FUS-DDIT3 samples
fd_scf <- Read10X(paste(data_directory,"/HT1080_Scf_FD/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
fd_scf <- CreateSeuratObject(fd_scf,project = "scf_FD", min.cells = 5)
fd_scf$group <- "scf-FD"

# Analyze WT samples
wt_scf <- Read10X(paste(data_directory,"HT1080_Scf_WT/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
wt_scf <- CreateSeuratObject(wt_scf,project = "scf_WT", min.cells = 5)
wt_scf$group <- "scf-WT"

#---------------- xenograft cultures -------------------
# Analyze WT samples
wt_xen <- Read10X(paste(data_directory,"HT1080_Xen_WT/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
wt_xen <- CreateSeuratObject(wt_xen,project = "xen_WT", min.cells = 5)
wt_xen$group <- "xenograft-WT"

# Analyze FUS-DDIT3 samples
fd_xen <- Read10X(paste(data_directory,"HT1080_Xen_FD/HT1080_FD_mouse_custom_hg38/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
fd_xen <- CreateSeuratObject(fd_xen,project = "xen_FD", min.cells = 5)
fd_xen$group <- "xenograft-FD"

# Merged dataset
merged_data <- merge(
  x = wt_scf,
  y = c(wt_xen, fd_scf, fd_xen),
  add.cell.ids = c('wt-scf', 'wt-xen', 'fd-scf', 'fd-xen'),
  project = 'scaffold'
)

unique(x = sapply(X = strsplit(x = colnames(
  x = merged_data), split = '_'), FUN = '[', 1))

table(merged_data$orig.ident)

# store mitochondrial percentage in object meta data
merged_data[["percent.mt"]] <- PercentageFeatureSet(
  object = merged_data, 
  pattern = "^MT-"
)
```




