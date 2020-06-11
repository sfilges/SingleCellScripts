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

## Required packages

The analysis depends on the following R packges. Analysing networks with Cytoscape also requires an installation of the cytoscape software as well the "StringApp" which can be installed using the cytoscape app manager.

```{r packages, message=FALSE, warning=FALSE}
# Packages for main analysis, clustering, data structures etc.
library(Seurat, quietly = TRUE))
library(mclust, quietly = TRUE))
library(biomaRt, quietly = TRUE))
library(msigdbr, quietly = TRUE))
library(clusterProfiler, quietly = TRUE))
library(destiny, quietly = TRUE))
library(SummarizedExperiment, quietly = TRUE))
require(clusterExperiment, quietly = TRUE))
library(MAST, quietly = TRUE)
library(DESeq2, quietly = TRUE))
library(sctransform, quietly = TRUE)

# Library for network analysis
library(RCy3, quietly = TRUE)

# Packages for data handeling and plotting
library(tidyverse, quietly = TRUE)
library(reshape, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(viridis, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(VennDiagram, quietly = TRUE)
library(DT, quietly = TRUE)

# Packages for pseudotime ordering
library(slingshot, quietly = TRUE)
library(monocle, quietly = TRUE)
library(SCORPIUS, quietly = TRUE)
```

## Detailed walkthrough

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

#### Quality control and filtering

Mitochondrial gene expression can be a confounding factor and very high
mtDNA expression might indicate dead cells which should be removed from further
analysis.

We plot nFeature_RNA (= number of genes per cell), nCount_RNA (= number of 
transcripts per cell) and percentage of mtDNA expression and filter cells
accordingly. When performing QC variables should be considered jointly.
For instance high mtDNa expression may also reflect cells with high 
respiratory activity rather than lysed cells.

```r
plot1 <- FeatureScatter(
  object = merged_data,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  group.by = "group"
) + theme(
  legend.position = "bottom",
  legend.title = element_blank()
)

plot2 <- FeatureScatter(
  object = merged_data,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "group"
) + theme(
  legend.position = "bottom",
  legend.title = element_blank()
)

plot1 + plot2
```

From the scatter plot of transcripts counts vs. mtDNA we can see that cells
with very high mtDNA content $>10\%$ also have very low transcript numbers and
are therefore probably lysed cells and should be removed.

```r
# Filter data to remove dead cells, outliers
merged_data <- subset(
  x = merged_data,
  subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 10
)

```

#### Cell cycle scoring

```r
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign cell cycle score to genes, which will be stored in the Seurat object
# metdata
merged_data <- Seurat::CellCycleScoring(
  object = merged_data,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

VlnPlot(
  object = merged_data, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  ncol = 3, 
  group.by = "Phase"
)

VlnPlot(
  object = merged_data, 
  features = c('S.Score', 'G2M.Score'), 
  ncol = 3, 
  group.by = "group"
)
```

#### Normalizing the data using regularized negative binomial regression

We normalise the data either using [sctransform](https://www.biorxiv.org/content/10.1101/576827v1)
or using the standard Seurat workflow depending on the setting of the 
sctransform variable.

sctransform models the expression of each gene as a negative binomial random variable with a mean that depends on other variables. Here the other variables can be used to model the differences in sequencing depth between cells and are used as independent variables in a regression model. In order to avoid overfitting, we will first fit model parameters per gene, and then use the relationship between gene mean and parameter values to fit parameters, thereby combining information across genes. Given the fitted model parameters, we transform each observed UMI count into a Pearson residual which can be interpreted as the number of standard deviations an observed count was away from the expected mean. If the model accurately describes the mean-variance relationship and the dependency of mean and latent factors, then the result should have mean zero and a stable variance across the range of expression. During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage.

We can also assign each cell a cell cycle score based on known cell-cycle
associated genes for S and G2/M phases. These are available through the Seurat
package.

We perform the normalization workflow and only regress out the variable 'percent.mt'.
We perform PCA on the scaled data with and withoutcell cycle regressed genes, however in both cases we regress out the mtDNAexpression. For this PCA we only use the annotated cell cycle genes! It is expected that cells cluster according to cell cycle stage without regression and that
regressed data shows much less separation according to cell cycle stage (but
not necessarily zero).

To filter cell cycle genes too set the filter_cell_cycle variable to TRUE.

This document will always show the effect of cell cycle regression by plotting the 
samples in PCA space with and without filtering, filter_cell_cycle only determines
which of the two will be used downstream.

```r
# Regress out mitochondrial expression and cell cycle stage
vars.to.regress = c('percent.mt', 'S.Score', 'G2M.Score')

# Decide which normalization workflow to use
if(sctransform){
  merged_data <- SCTransform(
    object = merged_data, 
    vars.to.regress = vars.to.regress,
    verbose = FALSE
  )
} else {
  merged_data <- Seurat::NormalizeData(
    object = merged_data,
    verbose = FALSE
  )
  
  merged_data <- Seurat::FindVariableFeatures(
    object = merged_data, 
    selection.method = 'vst', 
    nfeatures = 2000
  )
  
  # Scale data without regressing out cell cycle
  merged_data_cc <- Seurat::ScaleData(
    object = merged_data,
    vars.to.regress = 'percent.mt'
  )
  
  # Regress out cell cycle
  merged_data_no_cc <- Seurat::ScaleData(
    object = merged_data,
    vars.to.regress = vars.to.regress
  )
  
  # Decide to use data with/without cell cycle regression
  if(filter_cell_cycle){
    merged_data <- merged_data_no_cc
  } else {
    merged_data <- merged_data_cc
  }
  
  # Perform PCA on cc regressed samples
  merged_data_cc <- Seurat::RunPCA(
    object = merged_data_cc,
    features = c(s.genes, g2m.genes)
  )
  
  before_cc_correction <- Seurat::DimPlot(merged_data_cc)
  
  merged_data_no_cc <- Seurat::RunPCA(
    object = merged_data_no_cc,
    features = c(s.genes, g2m.genes)
  )
  
  after_cc_correction <- Seurat::DimPlot(merged_data_no_cc)
}
```


