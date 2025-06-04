# SingleCellScripts

Collection of scripts for processing and analysis of single cell RNA seq generated with 10x genomics systems or similar technologies.

The scripts were developed using *Seurat v4*. Detailed methods and results can be found in the [Ranji et al., Journal of Translational Medicine (2024)](https://link.springer.com/article/10.1186/s12967-024-05211-w#Sec2). The scRNA-Seq data from the study (10x Genomics) is available from GEO under accession [GSE191132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191132).

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
library(SummarizedExperiment, quietly = TRUE))
require(clusterExperiment, quietly = TRUE))
library(MAST, quietly = TRUE)
library(DESeq2, quietly = TRUE))
library(sctransform, quietly = TRUE)

# Library for network analysis
library(RCy3, quietly = TRUE)

# Packages for data handeling and plotting
library(tidyverse, quietly = TRUE)
library(viridis, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(DT, quietly = TRUE)

# Packages for pseudotime ordering
library(slingshot, quietly = TRUE)
library(monocle, quietly = TRUE)
library(SCORPIUS, quietly = TRUE)
```

## Detailed walkthrough

### (1) Preprocessing data

The preprocessing script first loads the Seurat and tidyverse packages and sets a few basic parameters for the analysis. In particular, the user must define a working directory.

```r
#----------------------- // USER defined variables //---------------------------
# define the working directory
working_directory <- "~/GitHub/SingleCellScripts/"

# run sctransform, this replaces the NormalizeData, ScaleData and FindVariableFeatures
sctransform = FALSE

# Should cell cycle filtering be applied or not?
filter_cell_cycle = FALSE
#--------------------------------------------------------------------------------------

# set other directories based on working directory
matrix_dir <- paste(working_directory,"data/count_matrices/",sep="")
rdata_dir <- paste(working_directory,"data/Rdata/",sep="")

# Check if directory for Rdata files exists, if not create it
ifelse(
  !dir.exists(rdata_dir), 
  dir.create(file.path(rdata_dir)), FALSE
)

# Store seurat object for quick re-use
seurat_object_save <- paste(rdata_dir,"seurat_object_preprocessed.Rdata",sep = "")

# Load libraries for preprocessing
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(Seurat, quietly = TRUE))
```

Next we load the count matrices for each individual sample and then merge them into an single objecy for combined analysis.

```r
#---------------- scaffold cultures -------------------
# Analyze FUS-DDIT3 samples
fd_scf <- Read10X(paste(matrix_dir,"HT1080_Scf_FD/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
fd_scf <- CreateSeuratObject(fd_scf,project = "scf_FD", min.cells = 5)
fd_scf$group <- "scf-FD"

# Analyze WT samples
wt_scf <- Read10X(paste(matrix_dir,"HT1080_Scf_WT/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
wt_scf <- CreateSeuratObject(wt_scf,project = "scf_WT", min.cells = 5)
wt_scf$group <- "scf-WT"

#---------------- xenograft cultures -------------------
# Analyze WT samples
wt_xen <- Read10X(paste(matrix_dir,"HT1080_Xen_WT/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
wt_xen <- CreateSeuratObject(wt_xen,project = "xen_WT", min.cells = 5)
wt_xen$group <- "xenograft-WT"

# Analyze FUS-DDIT3 samples
fd_xen <- Read10X(paste(matrix_dir,"HT1080_Xen_FD/outs/filtered_gene_bc_matrices/custom_egfp_hg38/",sep=""))
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


### (2) Dimensionality reduction

```r
# load object
seurat_object_save <- "~/GitHub/SingleCellScripts/data/Rdata/seurat_object_preprocessed.Rdata"
load(seurat_object_save)

# define output file
seurat_object_final <- "~/GitHub/SingleCellScripts/data/Rdata/seurat_object_final.Rdata"
```

#### Identification of highly variable features (feature selection)

Even after removing features found in only a few cells and removing artefacts
the data set may be very large with > 15,000 features. Thus, only genes highly 
variable across the data set are used for downstream analysis. According to
Luecken and Theis typical cut-offs are between 1000 and 5000 genes. Although
downstream analysis is robust as to the the number of genes chosen, they recommend
to use larger values. We used 2000 genes, the default in the Seurat package.

```r
# Identify the 10 most highly variable genes
top30 <- head(x = VariableFeatures(object = merged_data), 30)

# plot variable features with and without labels
plot1 <- Seurat::VariableFeaturePlot(object = merged_data)
Seurat::LabelPoints(plot = plot1, points = top30, repel = TRUE)
```

#### Dimensionality reduction & clustering

Despite containing many cells and features scRNA-seq data is inherently low-dimensional,
i.e. most biological variation in scRNA-seq data can be shown in few dimensions.

Dimensionality reduction techniques generate reduced representations through
linear or non-linear combination of gene expression vectors. Non-linear methods
like [tSNE](http://www.jmlr.org/papers/v9/vandermaaten08a.html) and
[UMAP](https://www.nature.com/articles/nbt.4314) result in better cell 
visualisations because they capture local similaity in high-dimensional data 
better than linear methods such as PCA. However, non-linear reduced dimensions 
are not straight-forward to interpret, wheras principal components and their 
associated gene loadings can be interpreted easily.

Dimensionality reduction techniques should be interpreted separatly: PCA or 
diffusion maps for general purpose summarization and trajectory inference, 
UMAP for exploratory visualization and UMAP with PAGA for complex data sets.

We first apply PCA to the merged data set and print the top five features (loadings) for
each of the first five principal components, both in positive and negative direction.

We visualise the PCA in two dimensions and colour cells by sample ID. The elbow plot showws the
PCs ranked by variance contribution, allowing us to determine a reduced dimensionality
for downstream processing. Typically we want to choose the number of PC around
a "hinge" in the elbow plot. Alternatively a statistical procedure called Jackstraw,
implemented in Seurat, can be used to determine the significance of PCs using
resampling. 

We also plot the top positive and negative loadings for the first six dimensions
in a heatmap.

```{r dimReduction1, message=FALSE} 
merged_data <- Seurat::RunPCA(
  object = merged_data,
  verbose = FALSE
)

elbow_plot <- ElbowPlot(
  object = merged_data
)

print(x = merged_data[['pca']], dims = 1:5, nfeatures = 5)

merged_data <- Seurat::JackStraw(
  object = merged_data,
  num.replicate = 100
)

merged_data <- Seurat::ScoreJackStraw(
  object = merged_data,
  dims = 1:15
)

jackstraw_plot <- Seurat::JackStrawPlot(
  object = merged_data,
  dims = 1:15
)

pca_dim_plot <- Seurat::DimPlot(
  object = merged_data,
  reduction = 'pca', 
  group.by = 'group'
)
pca_dim_plot

pca_dim_loadings <- Seurat::VizDimLoadings(
  object = merged_data,
  dims = 1:2, 
  reduction = 'pca',
  combine = TRUE
)
pca_dim_loadings

pca_dim_heatmap <- Seurat::DimHeatmap(
  object = merged_data,
  dims = 1:6,
  cells = 500,
  balanced = TRUE,
  fast = FALSE
)
pca_dim_heatmap
```

#### Neighbour embedding and clustering

As in PhenoGraph, we first construct a k-nearest neighbour (KNN) graph based on 
the euclidean distance in PCA space, and refine the edge weights between any two 
cells based on the shared overlap in their local neighborhoods (Jaccard similarity). 
This step is performed using the FindNeighbors function, and takes as input the 
previously defined dimensionality of the dataset (first 15 PCs).

To cluster the cells, we next apply modularity optimization techniques such as the 
Louvain or Leiden algorithms, to iteratively group cells together, with the goal of 
optimizing the standard modularity function. The FindClusters function implements this 
procedure, and contains a resolution parameter that sets the ‘granularity’ of the 
downstream clustering, with increased values leading to a greater number of clusters. 
We find that setting this parameter between 0.4-1.2 typically returns good results 
for single-cell datasets of around 3K cells. Optimal resolution often increases for 
larger datasets. The clusters can be found using the Idents function. A relatively
easy explanation of the Louvain and Leiden algorithms is bublished [here](https://www.nature.com/articles/s41598-019-41695-z).

```r
# Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset. We first 
# determine the k-nearest neighbors of each cell. We use this knn graph to construct the
# SNN graph by calculating the neighborhood overlap (Jaccard index) between every cell 
# and its k.param nearest neighbors.
merged_data <- FindNeighbors(
  object = merged_data,
  k.param = 20,
  dims = 1:15
)

# Find clusters by optimizing modularity of a shared nearest neigbour graph using
# community detection
# algorithm = 2 (Louvain)
# algorithm = 4 (Leiden)
# default resolution is 0.8; use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
merged_data <- FindClusters(
  object = merged_data,
  algorithm = 2,
  resolution = 1.2
)

merged_data <- RunUMAP(
  object = merged_data, 
  dims = 1:15
)
```

#### Cell cycle percentages

```r
dat <- merged_data@meta.data

ggplot(dat, aes(x= Phase,  group=group)) + 
    geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
    geom_text(aes( label = scales::percent(..prop..),
                   y= ..prop.. ), stat= "count", vjust = -.5) +
    labs(y = "Percent", fill="group") +
    facet_grid(~group) +
    scale_y_continuous(
      labels = scales::percent
      ) +
  theme_bw()
```

#### UMAP

```{r plots2, echo=FALSE, fig.width=6, fig.height=4}
umap1 <- Seurat::DimPlot(
  object = merged_data,
  reduction = 'umap',
  group.by = 'group'
) + theme(
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.x=element_blank(),
  axis.ticks.y=element_blank(),
  legend.position = 'right'
) + scale_color_brewer(palette = 'Dark2')
  
#scale_color_manual(values = c("#ff8668","gray70", "#5e8aab", "#ffce8e"))

umap2 <- Seurat::DimPlot(
  object = merged_data,
  reduction = 'umap',
  group.by = 'seurat_clusters',
  label = TRUE
) + theme(
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.x=element_blank(),
  axis.ticks.y=element_blank(),
  legend.position = 'none'
)

umap3 <- Seurat::DimPlot(
  object = merged_data,
  reduction = 'umap',
  group.by = 'Phase'
) + theme(
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.x=element_blank(),
  axis.ticks.y=element_blank(),
  legend.position = 'bottom'
)

umap1 | (umap2 / umap3)

save(merged_data, file = seurat_object_final)
```

