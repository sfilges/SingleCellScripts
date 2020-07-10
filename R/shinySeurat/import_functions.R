#----------------// Seurat shiny functions //-------------------
#
# Author: Stefan Filges
#
#
# 


#' Import Cell Ranger sample
#'
#'
#' @param path Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X
#'
#'
#' @export
#'
#' @importFrom Seurat Read10X CreateSeuratObject
#'
importCellRangerSample <- function(path, name, min_cells = 5){
  
  object <- Seurat::Read10X(data.dir = path)
  object <- Seurat::CreateSeuratObject(object, project = name, min.cells = min_cells)
  object$group <- name
  
  return(object)
}

#' Merge multiple Seurat object based on file paths
#'
#' @param df Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files provided by 10X
#'
#'
#' @export
#'
#' @import Seurat
#'
mergeSeuratObjects <- function(df, min_cells = 5){
  
  objects_to_merge <- list()
  
  for(i in 1:nrow(df)){
    row <- df[i,]
    
    #print(row)
    
    objects_to_merge <- append(
      x = objects_to_merge, 
      values = importCellRangerSample(
        path = row$path,
        name = row$name,
        min_cells = min_cells
      )
    )
  }
  
  #objects_to_merge <- c(wt_scf,wt_xen, fd_scf, fd_xen)
  # print(objects_to_merge)
  
  # Merged dataset
  merged_data <- merge(
    x = objects_to_merge[[1]],
    y = objects_to_merge[[2:length(objects_to_merge)]]
  )
  
  # store mitochondrial percentage in object meta data
  merged_data[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object = merged_data, 
    pattern = "^MT-"
  )
  
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- Seurat::cc.genes$s.genes
  g2m.genes <- Seurat::cc.genes$g2m.genes
  
  # Assign cell cycle score to genes, which will be stored in the Seurat object
  # metdata
  merged_data <- Seurat::CellCycleScoring(
    object = merged_data,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )
  
  
  unique(x = sapply(X = strsplit(x = colnames(
    x = merged_data), split = '_'), FUN = '[', 1))
  
  print(table(merged_data$orig.ident))
  
  return(merged_data)
  
}


qc_scatter <- function(object){
  
  plot1 <- Seurat::FeatureScatter(
    object = object,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    group.by = "group"
  ) + theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
  
  plot2 <- Seurat::FeatureScatter(
    object = object,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "group"
  ) + theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
  
  return(plot1 + plot2)
  
}


