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
    
    print(row)
    
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
  
  print(objects_to_merge)
  
  # Merged dataset
  merged_data <- merge(
    x = objects_to_merge[[1]],
    y = objects_to_merge[[2:length(objects_to_merge)]]
  )
  
  return(merged_data)
  
}




