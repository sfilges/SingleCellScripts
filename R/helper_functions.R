##---------------------------------// Header //------------------------------
## 
## Miscellaneous R functions
##
## This document contains various functions for data analysis and visualisation.
##
## Author: Stefan Filges
##
## Version 2020-06-16
##
##--------------------------------//  Main //----------------------------------



#' Custom doplot for gene set enrichments
#' 
#' Generates dotplots to visualise gene-set enrichments similar to enrichplot::dotplot.
#' 
#' @export
#' 
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr separate
#'
#' @param egmt object created with the clusterProfiler::enricher function
#' @param showCategory number of categories to plot, default is 20
#' @param x.axis which variable to show on the x axis, either 'qvalue' or 'GeneRatio'
#' @param font.size font size for labels
#' @param colour colour scheme to use
#' 
#' @author Stefan Filges
#' 
custom_dotplot <- function(
  egmt,
  showCategory = 20,
  x.axis = 'GeneRatio',
  font.size = 7,
  colour = 'magma'
){
  
  qval_cutoff <- egmt@qvalueCutoff
  pval_cutoff <- egmt@pvalueCutoff
  
  # filter results table
  data <- egmt@result %>%
    dplyr::arrange(qvalue) %>%
    dplyr::select(-geneID) %>%
    dplyr::filter(pvalue < pval_cutoff) %>%
    dplyr::filter(qvalue < qval_cutoff) %>%
    tidyr::separate(GeneRatio, c("top", "bottom"), sep = "/") %>%
    dplyr::mutate(GeneRatio = as.numeric(top)/as.numeric(bottom))
  
  # select top categories
  data <- head(data, showCategory)
  
  # generate plot
  data$category <- factor(data$ID, levels = rev(data$ID))
  if(x.axis == 'GeneRatio') {
    dp <- ggplot(
      data = data,
      mapping = aes(x = GeneRatio, y = category, color = -log10(qvalue))
    ) +
      geom_point(aes(size = Count)) +
      scale_color_continuous(
        low="red", high="blue",
        name = '-log10(qvalue)',
        guide=guide_colorbar(reverse=FALSE)
      ) +
      theme_minimal() +
      theme(
        axis.line = element_line(color="black", size = 0.2),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = font.size)
      ) + guides(size = guide_legend(reverse = TRUE))
  } else if (x.axis == 'qvalue'){
    dp <- ggplot(
      data = data,
      mapping = aes(x = -log10(qvalue), y = category, color = GeneRatio)
    ) +
      geom_point(aes(size = Count)) +
      #scale_color_viridis_c(option = colour, guide = guide_colorbar(reverse=FALSE)) +
      scale_color_continuous(low="red", high="blue",name = 'GeneRatio',guide=guide_colorbar(reverse=FALSE)) +
      theme_minimal() +
      theme(
        axis.line = element_line(color="black", size = 0.2),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = font.size)
      ) + guides(size = guide_legend(reverse = TRUE))
  }
  return(dp)
}


enricher_custom <- function(
  gene,
  pvalueCutoff,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  USER_DATA) {
  
  gene <- as.character(unique(gene))
  qExtID2TermID <- EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  
  if (is.null(qTermID)) {
    message("--> No gene can be mapped....")
    p2e <- get("PATHID2EXTID", envir = USER_DATA)
    sg <- unlist(p2e[1:10])
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, collapse = ","))
    message("--> return NULL...")
    return(NULL)
  }
  qExtID2TermID.df <- data.frame(
    extID = rep(names(qExtID2TermID),
                times = lapply(qExtID2TermID, length)),
    termID = qTermID
  )
  
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  qTermID2ExtID <- with(
    qExtID2TermID.df, split(as.character(extID),as.character(termID))
  )
  
  extID <- ALLEXTID(USER_DATA)
  
  if (missing(universe)) 
    universe <- NULL
  if (!is.null(universe)) {
    extID <- as.character(universe)
  }
  
  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)
  qTermID <- unique(names(qTermID2ExtID))
  termID2ExtID <- TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)
  geneSets <- termID2ExtID
  idx <- get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)
  
  if (sum(idx) == 0) {
    msg <- paste("No gene set have size >", minGSSize, "...")
    message(msg)
    message("--> return NULL...")
    return(NULL)
  }
  
  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))
  
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]
  N <- rep(length(extID), length(M))
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(
    numWdrawn = k - 1, 
    numW = M, 
    numB = N - M, 
    numDrawn = n
  )
  
  pvalues <- apply(
    args.df, 1, function(n) phyper(n[1], n[2], n[3], n[4], lower.tail = FALSE)
  )
  
  GeneRatio <- apply(
    data.frame(a = k, b = n), 1, function(x) paste(x[1], "/", x[2], sep = "", collapse = "")
  )
  
  BgRatio <- apply(
    data.frame(a = M, b = N), 1, function(x) paste(x[1], "/", x[2], sep = "", collapse = "")
  )
  
  Over <- data.frame(
    ID = as.character(qTermID), 
    GeneRatio = GeneRatio, 
    BgRatio = BgRatio, 
    pvalue = pvalues,
    stringsAsFactors = FALSE
  )
  
  p.adj <- p.adjust(Over$pvalue, method = pAdjustMethod)
  
  qobj <- tryCatch(qvalue(p = Over$pvalue, 
                          lambda = 0.05,
                          pi0.method = "bootstrap"),
                   error = function(e) NULL)
  
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  }
  else {
    qvalues <- NA
  }
  
  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse = "/"))
  geneID <- geneID[qTermID]
  
  Over <- data.frame(
    Over,
    p.adjust = p.adj,
    qvalue = qvalues,
    geneID = geneID, 
    Count = k, 
    stringsAsFactors = FALSE
  )
  
  Description <- TERM2NAME(qTermID, USER_DATA)
  
  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx, ]
  }
  
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)
  row.names(Over) <- as.character(Over$ID)
  
  x <- new("enrichResult", result = Over, pvalueCutoff = pvalueCutoff, 
           pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff, 
           gene = as.character(gene), universe = extID, geneSets = geneSets, 
           organism = "UNKNOWN", keytype = "UNKNOWN", ontology = "UNKNOWN", 
           readable = FALSE)
  
  return(x)
}

enricher_edit<-function(
  gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, minGSSize = 10, 
  maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE, TERM2NAME = NA
){
  USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
  enricher_custom(
    gene = gene, pvalueCutoff = pvalueCutoff, 
    pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
    maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = USER_DATA
  )
}

body(enricher)<-body(enricher_edit)

environment(enricher_custom)<-environment(clusterProfiler:::enricher_internal)


#' Gene set enrichment plot
#' 
#' 
#' @param data gmt object
#' @param id cluster id
#' 
#' @import monocle
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
enrichment_plot <- function(data, id){
  
  gene <- data %>% dplyr::filter(cluster == id)
  
  kegg <- read.gmt(
    gmtfile = paste(working_directory, "data/MSigDB/c2.cp.kegg.v7.1.symbols.gmt", sep = "")
  )
  cgp <- read.gmt(
    gmtfile =  paste(working_directory, "data/MSigDB/c2.cgp.v7.1.symbols.gmt", sep = "")
  )
  reactome <- read.gmt(
    gmtfile =  paste(working_directory, "data/MSigDB/c2.cp.reactome.v7.1.symbols.gmt", sep = "")
  )
  
  hgnc_symbols <- read.csv(
    file = paste(working_directory, "data/MSigDB/hgnc.txt", sep = ""),
    sep = '\t',
    header=TRUE
  )
  
  hgnc_symbols <- hgnc_symbols$Approved.symbol
  
  egmt <- clusterProfiler::enricher(
    gene$gene,
    TERM2GENE=cgp,
    pvalueCutoff = 0.1,
    pAdjustMethod = 'fdr',
    minGSSize = NA,
    maxGSSize = NA,
    qvalueCutoff = 0.05,
    universe = hgnc_symbols
  )
  cgp_dotplot <- custom_dotplot(egmt, x.axis = "qvalue")
  
  egmt <- clusterProfiler::enricher(
    gene$gene,
    TERM2GENE=kegg,
    pvalueCutoff = 0.1,
    pAdjustMethod = 'fdr',
    minGSSize = NA,
    maxGSSize = NA,
    qvalueCutoff = 0.05,
    universe = hgnc_symbols
  )
  kegg_dotplot <- custom_dotplot(egmt, x.axis = "qvalue")
  
  egmt <- enricher(
    gene$gene,
    TERM2GENE=cgp,
    pvalueCutoff = 0.1,
    pAdjustMethod = 'fdr',
    minGSSize = NA,
    maxGSSize = NA,
    qvalueCutoff = 0.05,
    universe = hgnc_symbols
  )
  cgp_dotplot <- custom_dotplot(egmt, x.axis = "qvalue")
  
  egmt <- enricher(
    gene$gene,
    TERM2GENE=reactome,
    pvalueCutoff = 0.1,
    pAdjustMethod = 'fdr',
    minGSSize = NA,
    maxGSSize = NA,
    qvalueCutoff = 0.05,
    universe = hgnc_symbols
  )
  reactome_dotplot <- custom_dotplot(egmt, x.axis = "qvalue")
  
  return((cgp_dotplot + kegg_dotplot))
}




pseudotime_hm <- function (cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2", 
          num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
          add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
          norm_method = c("log", "vstExprs"), scale_max = 3, 
          scale_min = -3, trend_formula = "~sm.ns(Pseudotime, df=3)", 
          return_heatmap = FALSE, cores = 1) 
{
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- monocle::genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
                       relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
      FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap::pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
                 cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
                 clustering_distance_rows = row_dist, clustering_method = hclust_method, 
                 cutree_rows = num_clusters, silent = TRUE, filename = NA, 
                 breaks = bks, border_color = NA, color = hmcols)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                         num_clusters)))
    
    annotation_colors = list(
      Cluster = RColorBrewer::brewer.pal(n = 4, name = "Set1")
    )
    
    names(annotation_colors$Cluster) <- unique(annotation_row$Cluster)
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
                                                               ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != 100) {
      stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
    }
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                                                      "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                                                       "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    if (!is.null(annotation_row)) 
      row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  if (!is.null(annotation_row)) 
    row.names(annotation_row) <- row_ann_labels
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  ph_res <- pheatmap::pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, show_rownames = show_rownames, 
                     show_colnames = F, clustering_distance_rows = row_dist,annotation_colors = annotation_colors,
                     clustering_method = hclust_method, cutree_rows = num_clusters, 
                     annotation_row = annotation_row, annotation_col = annotation_col, 
                     treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols, 
                     border_color = NA, silent = TRUE, filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap) {
    return(ph_res)
  }
}







