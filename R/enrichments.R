
##### REMEMBER TO SET YOUR WORKING DIRECTORY TO SOURCE FILE LOCATION IN RSTUDIO: Session-->Set Working Directory-->To source file location

library(dplyr)
library(VennDiagram)
library(msigdbr)
library(clusterProfiler)
library(xlsx)

## If the packages are not installed, try to run: install.packages("BiocManager")
## Then, run: BiocManager::install(c("dplyr","VennDiagram","msigdbr","clusterProfiler","xlsx"), update=TRUE, ask=FALSE)
## If there still are problems with packages installation, it is likely you need to install extra dependencies or libraries, contact me.

##### Load anntoation and background data:

GO_BP<-read.gmt("c5.bp.v7.0.symbols.gmt")
Reactome<-read.gmt("c2.cp.reactome.v7.0.symbols.gmt")
Hallmarks<-read.gmt("h.all.v7.0.symbols.gmt")
TFT<-read.gmt("c3.tft.v7.0.symbols.gmt")

hgnc_symbols<-read.csv("hgnc.txt", sep = "\t", header=T)
hgnc_symbols<-hgnc_symbols$Approved.symbol

##### Customize enricher function to match msigdbr results:

# We need to modifiy the enrichment analysis function to match more  closely the algorythm used in the msigdb website.
# The main difference is the definition of the universe (background), which in the GSEA web tool is defined as all the 
# genes in the NCBI database with a common gene name annotation, while in this package it is defined as all the genes that
# take part in any of the categories of the database analyzed.

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

enricher_edit<-function(gene, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe, 
                         minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, TERM2GENE, 
                         TERM2NAME = NA) 
{
  USER_DATA <- build_Anno(TERM2GENE, TERM2NAME)
  enricher_custom(gene = gene, pvalueCutoff = pvalueCutoff, 
                  pAdjustMethod = pAdjustMethod, universe = universe, minGSSize = minGSSize, 
                  maxGSSize = maxGSSize, qvalueCutoff = qvalueCutoff, USER_DATA = USER_DATA)
}

body(enricher)<-body(enricher_edit)

environment(enricher_custom)<-environment(clusterProfiler:::enricher_internal)

##### Load your data:

# You should have your data structured as .txt file, you can do this just by copying the column of gene names that you are interested in and pasting in wordpad or similar.
# If you want to run this script automatically, I recoomend you have all the files including the .R script, the annotation files and the .txt data named as data.txt file 
# in a folder. You can just copy and paste this folder multiple time for the different ananlysis and just substitute your data.txt file. Always add a identifier for your 
# sample or experiment in the first row of the data.txt file, which will be used for the output name files. If you don't want to have this in your output filenames, use a
# single dot as your first row of the data.txt file, OTHERWISE THE FIRST GENE ID IN THE LIST WILL BE USED FOR THE FILENAMES INSTEAD!!!!

data_table<-read.csv("data.txt", header = TRUE, check.names = F)
sample_name<- gsub(pattern = "^\\.$", replacement = "", x =names(data_table)) #Save sample name for output naming. In case you hd a single dot in the first row, no sample name is used in the output files. 
data<-as.character(data_table[,1])
data<-gsub(pattern = "\\s*",replacement = "",x = data, perl = T) #Remove spaces from ids in case there were any
data<-gsub(pattern = "\\..*",replacement = "",x = data, perl = T) #Remove isoforms annotation if needed


# The commands below automatically save the resulting plot in a .png format in the Results folder. You can change the resolution modifying the width and heigth parameters
# in the "png(...)" commands.

##### Hallmarks enrichment analysis:

Hallmarks_enrichment<- enricher(
  gene = data,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  minGSSize = NA,
  maxGSSize = NA,
  qvalueCutoff = 0.05,
  TERM2GENE = Hallmarks,
  universe = hgnc_symbols
)

df=as.data.frame(Hallmarks_enrichment)
df$GeneRatio<-sapply(df$GeneRatio, function(x) eval(parse(text=x)))
df$BgRatio<-sapply(df$BgRatio, function(x) eval(parse(text=x)))
term=df$ID[order(df$GeneRatio, decreasing=T)]
png(paste(if (sample_name==""){"Results/Hallmarks enrichment"}else{"Results/Hallmarks enrichment "},sample_name,".png",sep=""), width = 1080, height = 940)
print(dotplot(Hallmarks_enrichment, showCategory=term[1:20], title="Hallmarks enrichment", color="qvalue")) # You can change the range of categories plotted, by default it will be the top 20 categories with a higher GeneRatio. You can also change the title of the plot in the the "main" argument.
dev.off()

##### Reactome enrichment analysis:

Reactome_enrichment<- enricher(gene = data, pvalueCutoff = 0.05, pAdjustMethod = "fdr", minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.05, TERM2GENE = Reactome, universe = hgnc_symbols)

df=as.data.frame(Reactome_enrichment)
df$GeneRatio<-sapply(df$GeneRatio, function(x) eval(parse(text=x)))
df$BgRatio<-sapply(df$BgRatio, function(x) eval(parse(text=x)))
term=df$ID[order(df$GeneRatio, decreasing=T)]
png(paste(if (sample_name==""){"Results/Reactome enrichment"}else{"Results/Reactome enrichment "},sample_name,".png",sep=""), width = 1080, height = 940)
print(dotplot(Reactome_enrichment, showCategory=term[1:20], title="Reactome enrichment", color="qvalue")) # You can change the range of categories plotted, by default it will be the top 20 categories with a higher GeneRatio. You can also change the title of the plot in the the "main" argument.
dev.off()

##### Gene ontology (biological processes) enrichment analysis:

GO_BP_enrichment<- enricher(gene = data, pvalueCutoff = 0.05, pAdjustMethod = "fdr", minGSSize = NA, maxGSSize = NA, qvalueCutoff = 0.05, TERM2GENE = GO_BP, universe = hgnc_symbols)

df=as.data.frame(GO_BP_enrichment)
df$GeneRatio<-sapply(df$GeneRatio, function(x) eval(parse(text=x)))
df$BgRatio<-sapply(df$BgRatio, function(x) eval(parse(text=x)))
term=df$ID[order(df$GeneRatio, decreasing=T)]
png(paste(if (sample_name==""){"Results/GO BP enrichment"}else{"Results/GO BP enrichment "},sample_name,".png",sep=""), width = 1080, height = 940)
print(dotplot(GO_BP_enrichment, showCategory=term[1:20], title="Gene ontology (biological processes) enrichment", color="qvalue")) # You can change the range of categories plotted, by default it will be the top 20 categories with a higher GeneRatio. You can also change the title of the plot in the the "main" argument.
dev.off()

##### Save the results per category in a excel friendly format:

for (i in ls(pattern = "*.enrichment$")) {
  filename<-gsub(pattern = "_", replacement = " ", x = i)
  if (sample_name==""){
    write.xlsx(eval(parse(text = i)), paste("Results/",filename,".xls",sep=""), row.names = F)
  } else{
    write.xlsx(eval(parse(text = i)), paste("Results/",filename," ",sample_name,".xls",sep=""), row.names = F)
  }
}
