if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

library(tidyverse)
library(RColorBrewer)
library(RCy3)

# define the working directory
working_directory <- "~/GitHub/SingleCellScripts/"

cytoscapePing()
cytoscapeVersionInfo()

top_genes <- read_csv2(paste(working_directory, "data/gene_lists/diff_genes_over_pseudotime_top_1500_q_0.01.csv", sep = ""))

string.cmd = paste('string protein query cutoff=0.4 species="Homo sapiens" limit=0 query="',
  paste(top_genes$gene, collapse=","),'"',sep="")

commandsGET(string.cmd)
commandsRun("analyzer analyze directed=FALSE")

loadTableData(
  as.data.frame(top_genes),
  data.key.column = "gene",
  table.key.column = 'display name'
)

style.name = "dataStyle"
createVisualStyle(style.name)
setVisualStyle(style.name)

sizes = c(20, 60, 100, 130)
control.points = c(0, 0.01, 0.05, 0.1)

#node.colors <- c(rev(brewer.pal(4, "Set1")))
node.colors <- c('#7CAE00', '#F8766D', '#00BFC4', '#C77CFF')

setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
setNodeSizeDefault(60, style.name)
setNodeColorDefault("#AAAAAA", style.name)
setEdgeLineWidthDefault(2, style.name)
setNodeLabelMapping('display name', style.name)
setNodeSizeMapping ('BetweennessCentrality', control.points, sizes, style.name = style.name)
setEdgeLineWidthDefault(0.0001, style.name = style.name)
setNodeColorMapping("cluster", c(1,2,3,4), node.colors, style.name=style.name)

#---------------------// Generating sub networks //----------------------

clusters <- c(1,2,3,4)

for(i in clusters){
  
  gene_selection <- dplyr::filter(top_genes, cluster == i)
  node.color <- node.colors[i]
  
  string.cmd = paste('string protein query cutoff=0.4 species="Homo sapiens" limit=0 query="',
                     paste(gene_selection$gene, collapse=","),'"',sep="")
  
  commandsGET(string.cmd)
  commandsRun("analyzer analyze directed=FALSE")
  
  loadTableData(as.data.frame(gene_selection),data.key.column = "gene",table.key.column = 'display name')
  
  style.name = paste("dataStyle",i,sep="")
  createVisualStyle(style.name)
  setVisualStyle(style.name)
  
  sizes = c(20, 60, 100, 130)
  control.points = c(0, 0.01, 0.05, 0.1)
  
  setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
  setNodeSizeDefault(60, style.name)
  setNodeColorDefault("#AAAAAA", style.name)
  setEdgeLineWidthDefault(2, style.name)
  setNodeLabelMapping('display name', style.name)
  setNodeSizeMapping ('BetweennessCentrality', control.points, sizes, style.name = style.name)
  setEdgeLineWidthDefault(0.0001, style.name = style.name)
  setNodeColorMapping("cluster", i, node.color, style.name=style.name)
  
  
  file.name <- paste(working_directory, "figures/string_network_cluster_",i,".pdf", sep="")
  #response <- exportImage(file.name, type = "pdf")
  
}








