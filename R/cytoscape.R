if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}


library(tidyverse)
library(RColorBrewer)
library(RCy3)

cytoscapePing()
cytoscapeVersionInfo()

top_genes <- read_csv2("diff_genes_over_pseudotime.csv")

top_genes[top_genes$cluster == 5,]$cluster <- 2

#%>%
  #dplyr::filter(cluster %in% c(2,5))

string.cmd = paste('string protein query cutoff=0.4 species="Homo sapiens" limit=0 query="',
  paste(top_genes$gene, collapse=","),'"',sep="")

commandsGET(string.cmd)
#commandsRun(string.cmd)

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

node.colors <- c(rev(brewer.pal(4, "Set1")))


setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
setNodeSizeDefault(60, style.name)
setNodeColorDefault("#AAAAAA", style.name)
setEdgeLineWidthDefault(2, style.name)
setNodeLabelMapping('display name', style.name)
setNodeSizeMapping ('BetweennessCentrality', control.points, sizes, style.name = style.name)
setEdgeLineWidthDefault(0.0001, style.name = style.name)
setNodeColorMapping("cluster", c(1,2,3,4), node.colors, style.name=style.name)



initial_string_network_png_file_name <- file.path(
  getwd(),"networks", "images", "initial_string_network.png"
)

if(file.exists(initial_string_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  response <- file.remove(initial_string_network_png_file_name)
}

response <- exportImage(initial_string_network_png_file_name, type = "png")
