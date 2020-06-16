##-------------------------------------// HEADER //-------------------------------
##
## Cytoscape networks script
##
##
## PURPOSE
##
## Create networks using cytoscape based on a list of input genes. The script
## colours nodes in the network based on assigment of a cluster identity (or any other variable)
## to each gene.
##
##
## REQUIREMENTS
##
## Requires installation of cyctoscape from https://cytoscape.org/. This script uses
## the R package RCy3, which uses the CyREST API to communicate with Cytoscape using
## external scripts (R, Python, etc.). Therefore, Cytoscape needs to be installed
## and run in the background in order for the script to work. Further requires installtion
## of the StringApp with the Cytoscape app manager or the R command below.
## 
##
## Version: 16.06.2020
##
## Author: Stefan Filges
## Email: stefan.filges@gu.se 
##
##
##-------------------------------------// MAIN //-------------------------------
#
# Install RCy3 if necessary
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

# Load required packages
library(tidyverse)
library(RColorBrewer)
library(RCy3)

# StringApp can also be installed from R:
#installApp("stringApp") # http://apps.cytoscape.org/apps/stringapp

# define the working directory
working_directory <- "~/GitHub/SingleCellScripts/"

# check connection to cytoscape
RCy3::cytoscapePing()
RCy3::cytoscapeVersionInfo()

# Import gene list and cluster assignment
top_genes <- readr::read_csv2(paste(working_directory, "data/gene_lists/diff_genes_over_pseudotime_top_1500_q_0.01.csv", sep = ""))
top_genes

# Define command to run in cytoscape using default parameters:
# (confidence) cutoff=0.4
# (Maximum additional interactors) limit = 0
string.cmd = paste('string protein query cutoff=0.4 species="Homo sapiens" limit=0 query="',
  paste(top_genes$gene, collapse=","),'"',sep="")

# Execute command and analyze the networks
RCy3::commandsGET(string.cmd)
RCy3::commandsRun("analyzer analyze directed=FALSE")

RCy3::loadTableData(
  as.data.frame(top_genes),
  data.key.column = "gene",
  table.key.column = 'display name'
)

# Define visualization options
style.name = "dataStyle"
createVisualStyle(style.name)
setVisualStyle(style.name)

# Define node sizes
sizes = c(20, 60, 100, 130)
# Define centrality cut-off to choose which size a node will have
control.points = c(0, 0.01, 0.05, 0.1)

# Colour nodes based on clusters
node.colors <- c(brewer.pal(4, "Set1"))
#node.colors <- c('#7CAE00', '#F8766D', '#00BFC4', '#C77CFF')

setNodeShapeDefault("ellipse", style.name) #remember to specify your style.name!
setNodeSizeDefault(60, style.name)
setNodeColorDefault("#AAAAAA", style.name)
setEdgeLineWidthDefault(2, style.name)
setNodeLabelMapping('display name', style.name)
setNodeSizeMapping ('BetweennessCentrality', control.points, sizes, style.name = style.name) # node size depends on BetweennessCentrality
setEdgeLineWidthDefault(0.0001, style.name = style.name)
setNodeColorMapping("cluster", c(1,2,3,4), node.colors, style.name=style.name)


#---------------------// SUB NETWORKS //----------------------
#
# Generate networks using the gene list for each cluster independently.
#
# Define variables (clusters)
clusters <- c(1,2,3,4)

# Iterate through each cluster and generate the network
for(i in clusters){
  # Select genes belonging to the i-th cluster
  gene_selection <- dplyr::filter(top_genes, cluster == i)
  
  # Pick colour to match the complete network above
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
}








