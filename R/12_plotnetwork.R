library(tidyverse)
library(phyloseq)
library(igraph)
library(ggridges)
library('brainGraph')


set.seed(123)
### UPLOAD FILES ###
files <- list.files(path = "./data/networks/filtered/nonweighted", 
                    full.names = TRUE)
networklist <- list()
remove.0v <- function(ingraph) {
  outgraph <- ingraph
  cat("\nOriginal Graph Vertices: ", length(V(ingraph)))
  components <- igraph::clusters(ingraph, mode="weak")
  biggest_cluster_id <- which.max(components$csize)
  vert_ids <- V(ingraph)[components$membership != biggest_cluster_id]
  outgraph <- delete.vertices(outgraph, vert_ids) # remove them
  cat("\nFiltered Graph Vertices: ", length(V(outgraph)), "\n")
  return(outgraph)
}

for (file in files){
  networklist <- append(networklist, 
                        list(remove.0v(readRDS(file = file))))
}


names(networklist) <- c("B", "BA", "BAF", "BF", "Fu")
