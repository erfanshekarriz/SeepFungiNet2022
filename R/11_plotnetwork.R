library(tidyverse)
library(phyloseq)
library(igraph)

set.seed(123)
### UPLOAD FILES ###
files <- list.files(path = "./data/networks/filtered/nonweighted", 
                    full.names = TRUE)

filesshort <- gsub(".rds", "", list.files(path = "./data/networks/filtered/nonweighted", 
                                          full.names = FALSE))
networklist <- list()
remove.0v <- function(ingraph) {
  outgraph <- ingraph
  cat("\nOriginal Graph Vertices: ", length(V(ingraph)))
  vert_ids <- V(ingraph)[degree(ingraph)<1]
  outgraph <- delete.vertices(outgraph, vert_ids) # remove them
  cat("\nFiltered Graph Vertices: ", length(V(outgraph)), "\n")
  return(outgraph)
}

for (file in files){
  networklist <- append(networklist, 
                        list(readRDS(file = file)))
}
names(networklist) <- c("B", "BA", "BAF", "BF", "Fu")

for (i in 1:length(networklist)){
  testnet <- remove.0v(networklist[[i]])
  degree <- degree(testnet)
  shape <- data.frame(Domain=V(testnet)$Domain) %>%
    mutate(Domain= if_else(Domain=="Eukaryota", "square", "circle")) %>%
    pull()
  color <- data.frame(Domain=V(testnet)$Domain) %>%
    mutate(Domain= if_else(Domain=="Eukaryota", "red", Domain), 
           Domain= if_else(Domain=="Archaea", "yellow", Domain), 
           Domain= if_else(Domain=="Bacteria", "green", Domain)) %>%
    pull()
  
  tiff(paste0("./data/graphs/networkfig/", filesshort[i], ".tiff"), units="cm", width=10, height=10, res=1000)
  plot.igraph(testnet, 
       vertex.label=NA,
       vertex.color = color, 
       vertex.shape = shape,
       vertex.size=degree*1.5)
  dev.off()
}



