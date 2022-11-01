library(visNetwork)
library(igraph)
library(tidyverse)


igraphnetw <- readRDS("./data/networks/filtered/weighted/BAF_igraphSLR_igraphw.rds")
igraphnet <- readRDS("./data/networks/filtered/nonweighted/BAF_igraphSLR_igraph.rds")
igraphnet <- delete.vertices()

edgelist <- as_edgelist(igraphnet) %>% data.frame() %>%
  rename(from=X1, to=X2) %>% 
  cbind(get.edge.attribute(igraphnetw))

nodedf <- data.frame(id=as_ids(V(igraphnetw)), 
                     degree=degree(igraphnetw), 
                    betweenness=betweenness(igraphnet)) %>%
  cbind(as.data.frame(get.vertex.attribute(igraphnetw))) %>%
  mutate(shape= str_replace_all(Domain, "Eukaryota", "triangle"))



visNetwork(nodedf, edgelist, width = "100%") %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)




