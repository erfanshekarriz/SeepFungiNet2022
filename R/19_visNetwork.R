library(visNetwork)
library(igraph)
library(tidyverse)
library(phyloseq)


igraphnetw <- readRDS("./data/networks/filtered/weighted/BAF_igraphSLR_igraphw.rds")
igraphnet <- readRDS("./data/networks/filtered/nonweighted/BAF_igraphSLR_igraph.rds")

# filter disconnected nodes
igraphnetw <- delete.vertices(igraphnetw, degree(igraphnetw)<2)
igraphnet <- delete.vertices(igraphnet, degree(igraphnet)<2)

edgelist <- as_edgelist(igraphnet) %>% data.frame() %>%
  dplyr::rename(from=X1, to=X2) %>% 
  cbind(get.edge.attribute(igraphnetw)) %>%
  mutate(value=weight, 
         color=if_else(weight>0, "green", "red"))

# nodedf <- data.frame(id=as_ids(V(igraphnetw)), 
#                      degree=degree(igraphnetw), 
#                     betweenness=betweenness(igraphnet)) %>%
#   cbind(as.data.frame(get.vertex.attribute(igraphnetw))) %>%
#   mutate(shape= str_replace_all(Domain, c("Eukaryota"="triangle", 
#                                           "Bacteria"="dot", 
#                                           "Archaea"="square"))) %>%
#   mutate(color= str_replace_all(Domain, c("Eukaryota"="red", 
#                                           "Bacteria"="green", 
#                                           "Archaea"="yellow"))) %>%
#   mutate(size= degree*4, 
#          label=gsub("_", " ", Species))


#### The node table was manually edited to classify trophic modes
#### To verify taxa to higher resolution the represantative 
#### sequences were manually BLASTED against the NCBI database
# write.csv(nodedf, "./data/Data/BAFvisNetworknode.csv")
nodedf <- read.csv("./data/Data/BAFvisNetworknode.csv", row.names = 1, na.strings=c("","NA")) 
unique(nodedf$group)
nodedf <- nodedf %>% mutate(color = if_else(group=="Methane Metabolism", "red", "white"), 
                            color = if_else(group=="Methane Balance Fungi", "red", color),
                            color = if_else(group=="Fungi Keystone Hub", "green", color),
                            color = if_else(group=="Sulfur Metabolism", "magenta", color),
                            color = if_else(group=="Sulfur & Nitrogen Metabolism", "orange", color),
                            color = if_else(group=="Nitrogen Metabolism", "yellow", color),
                            color = if_else(group=="Hydrocarbon Metabolism", "purple", color),
                            color = if_else(group=="Organohalide Metabolism", "orange", color), 
                            color = if_else(is.na(group), "grey", color))  %>%
  filter(degree>=3) %>% mutate(font.size= if_else(color=="grey", 1, 35))
unique(nodedf$color)


net <- visNetwork(nodedf, edgelist, width = "100%") %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visOptions(selectedBy = "group",
             highlightNearest = list(enabled = T, degree = 1, hover = T),
             nodesIdSelection = TRUE) %>%
  # visOptions(selectedBy = "Domain") %>%
  visInteraction(dragNodes = TRUE, 
                 dragView = TRUE, 
                 zoomView = TRUE) %>%
  visLayout(randomSeed = 1234)  %>% 
  visEdges(smooth = list(enabled=T, roundness=0.1)) %>% visInteraction(navigationButtons = TRUE)




net %>% visSave("./data/networks/InteractiveNetworkBAF.html", selfcontained = TRUE, background = "white")
net %>% visSave("./data/graphs/InteractiveNetworkBAF.html", selfcontained = TRUE, background = "white")







