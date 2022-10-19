library(phyloseq)
library(tidyverse)





#### UPLOAD GRAPHS ####
files <- list.files(path = "./data/networks/raw",
                    pattern = "igraph",
                    recursive = TRUE,
                    full.names = TRUE)

graphlist <- lapply(files, readRDS)
names(graphlist) <- files

####
physeqcomb <- readRDS("./data/Data/physeqcombined.rds")



#### FIND OTUs in at least two sampling sites (from ROV1, ROV2, ROV3) ####
metadata <- data.frame(sample_data(physeqcomb)) %>%
  select(ROV)
otudf <- data.frame(otu_table(physeqcomb)) %>%
  mutate_all(function(x) if_else(x>0,
                                 1,
                                 0)) %>%
  cbind(metadata)

# count how many samples EACH ASV appears in EACH SITE 
sitecount <- otudf %>%
  group_by(ROV) %>%
  summarize_all(sum) %>%
  column_to_rownames(var = "ROV")

# select ASVs that appear in at least ROV1 and ROV2 (highest activity cold seeps)
nonprevTaxa <- sitecount %>%
  filter(rownames(.) !="ROV3") %>%
  mutate_all(function(x) if_else(x>0,
                                 1,
                                 0)) %>%
  summarise_all(sum) %>%
  select_if(function(col) sum(col) < 2) %>%
  colnames()




#### FILTER and WRITE Graphs to "filtered" networks directory ####
for (i in 1:length(graphlist)){
  intersectV <- intersect(nonprevTaxa, as_ids(V(graphlist[[i]])))
  cat("Before\n")
  cat(paste("Vertices:", length(V(graphlist[[i]]))))
  cat("\n")
  cat(paste("Edges:", length(E(graphlist[[i]]))))
  cat("\n")
  
  newgraph <- delete.vertices(graphlist[[i]], intersectV)
  cat("After\n")
  cat(paste("Vertices:", length(V(newgraph))))
  cat("\n")
  cat(paste("Edges:", length(E(newgraph))))
  cat("\n\n")
  
  saveRDS(newgraph, names(graphlist)[i] %>%
            str_replace("raw/.+/","filtered/"))
  
}
