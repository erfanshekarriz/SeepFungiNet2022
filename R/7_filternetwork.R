library(phyloseq)
library(tidyverse)





#### UPLOAD GRAPHS 
files <- list.files(path = "./data/networks/raw",
                    pattern = "igraph",
                    recursive = TRUE,
                    full.names = TRUE)

graphlist <- lapply(files, readRDS)
names(graphlist) <- files

####
physeqcomb <- readRDS("./data/Data/physeqcombined.rds")



#### FIND OTUs in at least two sampling sites (ROV1, ROV2, ROV3)
metadata <- data.frame(sample_data(physeqcomb))
otudf <- data.frame(otu_table(physeqcomb)) %>%
  if_else(.>0, 1, 0) %>%
  cbind(metadata)


twositeotus <- otudf %>%
  group_by(ROV) %>%
  



