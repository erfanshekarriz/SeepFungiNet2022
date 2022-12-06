library(phyloseq)
library(tidyverse)
library(metagMisc)


# SET PREV THERSHOLD
prevthreshold <- 0.2


#### Prepare data ####
inputdir <- getwd()
files <- list.files(path= './data/', 
                    pattern=".rds", 
                    full.names=TRUE, 
                    recursive=FALSE) 

physeqlist <- list()
for (file in files){
  physeqlist <- append(physeqlist, readRDS(file))
  
}
names(physeqlist) <- (files %>%
                        str_match(pattern = "//(.*).rds"))[,2]



#### Get a COMBINED metadata data frame ####
metadatacom <- lapply(physeqlist,
                      FUN = function(x) {data.frame(sample_data(x)) %>%
                          dplyr::mutate(SampleID = rownames(.))}
) %>% 
  reduce(inner_join, by = NULL) %>%
  as_tibble() %>%
  column_to_rownames(var = "SampleID") 

commonsamp <- rownames(metadatacom)


#### FILTER & MAKE PHYLOSEQ ####

physeqlistfilt <- list()
for (i in 1:length(physeqlist)){
  # common sample
  physeqlist[[i]] <- subset_samples(physeqlist[[i]] , 
                                    sample_names(physeqlist[[i]] ) %in% commonsamp)
  sample_data(physeqlist[[i]]) <- metadatacom
  
  # filter step that you can customize
  physeqfilt <- phyloseq_filter_prevalence(physeqlist[[i]], 
                                           prev.trh = prevthreshold, 
                                           abund.trh = NULL,
                                           threshold_condition = "OR", 
                                           abund.type = "total")
  physeqlistfilt <- append(physeqlistfilt, physeqfilt)
  names(physeqlistfilt)[i] <- names(physeqlist)[i]
  
  # save the phyloseq files with a SPECIFIC name for makeNetwork.R 
  physeqname <- paste(names(physeqlistfilt)[i], 
                      "_NetREADY_",
                      "taxprevl_", prevthreshold, 
                      ".rds",
                      sep='')
  saveRDS(physeqlistfilt, paste("./data/", physeqname, sep=""))
  
  # check how many taxa were filtered
  cat(paste("\n\nPrevelance filtering done on", names(physeqlistfilt)[i]))
  cat(paste("\nBefore filtering:", ntaxa(physeqlist[[i]]), "taxa..."))
  cat(paste("\nAfter filtering:", ntaxa(physeqfilt), "taxa..."))
}

