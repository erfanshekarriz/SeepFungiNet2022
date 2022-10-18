library(phyloseq)
library(tidyverse)


# input phyloseq and column to split on
path <- "./data/physeq18SFungi.rds"
physeq <- readRDS(path)

# retrieve split factors
metadata <- data.frame(sample_data(physeq)) %>%
  mutate(Seep = if_else(condition = (ROV %in% c("ROV1", "ROV2", "ROV3")), 
                        true="seep", false="non-seep"))

metadatasplit <- metadata %>%
  split(f = as.factor(.$Seep))  %>%  # change the factor here to split on!
  lapply(., FUN = function(x){rownames(x)})



# automatic splitting and saving 
for (i in seq_along(metadatasplit)){
  physeqsub <- subset_samples(physeq, 
                              sample_names(physeq) %in% metadatasplit[[i]])
  saveRDS(physeqsub, 
          paste(str_sub(path, end=-5), "_", names(metadatasplit)[i], ".rds", sep = ""))
}
