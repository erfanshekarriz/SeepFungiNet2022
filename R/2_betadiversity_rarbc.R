library(phyloseq)
library(tidyverse)
library(vegan)
library(reshape2)
library(ggpmisc)
library(corrplot)

#### THE ANALYSIS OF RAREFIED DISTANCES TAKES A LONG TIME



#### LOAD DATA ####
physeq16S <- readRDS('./data/Data/physeq16S.rds')
physeq18S <- readRDS('./data/Data/physeq18SFungi.rds')



# metadata 
sampledataplotraw <- data.frame(sample_data(physeq16S)) %>%
  as.data.frame() 

features <- colnames(sampledataplotraw) %>% as.data.frame()

### read dist matrix if already done
raredist16S <- readRDS('./data/Data/raredist_16S.rds')
raredist18S <- readRDS('./data/Data/raredist_18S.rds')

#### CALCULATE RAREDIST ####

#### 16S
# Rarefied braycurtis dist
otu16S <- data.frame(otu_table(physeq16S))
set.seed(1234)
raredist16S <- avgdist(otu16S, 
                        dmethod = "bray", 
                        sample=min(sample_sums(physeq16S)), 
                        iterations = 100)
df <- melt(as.matrix(raredist16S), varnames = c("row", "col"))
saveRDS(raredist16S, file = "raredist_16S.rds")
write.csv(df, "16S_rarefied_distmatrixOTU.csv")


#### 18S
# Rarefied braycurtis dist
otu18SFungi <- data.frame(otu_table(physeq18S))
set.seed(1234)
raredist18S <- avgdist(otu18SFungi, 
                        dmethod = "bray", 
                        sample=min(sample_sums(physeq18S)), 
                        iterations = 100)
df <- melt(as.matrix(raredist18S), varnames = c("row", "col"))
saveRDS(raredist18S, file = "raredist_18S.rds")
write.csv(df, "18S_rarefied_distmatrixOTU.csv")



#### ORDINATION & PLOT ####
### 16S
set.seed(123)
pcoa16S <- cmdscale(raredist16S, 
                    k=2, 
                    eig = TRUE)

plotdf16S<- as.data.frame(pcoa16S$points) %>%
  as_tibble(rownames = "SampleID") %>%
  mutate(PCo1 = V1, PCo2 = V2) %>%
  select(-V1, -V2) %>%
  inner_join(., 
             sampledataplotraw %>%
               rownames_to_column(var = "SampleID"), 
             by = "SampleID") 


plotdf16S %>%
  ggplot(aes(x = PCo1, y = PCo2, color = ROV)) +
  geom_point() + 
  theme_bw()


### 18S Fungi 
set.seed(123)
pcoa18S <- cmdscale(raredist18S, 
                    k=2, 
                    eig = TRUE)

plotdf18S <- as.data.frame(pcoa18S$points) %>%
  as_tibble(rownames = "SampleID") %>%
  mutate(PCo1 = V1, PCo2 = V2) %>%
  select(-V1, -V2) %>%
  inner_join(., 
             sampledataplotraw %>%
               rownames_to_column(var = "SampleID"), 
             by = "SampleID") 
plotdf18S %>%
  ggplot(aes(x = PCo1, y = PCo2, color = as.factor(ROV))) +
  geom_point() + 
  theme_bw()


