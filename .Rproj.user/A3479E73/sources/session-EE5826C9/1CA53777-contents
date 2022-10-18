library(phyloseq)
library(tidyverse)
library(vegan)
library(ggpubr)



#### LOAD DATA ####
physeq16S <- readRDS('./data/physeq16S.rds')
physeq18S <- readRDS('./data/physeq18SFungi.rds')



# metadata 
sampledataplotraw <- data.frame(sample_data(physeq16S)) %>%
  as.data.frame() %>%
  mutate(SampleID = rownames(.))
features <- colnames(sampledataplotraw) %>% as.data.frame()


#### ANALYSE ####
### 16S 
min(sample_sums(physeq16S))
otutable_16S <- data.frame(otu_table(physeq16S))
set.seed(1234)
rareAlpha16S <- vegan::rarefy(otutable_16S, 
                              min(sample_sums(physeq16S)))


### 18S Fungi 
min(sample_sums(physeq18S))
otutable_18S <- data.frame(otu_table(physeq18S))
set.seed(1234)
rareAlpha18S <- vegan::rarefy(otutable_18S, 
                              min(sample_sums(physeq18S)))



#### CLEAN & MERGE DATA ####
rareAlpha16Splot <- rareAlpha16S %>%
  as.data.frame() %>%
  rename(., richness= `.`) %>%
  rownames_to_column(var = "SampleID") %>%
  inner_join(sampledataplotraw, by = "SampleID") %>%
  group_by(ROV) %>%
  mutate(meanrichness = mean(richness))

rareAlpha18Splot <- rareAlpha18S %>%
  as.data.frame() %>%
  rename(., richness= `.`) %>%
  rownames_to_column(var = "SampleID") %>%
  inner_join(sampledataplotraw, by = "SampleID") %>%
  group_by(ROV) %>%
  mutate(meanrichness = mean(richness))


rareAlphaALL <- bind_rows(list(`16S`= rareAlpha16Splot, 
                               `18S`= rareAlpha18Splot), 
                          .id = "Group")
write_csv(rareAlphaALL, "./data/rareifiedAlphaDiversity.csv")
rm(rareAlpha16Splot, rareAlpha18Splot)

#### PLOT ####

# ROV level diversity 
rareAlphaALL %>%
  ggplot(aes(x = reorder(ROV, +meanrichness), 
             y = richness, 
             color = ROV)) +
  geom_jitter(size = 0.01, 
              color = "gray") + 
  # geom_violin(color = "gray") + 
  geom_boxplot(outlier.shape = NA, width = 0.3) + 
  stat_compare_means(method = "kruskal.test", size = 2.8, color = "black")+
  scale_color_brewer(palette = "Set1") + 
  theme_bw() + 
  facet_wrap(~Group, scales = 'free_y') + 
  xlab("Sample Collection Station")

# ROV level methane concentrations
rareAlphaALL %>%
  ggplot(aes(x = reorder(ROV, +Methane), 
             y = Methane, 
             fill = ROV)) +
  # geom_violin(color = "gray") + 
  geom_col(outlier.shape = NA, width = 0.3) + 
  stat_compare_means(method = "kruskal.test", size = 2.8, color = "black")+
  scale_fill_brewer(palette = "Set1") + 
  theme_bw() + 
  facet_wrap(~Group, scales = 'free_y') + 
  xlab("Sample Collection Station")



# Effect of Methane levels on Diversity
rareAlphaALL %>%
  ggplot(aes(x = reorder(ROV, +Methane), 
             y = Methane, 
             fill = ROV)) +
  # geom_violin(color = "gray") + 
  geom_col(outlier.shape = NA, width = 0.3) + 
  stat_compare_means(method = "kruskal.test", size = 2.8, color = "black")+
  scale_fill_brewer(palette = "Set1") + 
  theme_bw() + 
  facet_wrap(~Group, scales = 'free_y') + 
  xlab("Sample Collection Station")




rareAlphaALL %>%
  ggplot(aes(x = Depth, 
             y = richness)) +
  geom_jitter(size = 0.1, 
              aes(color = ROV)) + 
  theme_bw() + 
  geom_smooth(method="lm", 
              formula = y~x, 
              se = FALSE) +
  facet_wrap(~Group, scales = 'free_y') 

  
