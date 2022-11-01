library(phyloseq)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(tidytree)
library(microbiome)



physeq <- readRDS("./data/Data/physeq18SFungi.rds")
taxtable <- data.frame(tax_table(physeq))
# plot_tree(physeq, "treeonly")



#### CLEAN PHY TREE ####
physeqmerge <- microbiome::transform(merge_samples(tax_glom(physeq, "Phylum"), "ROV"), "clr")
taxdf <- data.frame(tax_table(physeqmerge))
otumerge <- data.frame(otu_table(physeqmerge)) %>% 
  mutate_all(as.numeric) %>% t() %>%
  as.data.frame() %>% rownames_to_column(var="SampleID") %>%
  mutate(SampleID = gsub("X", "", SampleID)) %>%
  merge((taxdf %>% rownames_to_column("SampleID") %>% select(Phylum, SampleID)), 
        by="SampleID") %>%
  rename(ID=SampleID)
treephylum <- phyloseq::phy_tree(physeqmerge)



# plot_tree(physeqmerge, color="ROV", 
#           label.tips="Phylum", 
#           size="abundance", plot.margin=0.6)

tree <- ggtree(treephylum,
               size=0.5, 
               ladderize = TRUE)
otumergeheat <- otumerge %>% pivot_longer(c(ROV1, ROV2, ROV3, ROV5), 
                                          names_to = "ROV", 
                                          values_to = "Abund") %>%
  select(-Phylum)
treephyldf <- otumerge %>% select(ID, Phylum)

# join tree & metadata
tree %<+% treephyldf + 
  geom_tiplab(aes(label=Phylum), 
              size=2.5, 
              hjust=-0.1, 
              fontface="bold.italic") + 
  geom_fruit(
    data=otumergeheat,
    geom=geom_tile,
    mapping=aes(x=ROV, y=ID, alpha=Abund),
    fill="green3", 
    color="white",
    pwidth=0.9,
    offset=1.1, 
    axis.params = list(axis = "x", text.angle = 90, 
                       text.size = 2, 
                       text.color ="green", 
                       hjust =-14.5, 
                       title = NULL, 
                       line.alpha = 0)
  ) + 
  scale_color_brewer(palette="Set1") + 
  scale_fill_brewer(palette="Set1") + 
  theme(legend.text=element_text(size = 6, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_text(size=8, face="bold"),
        legend.position = "right",
        panel.background=element_blank(), 
        plot.background=element_blank()) + 
  scale_y_continuous(expand=c(0,3)) + 
  labs(alpha='CLR \nAbundance\n') + 
  scale_alpha_continuous(range=c(0.0,1))


ggsave("./data/graphs/Fig1A_19phylumtree.tiff",
       width = 10,
       height = 11.7,
       units = "cm",
       dpi = 1000 )



