library(selbal)
library(phyloseq)
library(tidyverse)
library(metagMisc)
library(grid)
library(ggpubr)


set.seed(1234)

#### LOAD DATA ####
physeq <- readRDS("./data/Data/physeq18SFungi.rds")
physeq.filt <- phyloseq_filter_prevalence(physeq, 
                           prev.trh = 0.1, 
                           abund.trh = NULL,
                           threshold_condition = "OR", 
                           abund.type = "total")
taxdf <- data.frame(tax_table(physeq.filt))
otumat <- as.matrix(data.frame(otu_table(physeq.filt)))
sampledf <- data.frame(sample_data(physeq.filt)) %>% mutate(Methane = Methane/1000)
methane <- as.numeric(sampledf[,"Methane"])
# ggplot() + geom_histogram(aes(methane), bins=10, fill="blue", color="white") + theme_bw()


avgDepth0Meth <- sampledf %>%
  filter(Depth==0) %>%
  pull(Methane) %>%
  max()

methlev <- sampledf %>% 
  mutate(MethLev = if_else(Methane<=avgDepth0Meth, "Low", "High")) %>%
  pull(MethLev) %>% as.factor()

ROVcov <- sampledf %>% dplyr::mutate(ROV=as.factor(ROV)) %>% dplyr::select(ROV)


#### SELBAL MODELLING WITH DISCRETICIZED METHANE LEVELS ####

selbalmod <- selbal.cv(x = otumat,
                       y = methlev,
                       n.fold = 5,
                       n.iter = 2000,
                       logit.acc = "AUC")
saveRDS(selbalmod, "./data/Data/selbalMETHLev10000.rds")



#### VISUALIZE THE BEST FIT MODEL ####
selbalmod <- readRDS("./data/Data/selbalMETHLev.rds")
# selbalmod$accuracy.nvar
# grid.draw(selbalmod$global.plot)
# plot.tab(selbalmod$cv.tab)
# selbalmod$cv.accuracy; summary(selbalmod$cv.accuracy)
# selbalmod$var.barplot
# selbalmod$glm
# selbalmod$opt.nvar

globalbal <- selbalmod$global.balance %>%
  mutate(Taxa=gsub("^X", "", Taxa)) %>%
  merge(., taxdf %>% rownames_to_column("Taxa"), by="Taxa") 

# manually name by BLASTing representative sequence against the NCBI database
manualnames <- c("Unkown Fungal ASV 3", 
                 "Derxomyces sp.",
                 "Unkown Fungal ASV 1", 
                 "Moniliella mellis", 
                 "Candida-Lodderomyces sp.", 
                 "Moniliella mellis", 
                 "Unkown Fungal ASV 2",
                 "Fusarium solani",
                 "Fusarium oxysporum")

globalbal <- globalbal %>%
  mutate(Label=manualnames)


#### CLEAN PLOTTING & SAVING ####
FINALresdf <- arrange(selbalmod$var.barplot$data, -sel)[1:selbalmod$opt.nvar, ] %>%
  mutate(name=gsub("^X", "", name)) %>%
  rename(Taxa=name) %>%
  merge(., globalbal, by="Taxa") %>%
  column_to_rownames("Taxa")

selbalmod$glm$data %>%
  mutate(Methane = if_else(numy==1, "Low", "High"),
         Methane=factor(Methane, levels=c("Low", "High"))) %>%
  rename(Balance=V1) %>%
  ggplot(aes(x=Methane, y=Balance)) + 
  geom_jitter(width=0.1, size=0.1, color="grey30", height=0) + 
  geom_boxplot(aes(color=Methane), fill=NA, 
               size=0.6, width=0.2, outlier.shape = NA) + 
  stat_boxplot(geom ='errorbar', 
               aes(color=Methane),
               size=0.6, width=0.2) + 
  stat_compare_means(aes(label = ..p.signif..), 
                     comparisons = list(c("High", "Low")), 
                     size = 5,
                     vjust = 0.2,
                     method = "wilcox.test", 
                     color = "black", 
                     bracket.size = 0.5, 
                     label.y = 8.2) +
  scale_y_reverse() + 
  theme_bw() +
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 9, vjust = -2), 
        axis.title.y = element_text(size = 9, vjust= 4), 
        legend.text =   element_text(size = 9), 
        legend.title =   element_text(size = 9),
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_rect(fill = "grey95", color = "black"),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey99"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  labs(x = "Methane Concentration\nRelative to Surface Layer",
       y = "Cross Validated Balance Values") + 
  scale_color_manual(values=c("#377EB8", "#E41A1C")) + 
  guides(color=guide_legend(title="Relative \nMethane \nConcentration"))


ggsave("./data/graphs/Fig5A_21selbal.tiff",
       width = 10,
       height = 7.5,
       units = "cm",
       dpi = 1000 )



FINALresdf %>%
  mutate(Group = if_else(Group=="NUM", "Low", "High")) %>%
  ggplot(aes(x=sel, y = reorder(Label, +sel))) +
  geom_point(aes(color=Group), size=2) + 
  theme_bw() +
  theme(axis.text =  element_text(size = 10, face="italic"), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_blank(), 
        legend.text =   element_text(size = 9), 
        legend.title =   element_text(size = 9),
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_rect(fill = "grey95", color = "black"),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey99"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  labs(x = "Appearance in Balance across \n10,000 Cross-Validations (%)") + 
  scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + 
  scale_x_continuous(limits = c(0,80)) +
  guides(color=guide_legend(title="Relative \nMethane \nConcentration")) 


ggsave("./data/graphs/Fig5B_21selbal.tiff",
       width = 12,
       height = 8,
       units = "cm",
       dpi = 1000 )

