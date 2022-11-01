library(lefser)
library(mia)
library(RColorBrewer)

set.seed(1234)
### LOAD DATA
physeq18SFungi <- readRDS("./data/Data/physeq18SFungi.rds")
sampdata <- data.frame(sample_data(physeq18SFungi)) %>%
  mutate(Seep = if_else(ROV=="ROV5", "Control", "Seep"))
sample_data(physeq18SFungi) <- sampdata
taxdf <- data.frame(tax_table(physeq18SFungi))


# make input for LEFSE 
SExp <- makeTreeSummarizedExperimentFromPhyloseq(physeq18SFungi)
lefseRes <- lefser(SExp, 
                   groupCol = "Seep", 
                   #blockCol = "age_category", 
                   kruskal.threshold = 0.05,
                   wilcox.threshold = 0.05,
                   lda.threshold = 2)

lefserPlot(lefseRes)
lefseResDF <- lefseRes %>% 
  dplyr::left_join(., (taxdf %>%
              rownames_to_column(var="Names"))) %>%
  add_column(Label = NA) %>%
  mutate(Group = if_else(scores < 0, "Control", "Seep")) %>%
  arrange(scores)



### ADD NAMES MANUALLY
lefseResDF <- lefseResDF[-5,] # remove redundant taxa
lefseResDF[1,"Label"] <- "Cryphonectriaceae"
lefseResDF[2,"Label"] <- "Verticillium longisporum" 
lefseResDF[3,"Label"] <- "Unidentified sp. 1"
lefseResDF[4,"Label"] <- "Unidentified sp. 2" 
lefseResDF[5,"Label"] <- "Talaromyces purpureogenus"
lefseResDF[6,"Label"] <- "Archaeorhizomyces finlayi"
lefseResDF[7,"Label"] <- "Unidentified sp. 3"
lefseResDF[8,"Label"] <- "Verticillium longisporum"
lefseResDF[9,"Label"] <- "Papiliotrema aurea"
lefseResDF[10,"Label"] <- "Suillus lakei"
lefseResDF[11,"Label"] <- "Malassezia sp."
lefseResDF[12,"Label"] <- "Fusarium oxysporum"


brewer.pal(3, "Set1")
#### PLOT 
lefseResDF %>%
  ggplot(aes(x=scores, y=reorder(Label, +scores))) +
  geom_col(aes(fill = reorder(Group, -scores), 
               color = reorder(Group, -scores)), 
           width=0.85) + 
  theme_bw() + 
  xlab("LDA Effect Score (log10)") + 
  scale_fill_manual(values=c("green3", "grey")) + 
  scale_color_manual(values=c("green3", "grey")) + 
  theme(axis.title.x=element_text(face="bold", size=8, vjust=-3), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 6, face ="bold", color="black"), 
        axis.text.x=element_text(size = 6, face ="bold", color="black"), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_blank(),
        legend.position = "right",
        panel.border=element_rect(colour="black", size = 1), 
        panel.background=element_blank(), 
        plot.background=element_blank(), 
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines"))  + 
  scale_x_continuous(expand=c(0,2))


ggsave("./data/graphs/Fig1B_16LDAeffect.tiff",
       width = 9.5,
       height = 6.5,
       units = "cm",
       dpi = 1000 )

