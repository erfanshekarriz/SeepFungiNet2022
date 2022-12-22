library("MicrobiomeStat")
library(tidyverse)
library(phyloseq)

set.seed(1234)
### LOAD DATA
physeq18SFungi <- readRDS("./data/Data/physeq18SFungi.rds")
sampdata <- data.frame(sample_data(physeq18SFungi)) %>%
  mutate(Seep = if_else(ROV=="ROV5", "Control", "Seep"))
sample_data(physeq18SFungi) <- sampdata
taxdf <- data.frame(tax_table(physeq18SFungi))
otumatrix <- data.frame(t(otu_table(physeq18SFungi))) 


# PERFORM 
lindares <- linda(feature.dat = otumatrix,
                  meta.dat = sampdata, 
                  formula = '~Seep', 
                  alpha = 0.05,
                  feature.dat.type = c('count'),
                  prev.filter = 0.05, 
                  mean.abund.filter = 0, 
                  max.abund.filter = 0,
                  is.winsor = TRUE,
                  outlier.pct = 0.03,
                  adaptive = TRUE,
                  corr.cut = 0.1,
                  pseudo.cnt = 0.5,
                  p.adj.method = "BH", 
                  n.cores = 4)

# This indicates the variable the log fold is relative to 
lindares$variables

# extract results and label 
lindaresdf <- rownames_to_column(lindares$output$SeepSeep)  %>% 
  dplyr::left_join(., rownames_to_column(taxdf), by="rowname") %>%
  add_column(Label = NA) %>%
  mutate(Group = if_else(log2FoldChange < 0, "Control", "Seep")) 


# AUTOMATIC PLOTS
# lindaplots <- linda.plot(lindares, c('Smokey', 'Sexmale'), 
#            titles = c(''), alpha = 0.1, lfc.cut = 1,
#            legend = TRUE, directory = NULL, width = 11, height = 8) 

# lindaplots[[1]][[1]]
# lindaplots[[2]][[1]]


signifdf  <- lindaresdf %>%
  filter(padj<0.05) %>%
  mutate(Label=paste0(Genus, " sp."))
signifdf[2, "Label"] <- "Cryphonectriaceae sp."


signifdf %>%
  mutate(log2FoldChange = -log2FoldChange) %>%
  ggplot(aes(x=log2FoldChange, y=reorder(Label, +log2FoldChange))) +
  geom_point(aes(fill = reorder(Group, -log2FoldChange), 
               color = reorder(Group, -log2FoldChange))) + 
  geom_errorbarh(aes(xmin=log2FoldChange-lfcSE, 
                     xmax=log2FoldChange+lfcSE), 
                 linewidth=0.4, height=0.4) + 
  theme_bw() + 
  xlab("Log-fold Change (Log2)") + 
  scale_fill_manual(values=c("green3", "grey")) +
  scale_color_manual(values=c("green3", "grey")) +
  theme(axis.title.x=element_text(face="bold", size=8, vjust=-3), 
        axis.title.y=element_blank(),
        axis.text.y=element_text(size = 8, face ="bold.italic", color="black"), 
        axis.text.x=element_text(size = 6, face ="bold", color="black"), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_blank(),
        legend.position = "none",
        panel.border=element_rect(colour="black", linewidth = 1), 
        panel.background=element_blank(), 
        plot.background=element_blank(), 
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines")) + 
  scale_x_continuous(expand=c(0,0.3))

ggsave("./data/graphs/Fig2B_14LinDA.tiff",
       width = 7,
       height = 5.5,
       units = "cm",
       dpi = 300 )




