library(tidyverse)
library(ggVennDiagram)
library(RColorBrewer)


physeq18SFungi <- readRDS("./data/Data/physeq18SFungi.rds")
physeqfilt <- merge_samples(physeq18SFungi, "ROV")
otudf <- data.frame(otu_table(physeqfilt)) %>% t() %>% data.frame()
taxdf <- data.frame(tax_table(physeqfilt))
sampdf <- data.frame(sample_data(physeqfilt))


x <- list(
  ROV1 = otudf$ROV1[otudf$ROV1>5], 
  ROV2 = otudf$ROV2[otudf$ROV2>5], 
  ROV3 = otudf$ROV3[otudf$ROV3>5],
  ROV5 = otudf$ROV5[otudf$ROV5>5]
)


venn <- Venn(x)
data <- process_data(venn)
ggplot() +
  geom_sf(aes(alpha=count), data = venn_region(data), fill="skyblue") +
  geom_sf(aes(color=name), size = 0.8, lty = "solid", 
          data = venn_setedge(data), show.legend = F) +
  geom_sf_text(aes(label = name, color=name), data = venn_setlabel(data), size=3, 
               fontface="bold", vjust=-0.5) +
  geom_sf_text(aes(label=paste0(signif(count*100/sum(count),2), "%")), 
               fontface = "italic", data = venn_region(data), 
               size=3) +
  geom_sf_text(aes(label=count), 
               fontface = "italic", data = venn_region(data), 
               vjust=-1, size=3) +
  theme_void() + 
  scale_color_brewer(palette="Set1") + 
  scale_y_continuous(expand=c(0,0.1)) + 
  scale_x_continuous(expand=c(0,0.1)) +
  theme(legend.position = "none", 
        plot.margin = margin(0,0,0,0, unit="cm"), 
        plot.title = element_text(hjust = 0.5, vjust=5, face="italic", 
                                  color="grey40", size=12)) #+  ggtitle(label="Shared Fungal ASVs")
  

ggsave("./data/graphs/Fig1B_18venndiagram.tiff",
       width = 7,
       height = 7,
       units = "cm",
       dpi = 1000 )



