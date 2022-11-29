library(tidyverse)

df <- read.csv("./data/Data/fungalhitsmetatrans.csv")


scale <- 5
levels <- c("Total Gene Count", "CYP450s", "Hydrophobins", "Ligninolytic Enzymes")
df %>% 
  pivot_longer(-Sample) %>%
  mutate(name = str_replace_all(name, c(cyp450="CYP450s", 
                                    fungseq="Total Gene Count", 
                                    hydrophobin="Hydrophobins", 
                                    ligninolytic="Ligninolytic Enzymes")), 
         value=as.integer(value)) %>%
  mutate(name = factor(name, levels = levels)) %>%
  ggplot(aes(x=value, y=Sample, fill=name)) + 
  geom_bar(stat = "identity", 
           position="dodge", 
           color="black", 
           width=0.3, alpha=0.8) + 
  xlab("\nUnique Fungal BLAST Hits (>97% Identitiy)") + 
  facet_wrap(~name, nrow=1, scales="free_x") + 
  ylab("Metatranscriptomic Samples\n") + 
  scale_color_brewer(palette="Set1") + 
  scale_fill_brewer(palette="Set1") + 
  theme(strip.background=element_rect(size=0.6, fill="grey93", color="grey93"), 
        strip.text.x=element_text(face = "bold.italic", size = 7),
        axis.title = element_text(face="italic", size=9), 
        axis.text =  element_text(size = 9), 
        axis.text.y =  element_text(size = 9, angle=0, face="bold"), 
        legend.text = element_text(size = 6, face="bold.italic"), 
        legend.background = element_blank(),
        legend.position = "none",
        plot.margin = margin(10,10,10,10),
        legend.title = element_blank(),
        panel.border =element_rect(colour="black", fill=NA, size = 1.2), 
        panel.background = element_rect(fill="white"), 
        plot.background = element_blank(), 
        panel.grid = element_line(color = "grey95")) + 
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))


ggsave("./data/graphs/Fig6B_24metatransres.tiff", 
       units="cm", 
       height=5, 
       width = 15, 
       dpi=1000)
