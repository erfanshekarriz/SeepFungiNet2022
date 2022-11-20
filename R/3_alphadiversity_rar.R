library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)
library(iNEXT)
library(corrplot)
library(PerformanceAnalytics)
library(GGally)


### the iNEXT rarefaction curves take a very long time to perform, 
### so we have commented out analysis needing lots of time. 

### If you wish to recreate everything from scratch feel free to manually un-comment!


set.seed(1234)


#### LOAD DATA ####
physeq16S <- readRDS('./data/Data/physeq16S.rds')
physeq18S <- readRDS('./data/Data/physeq18SFungi.rds')
physeq18SROV <- merge_samples(physeq18S, "ROV")
physeq16SROV <- merge_samples(physeq16S, "ROV")

# otu16S <- data.frame(otu_table(physeq16S))
# otu18S <- data.frame(otu_table(physeq18S))
# otu16SROV <- data.frame(otu_table(physeq16SROV))
# otu18SROV <- data.frame(otu_table(physeq18SROV))


tax16S <- data.frame(tax_table(physeq16S))
tax18S <- data.frame(tax_table(physeq18S))

metadata <- data.frame(sample_data(physeq16S))

metadata %>%
  dplyr::group_by(ROV) %>%
  dplyr::summarize(n=n())

# metadata 
sampledataplotraw <- data.frame(sample_data(physeq16S)) %>%
  as.data.frame() %>%
  mutate(SampleID = rownames(.))
features <- colnames(sampledataplotraw) %>% as.data.frame()


#### ROV BASED RAREFACITON RICHNESS CURVES ####

### 16S 
# inext16S <- iNEXT(t(otu16SROV), 
#                   q=c(0,1,2),# hill number diversity order
#                   datatype="abundance", 
#                   nboot=100, 
#                   se = FALSE, 
#                   conf = 0.95)
# 
# saveRDS(inext16S, "./data/Data/ROVrarecurveiNEXT16S.rds")
inext16S <- readRDS("./data/Data/ROVrarecurveiNEXT16S.rds")

plotdf16S <- inext16S$iNextEst$size_based %>% mutate(m = m/10000, 
                                                     qD = qD/1000)  %>%
  filter(Order.q == 0)
observed16S <- plotdf16S %>% filter(Method == "Observed")
extrapolate16S <- plotdf16S %>% filter(Method == "Extrapolation")



### 18S
# inext18S <- iNEXT(t(otu18SROV), 
#                   q=c(0,1,2),# hill number diversity order
#                   datatype="abundance", 
#                   nboot=100, 
#                   se = TRUE, 
#                   conf = 0.95)
# saveRDS(inext18S, "./data/Data/ROVrarecurveiNEXT18S.rds")
inext18S <- readRDS("./data/Data/ROVrarecurveiNEXT18S.rds")

plotdf18S <- inext18S$iNextEst$size_based %>% mutate(m = m/10000, 
                                                     qD = qD/1000)  %>%
  filter(Order.q == 0)

observed18S <- plotdf18S %>% filter(Method == "Observed") 
extrapolate18S <- plotdf18S %>% filter(Method == "Extrapolation")

#### PLOT ROV BASED COMBINED ####
plotdf <- rbind(cbind(plotdf16S, Group="16S"), cbind(plotdf18S, Group="18S Fungi")) %>%
  mutate(Group = factor(Group, levels=c("18S Fungi", "16S")))
plotext <- rbind(cbind(extrapolate16S, Group="16S"), cbind(extrapolate18S, Group="18S Fungi")) %>%
  mutate(Group = factor(Group, levels=c("18S Fungi", "16S")))
plotobs <- rbind(cbind(observed16S, Group="16S"), cbind(observed18S, Group="18S Fungi")) %>%
  mutate(Group = factor(Group, levels=c("18S Fungi", "16S")))


scales_y <- list(
  `16S` = scale_y_continuous(limits = c(0, 27)),
  `18S Fungi` = scale_y_continuous(limits = c(0, 1.5))
)

plotdf %>%
  filter(Method == "Rarefaction") %>%
  mutate(Group = factor(Group, levels=c("18S Fungi", "16S"))) %>%
  ggplot(aes(x = m, y = qD, group = Assemblage)) + 
  geom_line(color = "grey", size = 0.40, linetype="solid") + 
  geom_line(data = plotext, aes(x = m, y=qD, group=Assemblage, color = Assemblage), 
            size = 0.90, linetype="dashed") + 
  geom_point(data = plotobs, aes(x = m, y = qD), size =1, 
             color = "blue") + 
  geom_text(data = plotobs, aes(x = m, y = qD, label = Assemblage, color = Assemblage), 
            size = 4, hjust=-1.5, vjust=-1, fontface="bold") + 
  geom_text(data = plotobs, aes(x = m, y = qD), 
            label = "OBSERVED", size = 2, vjust=-1) + 
  facet_wrap(~Group, scales = "free", nrow = 2, strip.position="top") + 
  theme_bw() + 
  theme(strip.background=element_rect(fill="grey95", color="white"), 
        strip.text=element_text(face = "bold.italic", size = 12), 
        axis.title.x=element_text(face="bold", size=11, vjust=-3), 
        axis.title.y=element_text(face="bold", size=11, vjust=4), 
        axis.text=element_text(size = 9, face ="bold"), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_text(size=6, face="bold"),
        legend.position="none",
        panel.border=element_rect(colour="black", size = 1), 
        panel.background=element_blank(), 
        plot.background=element_blank(), 
        plot.margin=unit(c(0.5,1,1,1), "cm"),
        panel.grid=element_line(color="grey95"), 
        panel.grid.minor=element_line(color="grey95"))  + 
  # expand_limits(y=c(0,27)) +
  scale_color_brewer(palette="Set1") + 
  xlab("Sequences Sampled (x10,000)") + 
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15)), 
                     position="left") + 
  ylab("Richness (x1000)") 


ggsave("./data/graphs/Fig2E_3alphadiversity.tiff",
       width = 9,
       height = 16,
       units = "cm",
       dpi = 1000 )




#### SAMPLE BASED COVERAGE CURVES ####

### 16S 
# inextSample16S <- iNEXT(t(otu16S), 
#                   q=c(0,1,2),# hill number diversity order
#                   datatype="abundance", 
#                   nboot=100, 
#                   se = FALSE, 
#                   conf = 0.95)
# 
# saveRDS(inextSample16S, "./data/Data/SAMPrarecurveiNEXT16S.rds")
inextSample16S <- readRDS("./data/Data/SAMPrarecurveiNEXT16S.rds")

plotdf16Ssamp <- inextSample16S$iNextEst$size_based %>% mutate(m = m/10000, 
                                                     qD = qD/1000)  %>% filter(Order.q == 0) %>%
  left_join(., metadata %>% rownames_to_column(var = "Assemblage"), by = "Assemblage" )
observed16Ssamp <- plotdf16Ssamp %>% filter(Method == "Observed")
extrapolate16Ssamp <- plotdf16Ssamp %>% filter(Method == "Extrapolation")

plotdf16Ssamp %>%
  filter(Method == "Rarefaction") %>%
  ggplot(aes(x = m, y = qD, group = Assemblage)) + 
  geom_line(color = "grey", size = 0.10, linetype="solid") + 
  geom_line(data = extrapolate16Ssamp, aes(x = m, y=qD, group=Assemblage, color=ROV), 
            size = 0.50, linetype="solid", color="red") + 
  geom_point(data = observed16Ssamp, aes(x = m, y = qD), size =1, 
             color = "blue") + 
  theme_bw() + 
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_text(size = 12, vjust= 4), 
        plot.title = element_text(size = 12, face = "italic"),
        legend.text =   element_text(size = 9), 
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey95"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        legend.position = c(0.87,0.84),
        legend.title = element_blank()) +
  expand_limits(x=c(0,22)) + 
  scale_color_brewer(palette = "Set1") + 
  xlab("Sequence Sampled (x10,000)") + 
  ylab("Richness (x1000)") + 
  ggtitle("Prokaryotic 16S Diversity")

ggsave("./data/graphs/supplementary/FigS1_3alphadiversity.png",
       width = 15,
       height = 13,
       units = "cm",
       dpi = 300 )



plotdf16Ssamp %>%
  filter(Method == "Rarefaction") %>%
  ggplot(aes(x = m, y = SC, group = Assemblage)) + 
  geom_line(color = "grey", size = 0.10, linetype="solid") + 
  geom_line(data = extrapolate16Ssamp, aes(x = m, y=SC, group=Assemblage), 
            size = 0.50, linetype="solid", color = "red") + 
  geom_point(data = observed16Ssamp, aes(x = m, y = SC), size =1, 
             color = "blue") + 
  theme_bw() + 
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_text(size = 12, vjust= 4), 
        plot.title = element_text(size = 12, face = "italic"),
        legend.text =   element_text(size = 9), 
        # legend.title =   element_text(size = 9),
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey95"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        legend.position = c(0.87,0.84),
        legend.title = element_blank()) + 
  xlab("Sequence Sampled (x10,000)") + 
  ylab("Proportion ASV Coverage") + 
  ggtitle("Prokaryotic 16S Diversity")

ggsave("./data/graphs/supplementary/FigS2_3alphadiversity.png",
       width = 15,
       height = 13,
       units = "cm",
       dpi = 300 )


### 18S
# inextSample18S <- iNEXT(t(otu18S), 
#                         q=c(0,1,2),# hill number diversity order
#                         datatype="abundance", 
#                         nboot=100, 
#                         se = FALSE, 
#                         conf = 0.95)
# 
# saveRDS(inextSample18S, "./data/Data/SAMPrarecurveiNEXT18S.rds")
inextSample18S <- readRDS("./data/Data/SAMPrarecurveiNEXT18S.rds")

plotdf18Ssamp <- inextSample18S$iNextEst$size_based %>% mutate(m = m/10000, 
                                                               qD = qD/1000)  %>% 
  filter(Order.q == 0) %>%
  left_join(., metadata %>% rownames_to_column(var = "Assemblage"), by = "Assemblage" )
  
observed18Ssamp <- plotdf18Ssamp %>% filter(Method == "Observed")
extrapolate18Ssamp <- plotdf18Ssamp %>% filter(Method == "Extrapolation")

plotdf18Ssamp %>%
  filter(Method == "Rarefaction") %>%
  ggplot(aes(x = m, y = qD, group = Assemblage, color=ROV)) + 
  geom_line(color = "grey", size = 0.10, linetype="solid") + 
  geom_line(data = extrapolate18Ssamp, aes(x = m, y=qD, group=Assemblage, color = ROV), 
            size = 0.50, linetype="solid", color="red") + 
  geom_point(data = observed18Ssamp, aes(x = m, y = qD), size =1, 
             color = "blue") + 
  theme_bw() + 
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_text(size = 12, vjust= 4), 
        plot.title = element_text(size = 12, face = "italic"),
        legend.text =   element_text(size = 9), 
        # legend.title =   element_text(size = 9),
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey95"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        legend.position = c(0.87,0.84),
        legend.title = element_blank()) + 
  scale_color_brewer(palette="Set1") +
  expand_limits(x=c(0,15)) + 
  xlab("Sequence Sampled (x10,000)") + 
  ylab("Richness (x1000)") + 
  ggtitle("Fungal 18S Diversity")

ggsave("./data/graphs/supplementary/FigS3_3alphadiversity.png",
       width = 15,
       height = 13,
       units = "cm",
       dpi = 300 )


plotdf18Ssamp %>%
  filter(Method == "Rarefaction") %>%
  ggplot(aes(x = m, y = SC, group = Assemblage)) + 
  geom_line(color = "grey", size = 0.10, linetype="solid") + 
  geom_line(data = extrapolate18Ssamp, aes(x = m, y=SC, group=Assemblage), 
            size = 0.50, linetype="solid", color = "red") + 
  geom_point(data = observed18Ssamp, aes(x = m, y = SC), size =1, 
             color = "blue") + 
  theme_bw() + 
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_text(size = 12, vjust= 4), 
        plot.title = element_text(size = 12, face = "italic"),
        legend.text =   element_text(size = 9), 
        legend.title =   element_text(size = 9),
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey95"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        legend.position = "none") + 
  xlab("Sequence Sampled (x10,000)") + 
  ylab("Proportion ASV Coverage") + 
  ggtitle("Fungal 18S Diversity")

ggsave("./data/graphs/supplementary/FigS4_3alphadiversity.png",
       width = 15,
       height = 13,
       units = "cm",
       dpi = 300 )



#### ALPHA DIVERSITY CORRELATION ####
numericdf <- observed18Ssamp %>%
  select_if(is.numeric) %>%
  select_if(~sum(!is.na(.x)) > 0) %>%
  select(-Order.q, -m, -SC) %>%
  dplyr::rename(Richness=qD) %>%
  drop_na()

corrdf <- cor(numericdf) 
myfunc <- function(i,j) mapply(function(a,b) cor.test(mtcars[[a]], mtcars[[b]])$p.value, i, j)
pvaldf <- outer(numericdf, numericdf, Vectorize(function(a, b) cor.test(a, b)$p.value)) 


# only label isgnificant 
png("./data/graphs/supplementary/FigS4_3alphadiversity.png",
       width = 10,
       height = 10,
       units = "cm",
     res = 300 )

corrplot(corrdf, type = "lower",method = 'number',
         tl.col = "black", tl.srt = 45, number.cex=0.5, tl.cex=0.5,
         p.mat = pvaldf, insig = "pch") 
dev.off()


ggpairs(numericdf) +  theme_classic()



# 16S
numericdf16S <- observed16Ssamp %>%
  select_if(is.numeric) %>%
  select_if(~sum(!is.na(.x)) > 0) %>%
  select(-Order.q, -m, -SC) %>%
  dplyr::rename(Richness=qD) %>%
  drop_na()

corrdf16S <- cor(numericdf16S) 
pvaldf16S <- outer(numericdf16S, numericdf16S, Vectorize(function(a, b) cor.test(a, b)$p.value)) 


# only label isgnificant 
png("./data/graphs/supplementary/FigS5_3alphadiversity.png",
    width = 10,
    height = 10,
    units = "cm",
    res = 300 )

corrplot(corrdf16S, type = "lower",method = 'number',
         tl.col = "black", tl.srt = 45, number.cex=0.5, tl.cex=0.5,
         p.mat = pvaldf16S, insig = "pch") 
dev.off()


ggpairs(numericdf16S) +  theme_classic()



