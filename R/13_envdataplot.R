library(plyr)
library(tidyverse)
library(dplyr)
library(phyloseq)

#### LOAD & CLEAN DATA ####
physeq <- readRDS("./data/Data/physeq18SFungi.rds")
envdatall <- data.frame(sample_data(physeq))

envdatall %>%
  select_if(is.numeric) %>%
  select(-Depth) %>%
  cbind(ROV=envdatall$ROV) %>%
  pivot_longer(is.numeric) %>%
  ggplot(aes(x=value)) +
  geom_histogram(bins=10, fill="white", aes(color=name, y=..ndensity..)) +
  facet_grid(ROV~name, scales="free_x") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.text.y=element_text(size = 7, face ="italic"), 
        axis.text.x=element_text(size = 7, face ="italic", angle=90), 
        axis.title = element_blank()) + 
  scale_color_brewer(palette="Set1")

ggsave("./data/graphs/supplementary/FigS3_15envdataplot.png",
       width = 25,
       height = 15,
       units = "cm",
       dpi = 300 )


envdatall %>%
  select_if(is.numeric) %>%
  select(-Depth) %>%
  cbind(ROV=envdatall$ROV) %>%
  pivot_longer(is.numeric) %>%
  ggplot(aes(x=value)) +
  geom_histogram(bins=15, fill="white", aes(color=name, y=..ndensity..)) +
  facet_grid(~name, scales="free") + 
  theme_classic() + 
  theme(legend.position = "none", 
        axis.text.y=element_text(size = 7, face ="italic"), 
        axis.text.x=element_text(size = 7, face ="italic", angle=90), 
        axis.title = element_blank()) + 
  scale_color_brewer(palette="Set1")

ggsave("./data/graphs/supplementary/FigS4_15envdataplot.png",
       width = 29,
       height = 6,
       units = "cm",
       dpi = 300 )



#### PLOT 


averagedf <- envdatall %>%
  drop_na() %>%
  dplyr::group_by(ROV) %>%
  dplyr::mutate(Methane= Methane/(9)) %>%
  dplyr::mutate(Methane = mean(Methane), 
         Sulfate = mean(Sulfate), 
         Sulfide = mean(Sulfide, na.rm=TRUE), 
         Ammonium = mean(Ammonium, na.rm=TRUE), 
         Phosphate = mean(Phosphate), 
         C13Methane = mean(C13Methane, na.rm=TRUE), 
         d15Nitrogen = mean(d15Nitrogen, na.rm=TRUE), 
         d13Carbon = mean(d13Carbon, na.rm=TRUE), 
         DIC = mean(DIC, na.rm=TRUE))%>%
  pivot_longer(colnames(envdatall)[9:17], 
               names_to = "Nutrient", values_to = "Concentration") 
  
averagedf %>%
  filter(Nutrient %in% c("Methane", "Sulfate", "Sulfide")) %>%
  dplyr::mutate(Nutrient = factor(Nutrient, levels=c("Methane", "Sulfate", "Sulfide"))) %>%
  ggplot(aes(x=reorder(ROV, Concentration), y=Concentration)) +
  geom_bar(aes(fill=ROV), width = 0.3, position="dodge", stat="identity") + 
  theme_bw() + 
  facet_wrap(~Nutrient, scales = "free_y", ncol = 1) +
  theme(strip.background=element_rect(size=0.6, fill="white"), 
        strip.text.x=element_text(face = "bold.italic", size = 7), 
        axis.title.x=element_text(face="italic", size=8, vjust=-3), 
        axis.title.y=element_text(face="italic", size=8, vjust=2), 
        axis.text.y=element_text(size = 7, face ="italic"), 
        axis.text.x=element_text(size = 7, face ="italic", angle=90), 
        legend.text=element_text(size = 6, face="italic"), 
        legend.position = "none",
        panel.border=element_rect(colour="black", size = 1), 
        panel.background=element_rect(fill="grey92"), 
        plot.background=element_blank(), 
        panel.grid.minor = element_line(color="grey98"),
        panel.grid.major = element_line(color="grey98"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines"))  +
  ylab("Concentration (mg/L)") +  
  xlab("Sampling Site")  + 
  guides(fill=guide_legend(title="Nutrient")) + 
  scale_fill_manual(values=c("blue2", "red2", "green", "grey"))


ggsave("./data/graphs/Fig1_15envdataplot.tiff",
       width = 5.5,
       height = 11,
       units = "cm",
       dpi = 1000 )

averagedf %>%
  ggplot(aes(x=reorder(ROV, Concentration), y=Concentration)) +
  geom_bar(aes(fill=ROV), width = 0.08, position="dodge", stat="identity") + 
  theme_bw() + 
  facet_wrap(~Nutrient, scales = "free_y", ncol = 3) +
  theme(strip.background=element_rect(size=0.6, fill="white"), 
        strip.text.x=element_text(face = "bold.italic", size = 7), 
        axis.title.x=element_text(face="italic", size=8, vjust=-3), 
        axis.title.y=element_text(face="italic", size=8, vjust=2), 
        axis.text.y=element_text(size = 7, face ="italic"), 
        axis.text.x=element_text(size = 7, face ="italic", angle=90), 
        legend.text=element_text(size = 6, face="italic"), 
        legend.position = "none",
        panel.border=element_rect(colour="black", size = 1), 
        panel.background=element_rect(fill="grey92"), 
        plot.background=element_blank(), 
        panel.grid.minor = element_line(color="grey98"),
        panel.grid.major = element_line(color="grey98"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines"))  +
  ylab("Concentration (mg/L)") +  
  xlab("Sampling Site")  + 
  guides(fill=guide_legend(title="Nutrient")) + 
  scale_fill_manual(values=c("blue2", "red2", "green", "grey"))


ggsave("./data/graphs/supplementary/FigS1_15envdataplot.png",
       width = 14,
       height = 14,
       units = "cm",
       dpi = 300 )



sumdepthROV <- 
  envdatall %>%
  drop_na() %>%
  dplyr::mutate(Methane= Methane/(9)) %>%
  dplyr::group_by(ROV, Depth) %>%
  dplyr::mutate(Methane = mean(Methane), 
                Sulfate = mean(Sulfate), 
                Sulfide = mean(Sulfide, na.rm=TRUE), 
                Ammonium = mean(Ammonium, na.rm=TRUE), 
                Phosphate = mean(Phosphate), 
                C13Methane = mean(C13Methane, na.rm=TRUE), 
                d15Nitrogen = mean(d15Nitrogen, na.rm=TRUE), 
                d13Carbon = mean(d13Carbon, na.rm=TRUE), 
                DIC = mean(DIC, na.rm=TRUE))%>%
  pivot_longer(colnames(envdatall)[9:17], 
               names_to = "Nutrient", values_to = "Concentration") %>%
  dplyr::filter(Depth %in% c(0,5,10,15,20,25,30))


depthdf <- envdatall %>%
  drop_na() %>%
  dplyr::mutate(Methane= Methane/(9)) %>%
  dplyr::filter(Depth %in% c(0,5,10,15,20,25,30)) %>%
  dplyr::group_by(ROV, Depth) %>%
  pivot_longer(colnames(envdatall)[9:17], 
               names_to = "Nutrient", values_to = "Concentration") 

sumdepthROVselect <- sumdepthROV %>%
  dplyr::filter(Nutrient %in% c("Methane", "Sulfate", "Sulfide")) %>%
  dplyr::mutate(Nutrient = factor(Nutrient, levels=c("Methane", "Sulfate", "Sulfide"))) 
depthdf %>%
  dplyr::filter(Nutrient %in% c("Methane", "Sulfate", "Sulfide")) %>%
  dplyr::mutate(Nutrient = factor(Nutrient, levels=c("Methane", "Sulfate", "Sulfide"))) %>%
  ggplot(aes(y=Concentration, x = Depth)) + 
  geom_point(data=sumdepthROVselect, aes(x=Depth, y=Concentration, color=ROV)) + 
  geom_smooth(data=sumdepthROVselect, aes(x=Depth, y=Concentration, color=ROV), se=FALSE, size=0.5) + 
  facet_wrap(~Nutrient, scales = "free_x", nrow=1) + 
  coord_flip() +   
  theme_bw() + 
  theme(strip.background=element_rect(size=0.6, fill="white"), 
        strip.text.x=element_text(face = "bold.italic", size = 7), 
        axis.title.x=element_text(face="bold", size=8, vjust=-3), 
        axis.title.y=element_text(face="bold", size=8, vjust=4), 
        axis.text.y=element_text(size = 7, face ="italic"), 
        axis.text.x=element_text(size = 7, face ="italic", angle=90), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_blank(),
        panel.border=element_rect(colour="black", size = 1), 
        panel.background=element_rect(fill="white"), 
        plot.background=element_blank(), 
        panel.grid.minor = element_line(color="grey98"),
        panel.grid.major = element_line(color="grey98"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(1, "lines"))  +
  ylab("Average Concentration (mg/L)") +  
  xlab("DBSF (cm)")  + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values=c("blue2", "red2", "green2", "grey")) + 
  geom_vline(xintercept = 13, color="orange", linetype="dashed") + 
  geom_text(x = 18, y=10, label = "SMTZ", color="red", size=5)


ggsave("./data/graphs/Fig1C_15envdataplot.tiff",
       width = 13,
       height = 9,
       units = "cm",
       dpi = 1000 )


depthdf %>%
  ggplot(aes(y=Concentration, x = Depth)) + 
  geom_point(data=sumdepthROV, aes(x=Depth, y=Concentration, color=ROV)) + 
  geom_smooth(data=sumdepthROV, aes(x=Depth, y=Concentration, color=ROV), se=FALSE, size=0.5) + 
  facet_wrap(~Nutrient, scales = "free_x", nrow=3) + 
  coord_flip() +   
  theme_bw() + 
  theme(strip.background=element_rect(size=0.6, fill="white"), 
        strip.text.x=element_text(face = "bold.italic", size = 7), 
        axis.title.x=element_text(face="bold.italic", size=8, vjust=-3), 
        axis.title.y=element_text(face="bold.italic", size=8, vjust=4), 
        axis.text.y=element_text(size = 7, face ="italic"), 
        axis.text.x=element_text(size = 7, face ="italic", angle=90), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_blank(),
        panel.border=element_rect(colour="black", size = 1), 
        panel.background=element_rect(fill="white"), 
        plot.background=element_blank(), 
        panel.grid.minor = element_line(color="grey98"),
        panel.grid.major = element_line(color="grey98"),
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(1, "lines"))  +
  ylab("Average Concentration (mg/L)") +  
  xlab("DBSF (cm)")  + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values=c("blue2", "red2", "green2", "grey")) + 
  geom_vline(xintercept = 13, color="orange", linetype="dashed") + 
  geom_text(x = 18, y=10, label = "SMTZ", color="red", size=5)


  

ggsave("./data/graphs/supplementary/FigS2_15envdataplot.png",
       width = 15,
       height = 20,
       units = "cm",
       dpi = 300 )





