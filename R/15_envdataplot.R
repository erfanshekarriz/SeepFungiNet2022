library(plyr)
library(tidyverse)
library(dplyr)

#### LOAD & CLEAN DATA ####
envdatall <- read.csv("./data/physeqraw/cold seep env.csv") %>%
  rename(Methane = CH4.mg.kg., 
         Phosphate = PO4.ppm., 
         Sulfide = sulfide_mg.L, 
         Sulfate = SO42...mg.L., 
         Ammonium = NH4...mg.L., 
         C13Methane = C13..CH4..AOM., 
         d15N = d15N..., 
         d13CV = d13CV.PDB.., 
         DIC = DIC.C..mmol.L.) %>%
  mutate(across(.cols = c(Methane, Phosphate, Sulfide, Sulfate, Ammonium, d15N,
                          d13CV, DIC), 
                ~ as.numeric(.x))) %>%
  select(-X, -X.1, -X.2, -X.3)
  # mutate(SampleID = paste("D", Samples, sep="") %>%
  #          str_replace_all("-", "") %>%
  #          str_replace("ROV0", "ROV"), 
  #        SampleIDnew = SampleID) %>% 
  # column_to_rownames(var="SampleID") %>%
  # mutate(ROV  = str_extract(SampleIDnew, "ROV."), 
  #        Depth = as.numeric(str_extract(SampleIDnew, "\\d?\\d$")))

envdat1 <- envdatall[,1:6] %>%
  mutate(Depth = as.numeric(str_extract(Samples, "\\d?\\d$")),
         Depth = plyr::round_any(Depth, 5, f = round),
         SampleID = paste("D", Samples, sep="") %>%
           str_replace_all("-", "") %>%
           str_replace("ROV0", "ROV"),
         SampleIDnew = SampleID) %>%
  column_to_rownames(var="SampleID") %>%
  mutate(ROV  = str_extract(SampleIDnew, "ROV."),
         ROV = str_replace_all(ROV, "O\\d", ""))
envdat2 <- envdatall[,7:9] %>%
  rename(SampleID = Samples.1) %>%
  mutate(SampleID = trimws(SampleID), 
         ROV  = str_extract(SampleID, "ROV[O0]."),
         ROV = str_replace_all(ROV, "(O)\\d", "1"), 
         ROV = str_replace_all(ROV, "0", ""),
         Depth = as.numeric(str_extract(SampleID, "\\d?\\d$")),
         Depth = plyr::round_any(Depth, 5, f = round)) %>% na.omit()


envdat3 <- envdatall[,10:12] %>%
  rename(SampleID = Samples.2) %>%
  mutate(SampleID = trimws(SampleID), 
         ROV  = str_extract(SampleID, "ROV[O0]."),
         ROV = str_replace_all(ROV, "(O)\\d", "1"), 
         ROV = str_replace_all(ROV, "0", ""),
         Depth = as.numeric(str_extract(SampleID, "\\d?\\d$")),
         Depth = plyr::round_any(Depth, 5, f = round)) %>% na.omit()


#### PLOT 


envdat1 %>%
  drop_na() %>%
  group_by(ROV) %>%
  mutate(Methane= Methane/(9)) %>%
  mutate(Methane = mean(Methane), 
         Sulfate = mean(Sulfate), 
         Sulfide = mean(Sulfide, na.rm=TRUE), 
         Ammonium = mean(Ammonium, na.rm=TRUE), 
         Phosphate = mean(Phosphate)) %>%
  pivot_longer(c(Methane, Sulfate, Sulfide, Ammonium, Phosphate), 
               names_to = "Nutrient", values_to = "Concentration") %>%
  filter(Nutrient!="Ammonium", 
         Nutrient!="Phosphate") %>%
  mutate(Nutrient = factor(Nutrient, levels=c("Methane", "Sulfate", "Sulfide", 
                           "Phosphate", "Ammonium"))) %>%
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



sumdepthROV <- 
  envdat1 %>%
  drop_na() %>%
  mutate(Methane= Methane/(9)) %>%
  group_by(ROV, Depth) %>%
  summarise(Methane = mean(Methane), 
            Sulfate = mean(Sulfate), 
         Sulfide = mean(Sulfide, na.rm=TRUE), 
         Ammonium = mean(Ammonium, na.rm=TRUE), 
         Phosphate = mean(Phosphate))  %>% 
  pivot_longer(c(Methane, Sulfate, Sulfide, 
                 Ammonium, Phosphate), 
               names_to = "Nutrient", values_to = "Concentration") %>%
  filter(Depth %in% c(0,5,10,15,20,25,30)) %>% 
  mutate(Nutrient = factor(Nutrient, levels=c("Methane", "Sulfate", "Sulfide", 
                                              "Ammonium", "Phosphate")))  %>%
  filter(Nutrient!="Ammonium",
         Nutrient!="Phosphate")


envdat1 %>%
  drop_na() %>%
  mutate(Methane= Methane/(9)) %>%
  filter(Depth %in% c(0,5,10,15,20,25,30)) %>%
  group_by(ROV, Depth) %>%
  pivot_longer(c(Methane, Sulfate, Sulfide, Ammonium, Phosphate), 
               names_to = "Nutrient", values_to = "Concentration") %>%
  mutate(Nutrient = factor(Nutrient, levels=c("Methane", "Sulfate", "Sulfide", 
                                              "Ammonium", "Phosphate"))) %>%
  filter(Nutrient!="Ammonium",
         Nutrient!="Phosphate") %>%
  ggplot(aes(y=Concentration, x = Depth)) + 
  # geom_point(aes(x=Depth, y=Concentration, color=ROV), size=0.1) + 
  geom_point(data=sumdepthROV, aes(x=Depth, y=Concentration, color=ROV)) + 
  geom_smooth(data=sumdepthROV, aes(x=Depth, y=Concentration, color=ROV), se=FALSE, size=0.5) + 
  facet_wrap(~Nutrient, scales = "free_x", nrow=1) + 
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
  xlab("Depth (cm)")  + 
  scale_x_reverse() + 
  scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values=c("blue2", "red2", "green2", "grey")) + 
  geom_vline(xintercept = 18, color="orange", linetype="dashed") + 
  geom_text(x = 18, y=10, label = "SMTZ", color="red", size=5)


ggsave("./data/graphs/Fig1C_15envdataplot.tiff",
       width = 13,
       height = 9,
       units = "cm",
       dpi = 1000 )

  





