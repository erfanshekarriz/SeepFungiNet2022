library(phyloseq)
library(tidyverse)
library(RColorBrewer)

set.seed(1234)

#### ALL EUKARYOTES ####
taxonALL <- read_delim("./data/physeqraw/taxonom_18s_all.txt") %>%
  filter(!grepl("D_1__Archaeplastida", Taxon)) %>%
  column_to_rownames(var=colnames(.)[1]) %>%
  mutate(Fungi= if_else(grepl("Fungi", Taxon), "Fungi", Taxon), 
         Fungi= if_else((Fungi=="D_0__Eukaryota"), "Unidentified", Fungi), 
         Fungi= if_else((Fungi!="Unidentified" & Fungi!="Fungi"), "Other", Fungi)) %>%
  mutate(Clade= str_match(Taxon, "D_2__(.+?);")[,2]) %>%
  mutate(Clade = na_if(Clade, "uncultured marine eukaryote"), 
         Clade = na_if(Clade, "uncultured phototrophic eukaryote"), 
         Clade = na_if(Clade, "metagenome"), 
         Clade = na_if(Clade, "uncultured"), 
         Clade = na_if(Clade, "uncultured Sarcosomataceae"), 
         Clade = na_if(Clade, "uncultured dinoflagellate"), 
         Clade = na_if(Clade, "uncultured rhodophyte"), 
         Clade = if_else((Fungi=="Fungi"), "Fungi", Clade), 
         Clade = if_else(is.na(Clade), "Unidentified", Clade))

unique(taxonALL$Clade)
otu18SFULL <- read_delim("./data/physeqraw/coldseep_18s_ASV_all eukaryotes.txt") %>%
  column_to_rownames(var=colnames(.)[1]) %>% 
  filter(rownames(.) %in% rownames(taxonALL)) %>%
  apply(., 2, function(x) x/sum(x)) %>%
  as.data.frame()



sampdata <- data.frame(sample_data(readRDS("./data/Data/physeq16S.rds")))

datamerge <- otu18SFULL %>%
  rownames_to_column(var="ASV") %>%
  merge(taxonALL %>%
          rownames_to_column(var="ASV"), 
        by="ASV") 


datagglom <- 
  datamerge %>%
  group_by(Fungi) %>%
  summarise(across(is.numeric, sum)) %>%
  ungroup() %>%
  column_to_rownames(var=colnames(.)[1]) %>%
  select(-Confidence) %>%
  t() %>% as.data.frame()



plotdfALL <- datagglom %>%
  rownames_to_column(var="SampleID") %>%
  merge(sampdata , 
        by = "SampleID") %>%
  group_by(ROV, Depth) %>%
  mutate(meanFung = mean(Fungi), 
         meanOth = mean(Other), 
         meanUni = mean(Unidentified)) %>%
  ungroup() %>%
  pivot_longer(c(meanFung, meanOth, meanUni), 
               names_to = "Domain", 
               values_to = "meanAbundance") %>%
  distinct(meanAbundance, .keep_all=TRUE) %>%
  add_column(Group="EUKARYOTES")

#### ALL PROTISTS ####
otu18SPROT <- read_delim("./data/physeqraw/coldseep_18S_ASV_protists.txt") %>%
  column_to_rownames(var=colnames(.)[1]) %>% 
  filter(rownames(.) %in% rownames(taxonALL)) %>%
  apply(., 2, function(x) x/sum(x)) %>%
  as.data.frame()

datamergePROT <- otu18SPROT %>%
  rownames_to_column(var="ASV") %>%
  merge(taxonALL %>%
          rownames_to_column(var="ASV"), 
        by="ASV")

unique(datamergePROT$Clade)

datagglomPROT <- datamergePROT %>%
  group_by(Fungi) %>%
  summarise(across(is.numeric, sum)) %>%
  ungroup() %>%
  column_to_rownames(var=colnames(.)[1]) %>%
  select(-Confidence) %>%
  t() %>% as.data.frame()

cladglomPROT <- datamergePROT %>%
  group_by(Clade) %>%
  summarise(across(is.numeric, sum)) %>%
  ungroup() %>%
  column_to_rownames(var=colnames(.)[1]) %>%
  select(-Confidence) %>%
  t() %>% as.data.frame()


plotdfPROT <- datagglomPROT %>%
  rownames_to_column(var="SampleID") %>%
  merge(sampdata,
        by = "SampleID") %>%
  group_by(ROV, Depth) %>%
  mutate(meanFung = mean(Fungi), 
         meanOth = mean(Other), 
         meanUni = mean(Unidentified)) %>%
  ungroup() %>%
  pivot_longer(c(meanFung, meanOth, meanUni), 
               names_to = "Domain", 
               values_to = "meanAbundance") %>%
  distinct(meanAbundance, .keep_all=TRUE) %>%
  add_column(Group="PROTISTS")

plotCLADEprot <- cladglomPROT %>%
  rownames_to_column(var="SampleIDX") %>%
  merge(sampdata %>%
          rownames_to_column(var="SampleIDX"), 
        by = "SampleIDX") %>%
  group_by(ROV, Depth) %>%
  summarise(across(is.numeric, mean)) %>%
  ungroup() %>%
  pivot_longer(unique(datamergePROT$Clade), 
               names_to = "Domain", 
               values_to = "meanAbundance") %>%
  distinct(meanAbundance, .keep_all=TRUE) %>%
  add_column(Group="PROTISTS")



#### PLOT ####
plotdf <- rbind(plotdfALL, plotdfPROT)%>%
  mutate(Kingdom = if_else(Domain=="meanFung", "Fungi", Domain), 
         Kingdom = if_else(Kingdom=="meanUni", "Unidentified", Kingdom), 
         Kingdom = if_else(Kingdom=="meanOth", "Other", Kingdom)) 

summarytab <- plotCLADEprot %>%
  group_by(Domain) %>%
  summarise(mean = mean(meanAbundance), 
            sd = sd(meanAbundance))

plotdf %>%
  filter(Group=="PROTISTS") %>%
  ggplot(aes(x=as.factor(Depth), y=meanAbundance, fill=Domain, group=ROV)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_y_continuous(expand = c(0, 0), position = "right") + 
  scale_fill_brewer(palette="Set1") +
  facet_grid(~ROV, scales = "free_x", switch="y") + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(face = "bold", size = 10), 
        axis.title.x=element_text(face="bold", size=8, vjust=-3), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 7, face ="bold"), 
        axis.text.x=element_text(size = 7, face ="bold"), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_text(face="bold"),
        legend.position = "right",
        panel.border=element_rect(colour="grey", size = 0), 
        panel.background=element_blank(), 
        plot.background=element_blank(), 
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid=element_blank(), 
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines"))  +
  xlab("Depth (cm)") +  
  ylab("Relative Abundance")  + 
  guides(fill=guide_legend(title="Mean Relative \nAbundance \n(18S Micreukaryote)\n"))


ggsave("./data/graphs/Fig2D_13percentfungi.tiff",
       width = 16,
       height = 6,
       units = "cm",
       dpi = 1000 )





plotCLADEprot %>%
  mutate(Domain= if_else((meanAbundance < 0.1 & Domain!="Unidentified"), 
                         "Other", Domain)) %>%
  mutate(Domain = factor(Domain, levels=c("Fungi", "Alveolata",
                                          "Rhizaria", "Stramenopiles",
                                          "Other", "Unidentified"))) %>%
  ggplot(aes(x=as.factor(Depth), y=meanAbundance, fill=Domain, group=ROV)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_y_continuous(expand = c(0, 0), position = "right") + 
  scale_fill_brewer(palette="Set3") +
  facet_grid(~ROV, scales = "free_x", switch="y") + 
  theme_bw() + 
  theme(strip.background=element_blank(), 
        strip.text.x=element_text(face = "bold", size = 10), 
        axis.title.x=element_text(face="bold", size=8, vjust=-3), 
        axis.title.y=element_blank(), 
        axis.text.y=element_text(size = 7, face ="bold"), 
        axis.text.x=element_text(size = 7, face ="bold"), 
        legend.text=element_text(size = 6, face="bold"), 
        legend.background=element_blank(),
        legend.title=element_blank(),
        legend.position = "right",
        panel.border=element_rect(colour="grey", size = 0), 
        panel.background=element_blank(), 
        plot.background=element_blank(), 
        plot.margin=unit(c(1,1,1,1), "cm"),
        panel.grid=element_blank(), 
        panel.spacing.y = unit(0.4, "lines"), 
        panel.spacing.x = unit(0.05, "lines"))  +
  xlab("Depth (cm)") +  
  ylab("Relative Abundance")  + 
  guides(fill=guide_legend(title="Mean Relative \nAbundance \n(18S Micreukaryote)\n"))


ggsave("./data/graphs/supplementary/FigS1_13percentfungi.png",
       width = 16,
       height = 6,
       units = "cm",
       dpi = 300 )





### CORELATION STUDY 


  



