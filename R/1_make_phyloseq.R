library(phyloseq)
library(tidyverse)
library(data.table)
library(ape)

##### 16S OTU #####

#### load OTU, Tree, and Taxonomy data ####
otu <- read.delim2("./data/physeqraw/coldseep_18S_ASV_protists.txt") %>%
  column_to_rownames(var = "X") %>%
  t() %>%
  as.data.frame()


x<- c(D_0__="D0_", 
      D_1__="D1_", 
      D_2__="D2_", 
      D_3__="D3_", 
      D_4__="D4_", 
      D_5__="D5_", 
      D_6__="D6_", 
      D_7__="D7_", 
      D_8__="D8_",
      D_9__="D9_",
      D_10__="D10_", 
      D_11__="D11_")

tax <- read.delim2("./data/physeqraw/taxonom_18s_all.txt") %>%
  column_to_rownames(var = "Feature_id") %>%
  filter(grepl("Fungi", Taxon)) %>%
  mutate(Taxon = str_replace_all(Taxon, c(D_0__="D0_", 
                                          D_1__="D1_", 
                                          D_2__="D2_", 
                                          D_3__="D3_", 
                                          D_4__="D4_", 
                                          D_5__="D5_", 
                                          D_6__="D6_", 
                                          D_7__="D7_", 
                                          D_8__="D8_",
                                          D_9__="D9_",
                                          D_10__="D10_", 
                                          D_11__="D11_"))) %>%
  separate(Taxon, 
           remove = FALSE,
           c("D0", 
             "D1", 
             "D2", 
             "D3", 
             "D4", 
             "D5", 
             "D6", 
             "D7", 
             "D8", 
             "D9", 
             "D10", 
             "D11"), 
           sep = ";")


for (irow in 1:nrow(tax)){
  tax[irow, "D0"] <- str_match(tax[irow, "Taxon"], "D0_(.+?);")[,2]
  tax[irow, "D1"] <- str_match(tax[irow, "Taxon"], "D1_(.+?);")[,2]
  tax[irow, "D2"] <- str_match(tax[irow, "Taxon"], "D2_(.+?);")[,2]
  tax[irow, "D3"] <- str_match(tax[irow, "Taxon"], "D3_(.+?);")[,2]
  tax[irow, "D4"] <- str_match(tax[irow, "Taxon"], "D4_(.+?);")[,2]
  tax[irow, "D5"] <- str_match(tax[irow, "Taxon"], "D5_(.+?);")[,2]
  tax[irow, "D6"] <- str_match(tax[irow, "Taxon"], "D6_(.+?);")[,2]
  tax[irow, "D7"] <- str_match(tax[irow, "Taxon"], "D7_(.+?);")[,2]
  tax[irow, "D8"] <- str_match(tax[irow, "Taxon"], "D8_(.+?);")[,2]
  tax[irow, "D9"] <- str_match(tax[irow, "Taxon"], "D9_(.+?);")[,2]
  tax[irow, "D10"] <- str_match(tax[irow, "Taxon"], "D10_(.+?);")[,2]
  tax[irow, "D11"] <- str_match(tax[irow, "Taxon"], "D11_(.+)")[,2]
}

taxclean <- tax %>%
  select(-Taxon, -D5, -D4) %>%
  rename(Domain = D0, 
         Supergroup = D1, 
         Superclade = D2, 
         Kingdom = D3, 
         Phylum = D6, 
         Class = D7, 
         Order = D8, 
         Family = D9,
         Genus = D10, 
         Species = D11)



tree <- read.tree("./data/physeqraw/rooted-tree-18s_protists.nwk")  

#### LOAD & merge environmental & sample data ####
envdat <- read.csv("./data/physeqraw/cold seep env.csv") %>%
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
  mutate(SampleID = paste("D", Samples, sep="") %>%
           str_replace_all("-", "") %>%
           str_replace("ROV0", "ROV"), 
         SampleIDnew = SampleID) %>% 
  column_to_rownames(var="SampleID") %>%
  mutate(ROV  = str_extract(SampleIDnew, "ROV."), 
         Depth = as.numeric(str_extract(SampleIDnew, "\\d?\\d$")))

sampledat <- read.delim2("./data/physeqraw/sample-metadata.txt") %>%
  mutate(SampleID = X %>%
           str_replace_all("\\.", ""), 
         Depth = as.numeric(Depth %>% str_extract_all("\\d+"))) %>% 
  column_to_rownames(var = "SampleID")

namesenv <- rownames(envdat) %>% data.frame()
namessamp <- rownames(sampledat) %>% data.frame()




evndatafilt <- envdat %>%
  filter(ROV %in% c("ROV1", "ROV2", "ROV3", "ROV5"), 
         Depth %in% c(0, 5, 10, 15, 20, 25, 30))

# manually change 0cm to correct form 
evndatafilt["D517ROV2PC2030", "Depth"] <- 0
evndatafilt["D519ROV32PC2010", "Depth"] <- 0
evndatafilt["D519ROV32LPC10", "Depth"] <- 0
evndatafilt["D519ROV32PC2030", "Depth"] <- 0
evndatafilt["D522ROV12PC2030", "Depth"] <- 0
evndatafilt["D522ROV12LPC10", "Depth"] <- 0


# merge repetitive samples
evndatafilt <- evndatafilt %>%
  mutate(ROVDepth = paste(ROV, Depth, sep = ''))



#### PLOT 
evndatafilt %>%
  pivot_longer(cols = c(Methane, Phosphate, Sulfide, Sulfate, 
                        Ammonium, d15N, d13CV, DIC, C13Methane), 
               names_to = "Data", 
               values_to = "Concentration") %>% 
  group_by(ROVDepth, Data) %>%
  mutate(meanConcentration = mean(Concentration, na.rm = TRUE)) %>%
  ggplot(aes(x = meanConcentration, y = Depth, color = ROV)) + 
  geom_point() + 
  facet_grid(~Data, scales = "free_x") + 
  theme_bw()

envdatfinal <- evndatafilt %>%
  pivot_longer(cols = c(Methane, Phosphate, Sulfide, Sulfate, 
                        Ammonium, d15N, d13CV, DIC, C13Methane), 
               names_to = "Data", 
               values_to = "Concentration") %>% 
  group_by(ROVDepth, Data) %>%
  mutate(meanConcentration = mean(Concentration, na.rm = TRUE)) %>%
  distinct(Concentration, .keep_all = TRUE) %>%
  select(-Concentration, -X, -X.1, -X.2, -X.3) %>%
  pivot_wider(names_from = Data, 
              values_from = meanConcentration) %>%
  distinct(ROVDepth, .keep_all = TRUE)


sampledatfinal <- merge(sampledat, envdatfinal) %>%
  column_to_rownames(var = "X")



#### FILTER otu table samples to match metadata ####
otufilt <- otu[rownames(sampledatfinal),]
otufilt <- otufilt[,rownames(taxclean)]
taxclean <- taxclean[colnames(otufilt), ]
treeclean <- keep.tip(tree, colnames(otufilt))

#### MAKE PHYLOSEQ #####
physeq <- phyloseq(otu_table((otufilt %>% as.matrix()), 
                             taxa_are_rows = FALSE), 
                   tax_table(taxclean %>% as.matrix()), 
                   sample_data(sampledatfinal), 
                   phy_tree(treeclean))

saveRDS(physeq, "./data/physeq18SFungi.rds")
