library(phyloseq)
library(tidyverse)
library(plyr)



physeq <- readRDS("./data/Data/physeq18SFungi.rds")
physeq16S <- readRDS("./data/Data/physeq18SFungi.rds")

AscomycotaFAM <- c("Metschnikowiaceae", "Dipodascaceae", "Aureobasidiaceae", 
                   "Arthopyreniaceae", "Pleosporaceae", "Torulaceae", 
                   "Tubeufiaceae", "Venturiaceae", "Botryosphaeriaceae", 
                   "Phyllostictaceae", "Eremomycetaceae", "Aspergillaceae", 
                   "Onygenaceae", "Ajellomycetaceae", "Gymnoascaceae", "Arthrodermataceae", 
                   "Trichocomaceae", "Extremaceae", "Stigmatodiscaceae", 
                   "Cladosporiaceae", "Leotiaceae", "Xylariaceae", "Sporocadaceae", 
                   "Apiosporaceae", "Ophiostomataceae", "Cryphonectriaceae", 
                   "Valsaceae", "Phomatosporaceae", "Plectosphaerellaceae", 
                   "Pleurotheciaceae", "Clavicipitaceae", "Cordycipitaceae", 
                   "Nectriaceae", "Glomerellaceae", "Microascaceae", 
                   "Halosphaeriaceae", "Graphiaceae", "Stachybotryaceae", 
                   "Lulworthiaceae", "Magnaporthaceae", "Coniochaetaceae", 
                   "Boliniaceae", "Chaetomiaceae", "Lasiosphaeriaceae", 
                   "Helotiaceae", "Sympoventuriaceae", "Sordariaceae", 
                   "Erysiphaceae", "Trichomonascaceae", "Taphrinaceae", 
                   "Pyxidiophoraceae", "Ceratostomataceae", "Orbiliaceae", 
                   "Herpotrichiellaceae", "Trichomeriaceae", "Pyronemataceae", 
                   "Ascobolaceae", "Rhizinaceae", "Lipomycetaceae", 
                   "Saccharomycetaceae", "Pichiaceae", "Debaryomycetaceae", 
                   "Phaffomycetaceae", "Saccharomycodaceae", "Saccharomycopsidaceae", 
                   "Archaeorhizomycetaceae")
BasidiomycotaFAM <- c("Tulasnellaceae", "Didymellaceae", "Schizophyllaceae", 
                      "Psathyrellaceae", "Hygrophoraceae", "Omphalotaceae", 
                      "Physalacriaceae", "Tricholomataceae", "Marasmiaceae", 
                      "Lyophyllaceae", "Amanitaceae", "Agaricaceae", 
                      "Strophariaceae", "Cortinariaceae", "Atheliaceae", 
                      "Geastraceae", "Suillaceae", "Coriolaceae", "Corticiaceae", 
                      "Chionosphaeraceae", "Agaricostilbaceae", "Fereydouniaceae", 
                      "Cryptobasidiaceae", "Entylomataceae", "Tilletiaceae", 
                      "Geminibasidiaceae", "Tritirachiaceae", "Cystobasidiaceae",
                      "Pucciniaceae", "Sporidiobolaceae", "Leucosporidiaceae", 
                      "Camptobasidiaceae", "Chrysozymaceae", "Filobasidiaceae", 
                      "Piskurozymaceae", "Cryptococcaceae", "Sirobasidiaceae", 
                      "Bulleribasidiaceae", "Trichosporonaceae", "Tetragoniomycetaceae", 
                      "Rhynchogastremataceae", "Bulleraceae", "Exobasidiaceae", 
                      "Brachybasidiaceae", "Mrakiaceae", "Cystofilobasidiaceae", 
                      "Auriculariaceae", "Stereaceae", "Russulaceae", "Dacrymycetaceae", 
                      "Hydnodontaceae", "Ceratobasidiaceae", "Serpulaceae", "Moniliellaceae", 
                      "Malasseziaceae", "Ustilaginaceae")
MucoromycotaFAM <- c("Claroideoglomus")
GlomeromycotaFAM <- c("Rhizophagus")
IncartaeSedisFAM <- c("Papulosaceae", "Incertae Sedis")
BlastocladiomycotaFAM  <- c("Physoderma")
AscomycotaORD <- c("Sordariales", "Agaricales", "Chaetothyriales")
BasidiomycotaORD <- c("Agaricales", "Spiculogloeales", "Cystofilobasidiales", 
                      "Auriculariales", "Russulales", "Thelephorales", 
                      "Polyporales")
BlastocladiomycotaORD <- c("Coelomomycetaceae", "Physodermataceae")
ZoopagomycotaORD <- c("Sigmoideomycetaceae", "Sigmoideomycetaceae")
ChytridiomycotaORD <- c("Synchytriaceae", "Gromochytriaceae", "Nowakowskiellaceae", 
                        "Chytridiaceae", "Monoblepharidaceae", "Spizellomycetaceae", 
                        "Lobulomycetaceae", "Rhizophydiaceae", "Chytriomycetaceae")
IncartaeSedisORD <- c("Incertae Sedis", "Powellomycetaceae", "Olpidiaceae")
NeocallimastigomycotaORD <- c("Neocallimastigaceae")
MucoromycotaORD <- c("Endogonaceae", "Mortierellaceae", "Umbelopsidaceae", 
                     "Cunninghamellaceae", "Lichtheimiaceae", "Backusellaceae", 
                     "Rhizopodaceae", "Mucoraceae")

taxtable <- data.frame(tax_table(physeq)) %>%
  dplyr::mutate(Phylum = if_else(Phylum %in% c("Pezizomycotina", "Pucciniomycotina"), 
                          "Basidiomycota", Phylum), 
         Phylum = str_replace(Phylum, "Glomeromycetes", "Glomeromycota"), 
         Phylum = str_replace(Phylum, "Chytridiomycetes", "Chytridiomycota"), 
         Phylum = if_else(Class=="Entomophthorales", "Entomophthoromycota", Phylum), 
         Phylum = if_else(Class=="Glomerales", "Glomeromycota", Phylum), 
         Phylum = if_else(Class=="Zoopagales", "Zygomycota", Phylum), 
         Phylum = if_else(Class %in% c("Dothideomycetes", "Eurotiomycetes", 
                                       "Arthoniomycetes",
                                       "Ascomycota",
                                       "Sordariomycetes", 
                                       "Leotiomycetes", 
                                       "Saccharomycetes"), "Ascomycota", Phylum), 
         Phylum = if_else(Class=="Agaricomycetes", "Basidiomycota", Phylum), 
         Phylum = if_else(Class %in% c("Blastocladiales", "Blastocladiomycetes"),
                          "Blastocladiomycota", Phylum), 
         Phylum = if_else(Class %in% c("Lobulomycetales","Chytridiales", 
                                       "Rhizophydiales", "Gromochytriales"), "Chytridiomycota", Phylum), 
         Phylum = if_else(Class %in% c("Mortierellales","Mucorales"), "Mucoromycota", Phylum), 
         Phylum = if_else(Family %in% AscomycotaFAM, "Ascomycota", Phylum), 
         Phylum = if_else(Family %in% BasidiomycotaFAM, "Basidiomycota", Phylum), 
         Phylum = if_else(Family %in% MucoromycotaFAM, "Mucoromycota", Phylum), 
         Phylum = if_else(Family %in% GlomeromycotaFAM, "Glomeromycota", Phylum), 
         Phylum = if_else(Family %in% IncartaeSedisFAM, "Incartae sedis", Phylum), 
         Phylum = if_else(Family %in% BlastocladiomycotaFAM, "Blastocladiomycota", Phylum), 
         Phylum = if_else(Order %in% AscomycotaORD, "Ascomycota", Phylum), 
         Phylum = if_else(Order %in% BasidiomycotaORD, "Basidiomycota", Phylum), 
         Phylum = if_else(Order %in% BlastocladiomycotaORD, "Blastocladiomycota", Phylum), 
         Phylum = if_else(Order %in% ChytridiomycotaORD, "Chytridiomycota", Phylum), 
         Phylum = if_else(Order %in% IncartaeSedisORD, "Incartae sedis", Phylum), 
         Phylum = if_else(Order %in% NeocallimastigomycotaORD, "Neocallimastigomycota", Phylum), 
         Phylum = if_else(Order %in% MucoromycotaORD, "Mucoromycota", Phylum), 
         Phylum = if_else(Order %in% ZoopagomycotaORD, "Zoopagomycota", Phylum),
         Phylum = if_else(is.na(Phylum), "Unidentified", Phylum),
         Confidence = as.numeric(Confidence))

table(taxtable$Phylum)
tax_table(physeq) <- taxtable %>% as.matrix()


### CLEAN METADATA
metadatadfraw <- read.csv("./data/physeqraw/cold seep env.csv")
metadatadf <- read.csv("./data/physeqraw/cold seep env.csv") %>%
  dplyr::rename(Methane = CH4.mg.kg., 
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
  select(1:6)
sampledf <- read.delim("./data/physeqraw/sample-metadata.txt")


sampledfclean <-
  sampledf %>%
  column_to_rownames(., var="X") %>%
  dplyr::rename(Coretype = Core.type) %>%
  dplyr::mutate(Depth = gsub("cm", "", Depth), 
         Depth = as.numeric(Depth), 
         Habitat = gsub("_mat", "", Habitat), 
         Habitat = gsub("_", "", Habitat)) %>% 
  dplyr::mutate_if(., is.character, as.factor)  %>%
  dplyr::mutate(SampleIdentifier = paste0(ROV, Depth)) %>%
  dplyr::mutate(Depth = plyr::round_any(Depth, 5, f = round), 
                SampleID = rownames(.))


meta1 <- metadatadf %>%
  dplyr::rename(Phosphate=colnames(.)[2], 
         Sulfide=colnames(.)[3], 
         Sulfate=colnames(.)[4], 
         Ammonium=colnames(.)[5]) %>%
  dplyr::mutate(Sulfide=as.numeric(Sulfide)) %>%
  column_to_rownames(var="Samples") %>%
  rownames_to_column(var="SampleID") %>%
  dplyr::mutate(Depth = as.numeric(str_match(SampleID, "\\d+$")), 
         ROV = as.factor(str_match(SampleID, "ROV[O0](\\d)")[,2]), 
         ROV = paste0("ROV", ROV))  %>%
  dplyr::mutate(Depth = plyr::round_any(Depth, 5, f = round)) %>%
  dplyr::mutate(SampleIdentifier = paste0(ROV, Depth))

meta2 <- metadatadfraw[,8:10] %>%
  dplyr::rename(SampleID = Samples.1, 
         C13Methane = `C13..CH4..AOM.`, 
         d15Nitrogen = `d15N...`) %>%
  dplyr::mutate(SampleID = trimws(SampleID)) %>%
  drop_na() %>%
  dplyr::mutate(Depth = as.numeric(str_match(SampleID, "\\d+$")), 
         ROV = as.factor(str_match(SampleID, "ROV[O0](\\d)")[,2]), 
         ROV = paste0("ROV", ROV)) %>% group_by(ROV)  %>%
  dplyr::mutate(Depth = plyr::round_any(Depth, 5, f = round)) %>%
  dplyr::mutate(SampleIdentifier = paste0(ROV, Depth))


meta3 <- metadatadfraw[,12:14] %>%
  dplyr::rename(SampleID = Samples.2, 
         d13Carbon = `d13CV.PDB..`, 
         DIC = `DIC.C..mmol.L.`) %>%
  drop_na() %>%
  dplyr::mutate(Depth = as.numeric(str_match(SampleID, "\\d+$")), 
         ROV = as.factor(str_match(SampleID, "ROV[O0](\\d)")[,2]), 
         ROV = paste0("ROV", ROV)) %>% group_by(ROV)  %>%
  dplyr::mutate(Depth = plyr::round_any(Depth, 5, f = round)) %>%
  dplyr::mutate(SampleIdentifier = paste0(ROV, Depth))


# MATCH SAMPLE NAMES
sampledfcleaner <- sampledfclean %>% 
  mutate(Phosphate=NA, 
         Sulfide=NA, 
         Sulfate=NA, 
         Ammonium=NA, 
         Methane=NA, 
         C13Methane=NA, 
         d15Nitrogen=NA, 
         d13Carbon=NA, 
         DIC=NA)
for (i in 1:nrow(sampledfcleaner)){
  identifier <- sampledfcleaner[i,"SampleIdentifier"]
  
  target1 <- meta1 %>% 
    filter(SampleIdentifier==identifier) %>%
    select(Phosphate, Sulfide, Sulfate, Ammonium, Methane)
  target1avg <- colMeans(target1, na.rm=TRUE)
  sampledfcleaner[i, 9:13] <- target1avg
  
  target2 <- meta2 %>% 
    ungroup() %>%
    filter(SampleIdentifier==identifier) %>%
    select(C13Methane, d15Nitrogen)
  target2avg <- colMeans(target2, na.rm=TRUE)
  sampledfcleaner[i, 14:15] <- target2avg
  
  target3 <- meta3 %>% 
    ungroup() %>%
    filter(SampleIdentifier==identifier) %>%
    select(d13Carbon, DIC)
  target3avg <- colMeans(target3, na.rm=TRUE)
  sampledfcleaner[i, 16:17] <- target3avg
  
  
}
sample_data(physeq) <- sampledfcleaner
sample_data(physeq16S) <- sampledfcleaner

saveRDS(physeq, "./data/Data/physeq18SFungi.rds")
saveRDS(physeq16S, "./data/Data/physeq16S.rds")

