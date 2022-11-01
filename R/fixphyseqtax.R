library(phyloseq)
library(tidyverse)


physeq <- readRDS("./data/Data/physeq18SFungi.rds")


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
  mutate(Phylum = if_else(Phylum %in% c("Pezizomycotina", "Pucciniomycotina"), 
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

saveRDS(physeq, "./data/Data/physeq18SFungi.rds")
