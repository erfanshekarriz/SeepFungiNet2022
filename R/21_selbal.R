library(selbal)
library(phyloseq)
library(tidyverse)
library(metagMisc)

set.seed(1234)

#### LOAD DATA ####
physeq <- readRDS("./data/Data/physeq18SFungi.rds")
physeq.filt <- phyloseq_filter_prevalence(physeq, 
                           prev.trh = 0.1, 
                           abund.trh = NULL,
                           threshold_condition = "OR", 
                           abund.type = "total")
taxdf <- data.frame(tax_table(physeq.filt))
otumat <- as.matrix(data.frame(otu_table(physeq.filt)))
sampledf <- data.frame(sample_data(physeq.filt)) %>% mutate(Methane = Methane/1000)
methane <- as.numeric(sampledf[,"Methane"])
ggplot() + geom_histogram(aes(methane), bins=10, fill="blue", color="white") + theme_bw()


avgDepth0Meth <- sampledf %>%
  filter(Depth==0) %>%
  pull(Methane) %>%
  max()

methlev <- sampledf %>% 
  mutate(MethLev = if_else(Methane<=avgDepth0Meth, "Low", "High")) %>%
  pull(MethLev) %>% as.factor()

ROVcov <- sampledf %>% dplyr::mutate(ROV=as.factor(ROV)) %>% dplyr::select(ROV)


#### SELBAL MODELLING WITH DISCRETICIZED METHANE LEVELS ####

selbalmod <- selbal.cv(x = otumat, 
                       y = methlev, 
                       n.fold = 10, 
                       n.iter = 100,
                       logit.acc = "AUC")



#### VISUALIZE THE BEST FIT MODEL ####
selbalmod$accuracy.nvar
grid.draw(selbalmod$global.plot)
# plot.tab(selbalmod$cv.tab)
selbalmod$cv.accuracy; summary(selbalmod$cv.accuracy)
selbalmod$var.barplot
selbalmod$glm
selbalmod$opt.nvar

globalbal <- selbalmod$global.balance %>%
  mutate(Taxa=gsub("^X", "", Taxa)) %>%
  merge(., taxdf %>% rownames_to_column("Taxa"), by="Taxa") %>%
  column_to_rownames(var="Taxa")



#### CLEAN PLOTTING & SAVING ####
grid.draw(selbalmod$global.plot2)







#### SELBAL MODEL SELECTION WITH METHANE CONC ####  
selbalmod2 <- selbal.cv(x = otumat, 
                        y = methane,
                        covar = ROVcov, 
                        n.fold = 5, 
                        n.iter = 100,
                        logit.acc = "AUC")

selbalmod2$accuracy.nvar
grid.draw(selbalmod2$global.plot)
# plot.tab(selbalmod2$cv.tab)
selbalmod2$cv.accuracy; summary(selbalmod2$cv.accuracy)
selbalmod2$var.barplot
selbalmod2$glm
selbalmod2$opt.nvar

globalbal2 <- selbalmod2$global.balance %>%
  mutate(Taxa=gsub("^X", "", Taxa)) %>%
  merge(., taxdf %>% rownames_to_column("Taxa"), by="Taxa") %>%
  column_to_rownames(var="Taxa")




