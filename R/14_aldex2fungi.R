library(ALDEx2)
library(tidyverse)
library(ggrepel)

#### LOAD DATA ####
physeq18SFun <- readRDS("./data/Data/physeq16S.rds")
# physeqROV <- subset_samples(physeq18SFun, ROV=="ROV3" | ROV=="ROV5")

otu18SFun <- data.frame(otu_table(physeq18SFun))
sampdata <- data.frame(sample_data(physeq18SFun)) %>%
  mutate(Seep = if_else(ROV == "ROV5", "Control", "Seep"))
taxadf <- data.frame(tax_table(physeq18SFun))

conditions <- sampdata$Seep
otu18SFunt <- data.frame(t(otu18SFun))



#### RUN ALDEX 2 ####

# Note that running aldex2 on default setting with this dataset faces the 
# issue of assymetry (look at bioconductor documentation) due to the sample size
# representation being very different between seep and non seep

# for this reason we run each sample site against ROV5 independently 
aldexdf <- aldex(otu18SFunt, conditions, mc.samples=100, test="t", effect=TRUE,
                 include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
aldex.plot(aldexdf, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(aldexdf, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")


aldexdfsg <- aldexdf %>%
  rownames_to_column(var = "ASVid") %>%
  mutate(ASVid = ASVid %>% str_replace_all("^X", ""))%>%
  filter(we.ep <= 0.05)  %>%
  left_join(., (taxadf %>% rownames_to_column(var = "ASVid"))) %>%
  column_to_rownames(var ="ASVid")

aldexdfnsg <- aldexdf %>%
  filter(we.ep > 0.05)

aldexdfsg %>%
  ggplot(aes(x = effect, y = we.ep)) + 
  geom_point(size = 0.5, color = "red") + 
  geom_label_repel(aes(label = Species), size = 2) + 
  geom_point(size = 0.5, data = aldexdfnsg, aes(x = effect, y = we.ep), color = "grey") + 
  geom_hline(yintercept = 0.05, color = "blue", linetype="dashed") + 
  scale_y_reverse() + 
  theme_bw()



