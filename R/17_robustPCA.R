library(phyloseq)
library(SpiecEasi)
library(tidyverse)
library(ggpubr)
library(corrplot)



### LOAD DATA
BAFslr <- readRDS("./data/networks/raw/BAF/BAF_SpiecEasi_network_nlambda20_lambda.min.ratio0.005_rep.num5_ncores4_stability0.01707_SLR.rds")
BAslr <- readRDS("./data/networks/raw/BA/BA_SpiecEasi_network_nlambda20_lambda.min.ratio0.005_rep.num5_ncores4_stability0.01607_SLR.rds")
Bslr <- readRDS("./data/networks/raw/BA/BA_SpiecEasi_network_nlambda20_lambda.min.ratio0.005_rep.num5_ncores4_stability0.01607_SLR.rds")

BAFphyseq <-  readRDS("./data/Data/physeqcombined.rds")
BAphyseq <-  readRDS("./data/Data/physeq16S_seep_NetREADY_taxprevl_0.5.rds")[[1]]
# Bphyseq <- readRDS("./data/Data/physeq18SFungi_seep_NetREADY_taxprevl_0.2.rds")[[1]]

sampdata <- data.frame(sample_data(BAFphyseq))


### PERFORM ROBUST PCA

# Bacterial-Archaeal-Fungal network
XBAF <- BAFslr$est$data
LBAF <- BAFslr$est$resid[[getOptInd(BAFslr)]]
rownames(LBAF) <- colnames(XBAF)
colnames(LBAF) <- colnames(XBAF)


Lsvd <- svd(LBAF)
ind <- Lsvd$d > 1e-09
loadings <- diag(sqrt(1/Lsvd$d[ind])) %*% t(Lsvd$v[, 
                                                     ind])
varexplained <- (Lsvd$d^2/sum(Lsvd$d^2))[Lsvd$d > 1e-09]
plot(Lsvd$d^2/sum(Lsvd$d^2), xlim = c(0, length(Lsvd$d)), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")
scores <- XBAF %*% t(loadings)
pcaBAF <- list(scores = scores, loadings = loadings, 
               varexpl = varexplained)

plotdfBAF <- data.frame(sampdata, 
                     `Sequencing Depth` = rowSums(otu_table(BAFphyseq)/1000), 
                     PCA = pcaBAF$scores) %>%
  rename_all(function(x) str_replace_all(x, "A\\.", "")) %>%
  rename(`Sampling Depth (m)` = Depth) %>%
  mutate(Methane= Methane/1000)

plotdfnumBAF <- plotdfBAF %>%
  select_if(is.numeric) %>%
  select(-PC13, -PC12) %>%
  drop_na()

corresBAF <- data.frame(cor(plotdfnumBAF)) %>%
  select_if(grepl("PC", names(plotdfnumBAF))) %>%
  filter(!grepl("PC", rownames(.))) %>%
  as.matrix()


corrplot(corresBAF, type = "full",method = 'circle',
         tl.col = "black", tl.srt = 45, number.cex=0.5, tl.cex=0.5,
         p.mat = corresBAF, insig = "blank")



# Bacterial-Archaeal network
XBA <- BAslr$est$data
LBA <- BAslr$est$resid[[getOptInd(BAslr)]]
rownames(LBA) <- colnames(XBA)
colnames(LBA) <- colnames(XBA)


Lsvd <- svd(LBA)
ind <- Lsvd$d > 1e-09
loadings <- diag(sqrt(1/Lsvd$d[ind])) %*% t(Lsvd$v[, 
                                                   ind])
varexplained <- (Lsvd$d^2/sum(Lsvd$d^2))[Lsvd$d > 1e-09]
plot(Lsvd$d^2/sum(Lsvd$d^2), xlim = c(0, length(Lsvd$d)), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")
scores <- XBA %*% t(loadings)
pcaBA <- list(scores = scores, loadings = loadings, 
               varexpl = varexplained)

plotdfBA <- data.frame(sampdata, 
                        `Sequencing Depth` = rowSums(otu_table(BAphyseq)/1000), 
                        PCA = pcaBA$scores) %>%
  rename_all(function(x) str_replace_all(x, "A\\.", "")) %>%
  rename(`Sampling Depth (m)` = Depth) %>%
  mutate(Methane= Methane/1000)

plotdfnumBA <- plotdfBA %>%
  select_if(is.numeric) %>%
  select(-PC13, -PC12) %>%
  drop_na()

corresBA <- data.frame(cor(plotdfnumBA)) %>%
  select_if(grepl("PC", names(plotdfnumBA))) %>%
  filter(!grepl("PC", rownames(.))) %>%
  as.matrix()


corrplot(corresBA, type = "full",method = 'circle',
         tl.col = "black", tl.srt = 45, number.cex=0.5, tl.cex=0.5,
         p.mat = corresBA, insig = "blank")




#### PLOT 

plotdfBAF %>%
  ggplot(aes(x=PC1, y = PC2)) +
  geom_point(aes(color=ROV)) +
  # scale_color_brewer(palette="Set1") +
  theme_bw()


plotdfBAF %>%
  ggplot(aes(x=PC1, y = `Sampling Depth (m)` )) +
  geom_point(size=1, color="grey") +
  geom_smooth(formula = y~x, method="glm", color="black", se=FALSE) +
  stat_cor(aes(label=..rr.label..)) +
  stat_cor(aes(label=..p.label..),
           p.accuracy =  0.01,
           vjust=3) +
  theme_bw() +
  ylab("Sequencing Depth") +
  xlab("PC1")



# logistic regression model
testdf <- plotdfBA %>%
  mutate(MethLev = if_else(Methane<20, 0, 1))%>%
  select(MethLev, PC3)


model <- glm(data=testdf, formula=MethLev~PC3, family = "binomial") 
R2 <- with(summary(model), 1 - deviance/null.deviance)
R2

plotdfLR <- testdf %>% select(PC3)
plotdfLR %>%
  mutate(MethLev = predict(model, plotdfLR, type="response")) %>%
  ggplot(aes(x=PC3, y=MethLev)) +
  geom_point(data=testdf, aes(x=PC3, y=MethLev, color=MethLev)) + 
  stat_smooth(method="glm", se=FALSE, method.args = list(family=binomial)) + 
  scale_x_reverse() + 
  theme_bw()




plotdfBAF %>%
  mutate(MethLev = if_else(Methane<70, "Low", "High"), 
         SulfateLev = if_else(Sulfate<1000, "Low", "High")) %>%
  ggplot(aes(x=SulfateLev, y=PC3)) +
  geom_point(aes(color=ROV), size=1, height = 10) +
  geom_smooth(formula = y~x, method="glm", color="black", se=FALSE) +
  stat_compare_means(method = "anova") + 
  theme_bw()

