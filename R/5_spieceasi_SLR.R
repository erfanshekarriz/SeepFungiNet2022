library(phyloseq)
library(SpiecEasi)
library(tidyverse)
library(igraph)
library(ggpubr)
library(gridExtra)

# name: construct spieceasi-slr network on server
#
# input: input directory with filtered phyloseq RDS files
# output: output directory for saving network


                                                       
#### INPUT PARAMETERS ####
# methods <- c('slr', 'glasso', 'mb')
methods <- c('slr')
nlambda.x <- 20
lambda.min.ratio.x <- 0.005
thresh.x <- 0.05
rep.num.x <- 5
ncores.x <- 4

# the maximum integer value for β while it begins from 0 and stepwise increases
# the model with the lowest eBIC is selected

rankmin <- 2
rankmax <- 50       # Maximum equal to the number of your ASVs
                      # By default set to ASV-1 count if you leave this variable blank
ranksteps <- 20

#### SETUP FILE DIRECTORIES ####
inputdir <- getwd()
dir.create("./data/graphs", showWarnings = FALSE)
dir.create("./data/networks", showWarnings = FALSE)



#### READ PHYSEQ in ./DATA DIRECTORY ####
files <- list.files(path= './data/Data', 
                    pattern="^.*_NetREADY_.*$", 
                    full.names=TRUE, 
                    recursive=FALSE) 

physeqlist <- list()
for (file in files){
  physeqlist <- append(physeqlist, readRDS(file))
}
names(physeqlist) <- files %>%
  str_extract(pattern = "(physeq...)")


#### FILTER FOR COMMON SAMPLES ####
commonsamp <- unlist(unique(lapply(physeqlist, 
                                   FUN = function(x) sample_names(x))))
physeqlistcom <- list()
for(i in 1:length(physeqlist)) {        
  physeqlistcom[[i]] <- subset_samples(physeqlist[[i]], sample_names(physeqlist[[i]]) %in% commonsamp) 
  # physeqlistcom[[i]] <- subset_taxa(physeqlistcom[[i]],
  #                                   (Domain %in% c("Bacteria", "Eukaryota")))
}
names(physeqlistcom) <- names(physeqlist)


#### MAKE COMBINED PHYLOSEQ ####
combinedotu <- bind_cols(lapply(physeqlistcom, FUN = function(x) as.data.frame(otu_table(x))))
combinedtax <- bind_rows(lapply(physeqlistcom, FUN = function(x) as.data.frame(tax_table(x))))
metadata <- as.data.frame(sample_data(physeqlistcom[[1]]))

physeqcombined <- phyloseq(otu_table(combinedotu%>% as.matrix(), taxa_are_rows = FALSE), 
                           tax_table(combinedtax %>% as.matrix()),
                           sample_data(metadata))

saveRDS(physeqcombined, paste("./data","physeqcombined.rds", sep = '/'))

cat(paste("Combined phyloseq object made..."))



#### INFER SPIEC-EASI NETWORK ####
pulsar.params <- list(thresh = thresh.x, 
                      rep.num= rep.num.x, 
                      ncores = ncores.x, 
                      seed = 1234)

cat(paste("Started making the network at", Sys.time(), 'HKT'))

if (!exists("rankmax")){
  rankmax <- ncol(combinedotu) - 1
}

#### SET RANK for BETA Estimation 
ranks <- unique(round(exp(seq(log(rankmin), log(rankmax), len=ranksteps))))
cat(paste0("Selecting with the following ranks: ", paste(ranks, collapse = " ")))
SElist <- list()
for (i in 1:length(methods)){
  i <- 1
  netname <- paste("SE", methods[i], sep = "")
  cat(paste("\nStarted making the", methods[i],  "network at", Sys.time(), 'HKT\n'))

  # perform model selection for SLR β parameter
  if (methods[i]=="slr"){
    se.slr <- spiec.easi(physeqlistcom,
                         method= methods[i],
                         r = ranks, 
                         lambda.log=TRUE,
                         nlambda= nlambda.x,
                         lambda.min.ratio= lambda.min.ratio.x,
                         pulsar.params = pulsar.params)
    ebic <- sapply(se.slr, function(x)
      ebic(x$refit$stars, x$est$data,
           x$est$loglik[x$select$stars$opt.index]))
    ebic <- ebic[ebic!=0]
    se.slr$ebic <- ebic
    ebicminInd <- names(which.min(ebic))
    SElist[i] <- se.slr[ebicminInd]
    SElist[[i]]$ebictested <- ebic
    SElist[[i]]$eBIC <- min(ebic)
    
    }
  else{SElist[i] <- list(spiec.easi(physeqlistcom,
                                   method= methods[i],
                                   nlambda= nlambdatemp,
                                   lambda.log=TRUE,
                                   lambda.min.ratio= lambda.min.ratio.x,
                                   pulsar.params = pulsar.params))}

  stability <- signif(getStability(SElist[[i]]), digits = 4)
  names(SElist)[i] <- netname
  
  cat(paste("Finished making the", methods[i],  "network at", Sys.time(), 'HKT'))
  
  networkname <- paste("SpiecEasi_network_", 
                       "nlambda", nlambda.x,  "_",
                       'lambda.min.ratio', lambda.min.ratio.x, "_",
                       'rep.num', rep.num.x, "_",
                       'ncores', ncores.x, "_",
                       'stability', stability, 
                       '_', toupper(methods[i]),
                       sep = "")
  
  saveRDS(SElist[[i]], paste("./data/networks/", networkname, ".rds", sep = ''))
  
  
}

cat(paste("\nDONE! Finished making ALL the networks at", Sys.time(), 'HKT'))


cat(paste("\nNow calculating parameters and plotting graphs..."))





#### MANNUALLY UPLOAD SE IF NECESSSARY #### 
# SEfiles <- list.files(path = './data/networks/', 
#                       pattern = "^SpiecEasi.*$", 
#                       full.names = TRUE)
# SEfilesshort <- list.files(path = './data/networks/', 
#                       pattern = "^SpiecEasi.*$", 
#                       full.names = FALSE)
# SElist <- list()
# for (i in 1:length(SEfiles)){
#   SElist <- append(SElist, list(readRDS(SEfiles[i])))
#   names(SElist)[i] <- SEfilesshort[i]
#   
# }
# physeqlistcom <- list(readRDS("./data/physeqcombined.rds"))





#### Make igraph objects ####
igraphlist <- list()
igraphlistw <- list()

if("slr" %in% methods){
  SEslr <- SElist$SEslr
  icov <- SEslr$select$est$icov[[getOptInd(SEslr)]]
  cov <- solve(icov)
  
  igraph <- adj2igraph(getRefit(SEslr), vertex.attr=(list(name=taxa_names(physeqcombined))))
  igraphw <- graph.adjacency(cov*getRefit(SEslr), mode = "undirected", weighted = TRUE, diag = FALSE)
  igraphw <- set.vertex.attribute(igraphw, 
                                  "name", 
                                  value= as_ids(V(igraph)))
  
  igraphlist <- append(igraphlist, list(igraph))
  names(igraphlist)[length(igraphlist)] <- "igraphSLR"
  igraphlistw <- append(igraphlistw, list(igraphw))
  names(igraphlistw)[length(igraphlistw)] <- "igraphSLR"
}
if("mb" %in% methods){
  SEmb <- SElist$SEmb
  igraph <- adj2igraph(getRefit(SEmb), 
                       vertex.attr=(list(name=taxa_names(physeqcombined))))
  igraphw <- adj2igraph(symBeta(getOptBeta(SEmb), mode='maxabs'), 
                        vertex.attr=(list(name=taxa_names(physeqcombined))))
  
  igraphlist <- append(igraphlist, list(igraph))
  names(igraphlist)[length(igraphlist)] <- "igraphMB"
  igraphlistw <- append(igraphlistw, list(igraphw))
  names(igraphlistw)[length(igraphlistw)] <- "igraphMB"
}
if("glasso" %in% methods){
  SEglasso <- SElist$SEglasso
  igraph <- adj2igraph(getRefit(SEglasso), 
                       vertex.attr=(list(name=taxa_names(physeqcombined))))
  igraphw  <- adj2igraph(cov2cor(apply(getOptCov(SEglasso), 2, as.numeric)), 
                         vertex.attr=(list(name=taxa_names(physeqcombined))))
  igraphlist <- append(igraphlist, list(igraph))
  names(igraphlist)[length(igraphlist)] <- "igraphGLASSO"
  igraphlistw <- append(igraphlistw, list(igraphw))
  names(igraphlistw)[length(igraphlistw)] <- "igraphGLASSO"
}



#### Add taxonomy information ####
paste("Saving igraph objects...")
options(warn = -1)
for (i in 1:length(igraphlist)){
  for (j in 1:length(colnames(combinedtax))){
    igraphlist[[i]] <- set_vertex_attr(graph = igraphlist[[i]],
                                       name = colnames(combinedtax)[j], 
                                       value = combinedtax[, colnames(combinedtax)[j]]) 
    igraphlistw[[i]] <- set_vertex_attr(graph = igraphlistw[[i]],
                                        name = colnames(combinedtax)[j], 
                                        value = combinedtax[, colnames(combinedtax)[j]]) 
    saveRDS(igraphlist[[i]], paste("./data/networks/", names(igraphlist)[i], "_igraph.rds", sep = ''))
    saveRDS(igraphlistw[[i]], paste("./data/networks/", names(igraphlistw)[i], "_igraphw.rds", sep = ''))
    
  }
}
options(warn = 0)


#### CALCULATE GRAPH BASED PARAMETERS
paste("Calculating network parameters...")

#### CALCULATE GLOBAL PARAMETERS
param_mat<- matrix(nrow = length(igraphlist), ncol = 7)
rownames(param_mat) <- names(SElist)
colnames(param_mat) <- c("Stability", "Largest Connected Component", 
                         "Average Degree", "Average Distance", "Average Diameter", 
                         "Average Betweenness Centrality", "Transativity")
customtheme <- theme_bw() + 
  theme(plot.margin = margin(1,1,1,1, "cm"), 
        # panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(face = "italic"),
        axis.title.x = element_text(margin=margin(t=10), size = 9, face = "italic"),
        axis.title.y =  element_text(margin=margin(r=5), size = 8, face = "italic"),
        axis.text.x = element_text(size = 10, face = "italic"),
        axis.text.y = element_text(size = 10, margin=margin(r=5), face = "italic"),
        title = element_text(size = 8, margin=margin(r=5), face = "italic"), 
        legend.text = element_text(size = 7),
        legend.position = "right",
        legend.background =  element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8))

for (i in 1:length(igraphlist)){
  logcon <- ifelse(components(igraphlist[[i]])$no == 1, FALSE, TRUE)
  cc <- components(igraphlist[[i]])
  
  param_mat[i,1] <- getStability(SElist[[i]]) 
  param_mat[i,2] <- cc$csize[1]
  param_mat[i,3] <- mean(degree(igraphlist[[i]]))
  param_mat[i,4] <- mean_distance(igraphlist[[i]], directed = FALSE, unconnected = logcon) #nw
  param_mat[i,5] <- diameter(igraphlist[[i]], directed=F, unconnected = logcon) #nw
  param_mat[i,6] <- centr_betw(igraphlist[[i]], directed=F, normalized = T)$centralization #
  param_mat[i,7] <- transitivity(igraphlist[[i]], type="global")
  
  
  # connected components histogram
  cctab <- cc$csize %>% 
    table() %>%
    as.data.frame() %>%
    rename(., 
           "Size of Connected Component" = `.`, 
           "Count" = Freq)
  
  png(paste("./data/graphs/", names(igraphlist)[i], ".png", sep=""), height = 500*nrow(cctab), width = 2000*ncol(cctab), res=600)
  cctab %>%
    tableGrob() %>%
    grid.arrange(top = names(igraphlist)[i], .)
  
  dev.off()
  
  
  
  # modularity analysis 
  # mod <- modularity(igraphw)
  
  
  #### GRAPH LOCAL PARAMETERS ####
  
  # edge distribution
  edgedis <-  E(igraphlistw[[i]])$weight
  p1 <- edgedis %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   binwidth = (max(edgedis)-min(edgedis))/length(edgedis)*100, 
                   fill = "Blue") + 
    geom_density(aes(y = ..density..), color = "Red") + 
    customtheme + 
    ggtitle("Edge Weight Distribution") + 
    labs(x = "Density distribution", 
         y = "Edge Weight")
  
  # degrees
  degs <- degree(igraphlist[[i]])
  p2 <- degs %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    filter(X != 0) %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   binwidth = (max(degs)-min(degs))/length(degs)*10, 
                   fill = "Blue") + 
    geom_density(aes(y = ..density..), 
                 color = "Red") + 
    customtheme + 
    ggtitle("Degree Distribution") + 
    labs(x = "Degree", 
         y = "Density distribution") 
  
  
  # weighted vertex degree (strength)
  weightdeg <- strength(igraphlistw[[i]], mode = "all")
  p3 <- weightdeg %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    filter(X != 0) %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   binwidth = (max(weightdeg)-min(weightdeg))/length(weightdeg)*10, 
                   fill = "Blue") + 
    geom_density(aes(y = ..density..), 
                 color = "Red") + 
    customtheme + 
    ggtitle("Weighted Vertex Degree Distribution") + 
    labs(x = "Weighted Vertex Degree (Strength)", 
         y = "Density distribution")
  
  
  
  # edge betweenness
  edgebet <- edge_betweenness(igraphlist[[i]], 
                              directed = FALSE)
  p4 <- edgebet %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    filter(X != 0) %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   bins = (max(edgebet)-min(edgebet))/10, 
                   fill = "Blue") + 
    geom_density(aes(y = ..density..), 
                 color = "Red") + 
    customtheme + 
    ggtitle("Edge Betweenness (Unweighted)") + 
    labs(x = "Edge Betweenness", 
         y = "Density distribution")
  
  # nodal betweenness centrality 
  nbc <- (centr_betw(igraphlist[[i]], normalized = F)$res/
            centr_betw(igraphlist[[i]], normalized = F)$theoretical_max)
  p5 <- nbc %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    filter(X != 0) %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   fill = "Blue", 
                   binwidth = (max(nbc)-min(nbc))/length(nbc)*10) + 
    geom_density(aes(y = ..density..),
                 color = "red") +
    customtheme + 
    ggtitle("Nodal Betweenness Centrality (Unweighted)") + 
    labs(x = "Normalized Betweenness Centrality", 
         y = "Density distribution")
  
  
  # nodal eigenvector centrality 
  centreigen <- centr_eigen(igraphlistw[[i]], directed=F, normalized = T)$vector
  p6 <- centreigen %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    filter(X != 0) %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   binwidth = (max(centreigen)-min(centreigen))/length(centreigen)*10,  
                   fill = "Blue") + 
    geom_density(aes(y = ..density..), 
                 color = "Red") + 
    customtheme + 
    ggtitle("Nodal Eigenvector Centrality (Unweighted)") + 
    labs(x = "Normalized Eigenvector Centrality", 
         y = "Density distribution")
  
  # nodal closeness centrality 
  centrclos <- centr_clo(igraphlistw[[i]], normalized = T)$res %>% na.omit()
  p7 <- centrclos %>%
    as.data.frame() %>%
    rename(X = ".") %>%
    filter(X != 0) %>%
    ggplot(aes(x = X)) +
    geom_histogram(aes(y = ..density..), 
                   binwidth = (max(centrclos)-min(centrclos))/length(centrclos)*10,  
                   fill = "Blue") + 
    geom_density(aes(y = ..density..), 
                 color = "Red") + 
    customtheme + 
    ggtitle("Nodal Closeness Centrality (Unweighted)") + 
    labs(x = "Normalized Closeness Centrality", 
         y = "Density distribution")
  
  
  ## ARRANGE PLOTS 
  pall <- ggarrange(p1, p2, p3, p4, p5, p6, p7, 
                    labels = c("A", "B", "C", "D", "E", "F", "G"),
                    ncol = 3, 
                    nrow = 3)
  annotate_figure(pall, top = text_grob("SLR", 
                                        color = "red", face = "bold", size = 14))
  ggsave(paste("./data/graphs/", names(igraphlist)[i], "parameters.png", sep=""), 
         dpi = 300, 
         height = 10, 
         width = 10)
  
}




#### #####

