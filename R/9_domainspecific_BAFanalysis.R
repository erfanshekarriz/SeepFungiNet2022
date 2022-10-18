library(tidyverse)
library(phyloseq)
library(igraph)
library(brainGraph)
library(binr)

set.seed(1234)

# function to remove disconnected nodes
remove.0v <- function(ingraph) {
  library(igraph)
  outgraph <- ingraph
  cat("\nOriginal Graph Vertices: ", length(V(ingraph)))
  V(outgraph)$names <- paste0("v", 1:length(V(outgraph))) #assign unique names to vertices to keep track
  disconnected.V <- which(degree(outgraph)<2) # index completely disconnected vertices
  outgraph <- delete.vertices(outgraph, disconnected.V) # remove them
  cat("\nFiltered Graph Vertices: ", length(V(outgraph)), "\n")
  return(outgraph)
}

igraphw <- remove.0v(readRDS("./data/networks/BAF/BAF_igraphSLR_igraphw.rds"))
igraph <- remove.0v(readRDS("./data/networks/BAF/BAF_igraphSLR_igraph.rds"))



physeq <- readRDS("./data/rawanalysis/physeqcombined.rds")
taxa <- data.frame(tax_table(physeq))[as_ids(V(igraph)),]


#### DEGREE ####
degrees <- data.frame(Value = degree(igraphw, normalized = FALSE)) %>%
  cbind(., taxa) %>%
  mutate(Domain = if_else(Domain=="Eukaryota", "Fungi", Domain)) %>%
  select(Value, Domain) %>%
  add_column(Factor = "Normalized Degree") %>%
  filter(Value>0.00)


#### EFFICIENCY #### 
efficiency <- data.frame(Value = efficiency(igraph)) %>%
  cbind(., taxa) %>%
  mutate(Domain = if_else(Domain=="Eukaryota", "Fungi", Domain)) %>%
  select(Value, Domain) %>%
  add_column(Factor = "Efficiency") %>%
  filter(Value>0.00)


#### BETWEENNESS #### 
betweenness <- data.frame(Value = betweenness(igraph, normalized = FALSE)) %>%
  cbind(., taxa) %>%
  mutate(Domain = if_else(Domain=="Eukaryota", "Fungi", Domain)) %>%
  select(Value, Domain) %>%
  add_column(Factor = "Betweenness") %>%
  filter(Value>0.00)


#### WEIGHTED AVERAGE ####
edgeweight <- as_edgelist(igraphw) %>%
  as.data.frame() %>%
  rename(from = V1, to = V2) %>%
  mutate(fromtaxa = taxa[from,"Domain"], 
         totaxa = taxa[to,"Domain"], 
         weights = E(igraphw)$weight) %>%
  add_column(Fungi = NA, 
             Bacteria = NA, 
             Archaea = NA) %>%
  mutate(Fungi = if_else(fromtaxa == "Eukaryota" | totaxa == "Fungi", "Yes", "No"), 
         Bacteria = if_else(fromtaxa == "Bacteria" | totaxa == "Bacteria", "Yes", "No"), 
         Archaea = if_else(fromtaxa == "Archaea" | totaxa == "Archaea", "Yes", "No")) %>%
  pivot_longer(cols = c(Fungi, Bacteria, Archaea), 
               names_to = "Domain", 
               values_to = "Weight") %>%
  mutate(Value = if_else(Weight == "Yes", abs(weights), 0)) %>%
  filter(Value != 0) %>%
  select(Domain, Value) %>%
  add_column(Factor = "Absolute Edge Weight") %>%
  filter(Value!=0.0)





#### MERGE & PLOT ####
# matchrow <- intersect(rownames(degrees), 
#                       rownames(efficiency))
# degeff <- inner(degrees[matchrow, ], efficiency[matchrow,])
plotdf <- rbind(betweenness, degrees, efficiency, edgeweight) %>%
  mutate(Factor = factor(Factor, levels = c("Absolute Edge Weight",
                                            "Efficiency", 
                                            "Betweenness", 
                                            "Normalized Degree")), 
         Domain = factor(Domain, levels = c("Fungi", "Bacteria", "Archaea")))  %>%
  filter(Factor!= "Efficiency")
  


rankcut <- 10
selectcolor <- plotdf %>%
  group_by(Factor) %>%
  mutate(rank = rank(Value)) %>%
  filter(rank > (max(rank)-rankcut))  

horizline <- selectcolor %>%
  group_by(Factor) %>%
  summarize(min = min(Value)) %>%
  add_column("text" = NA) %>%
  mutate(text = as.character(min)) 
dataMedian <- summarise(group_by(plotdf, Factor, Domain), 
                        MD = signif(median(Value), 2), 
                        Pos = quantile(Value, prob=0.80)) 

selectnocolor  <-  plotdf %>%
  group_by(Factor) %>%
  mutate(rank = rank(Value)) %>%
  filter(rank < (max(rank)-rankcut))  

my_comparisons <- list(c("Bacteria", "Archaea"),
                       c("Bacteria", "Fungi"),
                       c("Fungi", "Archaea"))

selectnocolor %>%
  ggplot(aes(x = Domain, y = Value, color = Domain)) +
  geom_jitter(size = 0.01, width=0.1, color = "grey", alpha = 0.4) + 
  geom_jitter(data=selectcolor, aes(x = Domain, y = Value), 
              color = "Red", size = 0.5, width=0.1, shape = 4) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # geom_hline(data = horizline, aes(yintercept = min),
  #            color = "red", size =0.2, linetype = "dashed") +
  stat_compare_means(aes(label = ..p.signif..), 
                     comparisons = my_comparisons,
                     size = 3, 
                     method = "wilcox.test", 
                     color = "gray70", 
                     vjust = 1.6, bracket.size = 0.5) + 
  geom_text(data = dataMedian, aes(Domain, Pos, label = paste("m =", MD), color = Domain), 
            size = 2, vjust = -2, fontface = "bold") + 
  facet_wrap(~Factor, ncol = 1, scales = "free_y", strip.position="top") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="grey99", color="grey", size =1), 
        strip.text = element_text(face = "italic", size = 6), 
        axis.title = element_blank(), 
        axis.text =  element_text(size = 5), 
        legend.text =   element_text(size = 6), 
        legend.title =   element_text(size = 6),
        panel.border =element_rect(colour="grey", size = 1)) + 
  expand_limits(y = 1) 


ggsave("./data/graphs/VAFanalysis.png",
       width = 7,
       height = 13,
       units = "cm",
       dpi = 1000 )




