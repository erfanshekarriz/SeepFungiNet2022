library(tidyverse)
library(phyloseq)
library(igraph)
library(brainGraph)
library(binr)
library(ggpubr)
library(RColorBrewer)

set.seed(1234)

# function to remove disconnected nodes
remove.0v <- function(ingraph) {
  outgraph <- ingraph
  cat("\nOriginal Graph Vertices: ", length(V(ingraph)))
  components <- igraph::clusters(ingraph, mode="weak")
  biggest_cluster_id <- which.max(components$csize)
  vert_ids <- V(ingraph)[components$membership != biggest_cluster_id]
  outgraph <- delete.vertices(outgraph, vert_ids) # remove them
  cat("\nFiltered Graph Vertices: ", length(V(outgraph)), "\n")
  return(outgraph)
}

igraphw <- remove.0v(readRDS("./data/networks/filtered/weighted/BAF_igraphSLR_igraphw.rds"))
igraph <- remove.0v(readRDS("./data/networks/filtered/nonweighted/BAF_igraphSLR_igraph.rds"))



physeq <- readRDS("./data/Data/physeqcombined.rds")
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
                                            "Normalized Degree",
                                            "Efficiency", 
                                            "Betweenness")), 
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
                        Pos = quantile(Value, prob=0.5)) 

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
                     size = 2.5, 
                     method = "wilcox.test", 
                     color = "orange", 
                     vjust = 1.6, bracket.size = 0.5) + 
  geom_text(data = dataMedian, aes(Domain, Pos, label = paste("m =", MD), color = Domain), 
            size = 1.5, vjust = -4, hjust= -0.15, fontface = "bold", angle=0) + 
  facet_wrap(~Factor, ncol = 1, scales = "free_y", strip.position="top") + 
  theme_bw() +
  theme(strip.background = element_rect(fill="grey99", color="grey", size =1), 
        strip.text = element_text(face = "bold.italic", size = 8), 
        axis.title = element_blank(), 
        axis.text =  element_text(size = 7), 
        axis.text.x =  element_text(color = brewer.pal(n=3,"Set1"), face ="bold"),
        legend.text =   element_text(size = 6), 
        legend.background = element_blank(),
        legend.title =   element_text(size = 6),
        panel.border =element_rect(colour="grey", size = 1), 
        panel.background = element_rect(fill = "grey90"), 
        plot.background = element_blank(), 
        panel.grid = element_line(color = "grey95"), 
        legend.position = "none") + 
  expand_limits(y = 1) + 
  scale_color_brewer(palette = "Set1")
  


ggsave("./data/graphs/Fig1D_10domainspecBAF.tiff",
       width = 4.5,
       height = 13,
       units = "cm",
       dpi = 1000 )




