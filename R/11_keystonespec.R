library(tidyverse)
library(phyloseq)
library(igraph)
library(brainGraph)
library(binr)
library(ggrepel)

set.seed(1234)


#### LOAD DATA ####
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


#### CALCULATE DEGREE & BETWEENNESS #### 
plotdf <- data.frame(Degree =degree(igraph, normalized = TRUE), 
           Betweenness = betweenness(igraph, normalized = TRUE)) %>%
  cbind(taxa) %>% 
  mutate(Domain = if_else(Domain=="Eukaryota", "Fungi", Domain))

quant <- 0.8
plotdfcol <- plotdf %>%
  group_by() %>%
  mutate(quant90deg = quantile(Degree, quant), 
         quant90bet = quantile(Betweenness, quant)) %>%
  filter(Degree >= quant90deg, 
         Betweenness >= quant90bet) %>%
  add_column(label = NA) %>%
  mutate(label = Species)

plotdfncol <- plotdf %>%
  group_by() %>%
  mutate(quant90deg = quantile(Degree, quant), 
         quant90bet = quantile(Betweenness, quant)) %>%
  filter(Degree < quant90deg, 
         Betweenness < quant90bet)

degquant <- unique(plotdfcol$quant90deg)
betquant <- unique(plotdfcol$quant90bet)


#### MANUALLY ADD LABELS ####
plotdfcol[1,"label"] <- "Acinetobacter sp."
plotdfcol[2,"label"] <- "SEEP-SRB1 sp."
plotdfcol[3,"label"] <- "SEEP-SRB1 sp."
plotdfcol[4,"label"] <- "Sulfurovum sp."
plotdfcol[5,"label"] <- "Aminicenantales sp."
plotdfcol[6,"label"] <- "Thiomicrospiraceae"
plotdfcol[7,"label"] <- "Corynebacterium tuberculostearicum"
plotdfcol[8,"label"] <- "Desulfatiglans sp."
plotdfcol[9,"label"] <- "Hydrogenophilus sp."
plotdfcol[10,"label"] <- "JS1 sp."
plotdfcol[11,"label"] <- "JS1 sp."
plotdfcol[12,"label"] <- "JS1 sp."
plotdfcol[13,"label"] <- "Gemmataceae"
plotdfcol[14,"label"] <- "Staphylococcus sp."
plotdfcol[15,"label"] <- "ANME-1b sp."
plotdfcol[16,"label"] <- "Marine Benthic Group D and DHVEG-1"
plotdfcol[17,"label"] <- "Metschnikowia sp."
plotdfcol[18,"label"] <- "Cladosporiaceae"
plotdfcol[19,"label"] <- "Arthrinium sp."
plotdfcol[20,"label"] <- "Fusarium oxysporum"
plotdfcol[21,"label"] <- "Derxomyces"
plotdfcol[22,"label"] <- "Suillus lakei"
plotdfcol[23,"label"] <- "Fusarium solani" # verified by BLAST search
plotdfcol[24,"label"] <- "Malassezia sp."



#### PLOT ####

plotdfncol %>%
  ggplot(aes(x = Degree, y = Betweenness, shape = Domain)) +
  geom_hline(yintercept = betquant-0.005, color = "orange", linetype="dashed", size = 0.4) + 
  geom_vline(xintercept = degquant-0.005, color = "orange", linetype="dashed", size = 0.4) +
  geom_jitter(color = "grey", size = 1) +
  geom_jitter(data = plotdfcol, aes(x=Degree, y=Betweenness, color = Domain), size=1) + 
  geom_label_repel(data = plotdfcol, 
                   max.overlaps=20,
                   force = 10,
                  min.segment.length = Inf,
             aes(x=Degree, y=Betweenness, label = label, color = Domain), 
             size = 2, 
             fill = NA,
             fontface = "bold") + 
  scale_color_brewer(palette = "Set1") + 
  theme_bw() + 
  xlab("Normalized Degree")+ 
  ylab("Normalized Betweenness")  + 
  theme(axis.text =  element_text(size = 10), 
        axis.title.x = element_text(size = 12, vjust = -2), 
        axis.title.y = element_text(size = 12, vjust= 4), 
        legend.text =   element_text(size = 9), 
        legend.title =   element_text(size = 9),
        legend.background = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.border =element_rect(colour="black", size = 1.3), 
        panel.grid = element_line(color = "grey95"),        
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(x = 0.08, y = 0.165, label = paste0("TOP ", 100-(quant*100), "th PERCENTILE"), size = 3, 
            fontface = "bold.italic", color = "orange")
  

  


ggsave("./data/graphs/Fig4C_11keystonespec.tiff",
       width = 15,
       height = 12,
       units = "cm",
       dpi = 1000 )

