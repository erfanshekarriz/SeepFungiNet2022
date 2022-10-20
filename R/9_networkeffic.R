library(tidyverse)
library(phyloseq)
library(igraph)
library(ggridges)
library('brainGraph')


set.seed(123)
### UPLOAD FILES ###
files <- list.files(path = "./data/networks/filtered/nonweighted", 
                    full.names = TRUE)

files <- files[files !="./data/networks/filtered/nonweighted/Fu_igraphSLR_igraph.rds"]
networklist <- list()
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

for (file in files){
  networklist <- append(networklist, 
                        list(remove.0v(readRDS(file = file))))
}


names(networklist) <- c("B", "BA", "BAF", "BF")



#### EFFICIENCY ANALYSIS ####

efficlist <- list()
for (i in 1:length(networklist)){
  efficlist <- append(efficlist, 
                      list(as.data.frame(as.list(efficiency(networklist[[i]], 
                                            type = "nodal", 
                                            weights = NULL, use.parallel = TRUE)))))
}
names(efficlist) <- c("B", "BA", "BAF", "BF")

# combine data
effic.plot.tax <- bind_rows(efficlist, .id = "Network") %>% t() 
colnames(effic.plot.tax) <- effic.plot.tax[1,]
effic.plot.tax <- data.frame(effic.plot.tax[-1,])


effic.plot <- effic.plot.tax %>%
  remove_rownames() %>%
  pivot_longer(everything()) %>%
  mutate(value = as.numeric(value)) %>%
  drop_na() %>%
  rename(Network = name, 
         Efficiency = value) %>%
  group_by(Network) %>%
  mutate(meanEff = mean(Efficiency))

efficorder <- effic.plot %>%
  arrange(-meanEff) %>% 
  pull(Network) %>% 
  unique()

effic.plot <- effic.plot %>%
  mutate(Network = factor(Network, levels = efficorder))

effic.plottext <- effic.plot %>%
  group_by(Network) %>%
  summarize(mean = signif(mean(Efficiency), 2), 
            sd = signif(sd(Efficiency), 2))

#### PLOT ####
effic.plot %>% 
  ggplot(aes(x = as.numeric(Efficiency), 
             y= reorder(Network, -meanEff), 
             fill =Network, color = Network)) +
  geom_density_ridges(#fill ="white",
    scale = 1.4,
    bandwidth = 0.02,  
    jittered_points = TRUE, 
    point_size = 0.07, 
    point_alpha = 1,
    alpha = 0.7, 
    quantile_lines = TRUE, 
    quantiles = 2)+
  geom_text(data = effic.plottext, aes(y=Network, x=0.42, 
                                       label =paste("μ", "=", mean, "σ", "=", sd)), 
            fontface = "bold", size = 1.5, vjust=-1, parse = FALSE) + 
  scale_x_continuous(limits = c(0.05, 0.5)) + 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), angle = 90, vjust = 0.5, hjust=0)) + 
  theme_bw() + 
  scale_y_discrete(expand = expansion(add = c(0.5, 1.6))) + 
  # labs(title = "The Effect of Different Domains on Network Efficiency") +
  ylab("Network") + 
  xlab("Information Transfer Efficiency") + 
  theme() + 
  scale_color_hue(h = c(180, 350)) + 
  scale_fill_hue(h = c(180, 350)) + 
  theme(text = element_text(size = 10)) + 
  # scale_x_continuous(expand = c(0.008, 0.008)) +
  guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE)) + 
  theme(axis.text = element_text(size = 7), 
        axis.title.x = element_text(size = 8, vjust = -1), 
        axis.title.y = element_text(size = 8, vjust= 1), 
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"), 
        plot.background = element_blank(), 
        legend.background = element_blank(),
        panel.border =element_rect(colour="black", size = 0.7), 
        legend.position = "none")  + 
  scale_color_brewer(palette="Set1") + 
  scale_fill_brewer(palette="Set1")

ggsave("./data/graphs/Fig4B_9networkeffic.tiff",
       width = 8 ,
       height = 5,
       units = "cm",
       dpi = 1000 )




### T-TEST CHECKS
library(tidyverse)
library(ggpubr)
library(rstatix)

ggqqplot(effic.plot, x = "Efficiency", facet.by = "Network")

Efficiency <- effic.plot$Efficiency
Network <- effic.plot$Network
mod <- aov(Efficiency~Network)
summary(mod)

TukeyHSD(mod, conf.level = 0.99)

