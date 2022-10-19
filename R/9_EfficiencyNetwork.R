library(tidyverse)
library(phyloseq)
library(igraph)
# library(fantaxtic)
library(ggridges)
library('brainGraph')


### UPLOAD FILES ###
files <- list.files(path = "./data/networks/attackanalysis", 
                    full.names = TRUE)
networklist <- list()
remove.0v <- function(ingraph) {
  library(igraph)
  outgraph <- ingraph
  cat("\nOriginal Graph Vertices: ", length(V(ingraph)))
  V(outgraph)$names <- paste0("v", 1:length(V(outgraph))) #assign unique names to vertices to keep track
  disconnected.V <- which(degree(outgraph)==0) # index completely disconnected vertices
  outgraph <- delete.vertices(outgraph, disconnected.V) # remove them
  cat("\nFiltered Graph Vertices: ", length(V(outgraph)), "\n")
  return(outgraph)
}

for (file in files){
  networklist <- append(networklist, 
                        list(remove.0v(readRDS(file = file))))
}


names(networklist) <- c("B", "BA", "BAF", "BF", "Fu")



#### EFFICIENCY ANALYSIS ####

efficlist <- list()
for (i in 1:length(networklist)){
  efficlist <- append(efficlist, 
                      list(as.data.frame(as.list(efficiency(networklist[[i]], 
                                            type = "nodal", 
                                            weights = NULL, use.parallel = TRUE)))))
}
names(efficlist) <- c("B", "BA", "BAF", "BF", "Fu")

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


#### PLOT ####
effic.plot %>% 
  filter(Efficiency > 0.2) %>% 
  ggplot(aes(x = as.numeric(Efficiency), 
             y= reorder(Network, -meanEff), 
             fill =Network, color = Network)) +
  geom_density_ridges(#fill ="white",
    scale = 1.4,
    bandwidth = 0.04,  
    jittered_points = TRUE, 
    point_size = 0.07, 
    point_alpha = 1,
    alpha = 0.7, 
    quantile_lines = TRUE, 
    quantiles = 2)+
  scale_x_continuous(limits = c(0.0, 0.7)) + 
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
  guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE))  


ggsave("./data/graphs/efficiency.png",
       width = 10,
       height = 7,
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

