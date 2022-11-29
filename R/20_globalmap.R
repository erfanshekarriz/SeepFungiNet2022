library(tidyverse)
library(RColorBrewer)
library(ggrepel)


world <- map_data("world")


df <- read.csv("./data/Data/globaldataset.csv") %>%
  distinct(Site, .keep_all = TRUE)

nb.cols <- nrow(df)
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

seedp <- 420

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "grey95", fill = "grey93", size = 0.1) + 
  geom_point(data = df,
    aes(degree_E, degree_N, color = Site),
    alpha = 0.7) +
  geom_label_repel(data = df,
                   aes(degree_E, degree_N, fill = Site, label=Site), 
                   color="black", size=3, fontface="bold.italic", alpha=0.7, 
                   seed = seedp) + 
  geom_label_repel(data = df,
                   aes(degree_E, degree_N,label=Site), 
                   color="black", size=3, fontface="bold.italic", fill=NA, 
                   seed = seedp) + 
  scale_color_manual(values = mycolors) + 
  scale_fill_manual(values = mycolors) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(panel.background = element_rect(fill="grey99"), 
        panel.border = element_rect(color="black", fill=NA, size=1.2, 
                                    linetype = "solid"), 
        axis.title = element_text(face="italic"), 
        axis.text = element_text(face="italic"), 
        legend.position =  c(0.5,0.24), 
        legend.direction="horizontal",
        legend.background = element_blank()) + 
  scale_discrete_identity(
    aesthetics = "label",
    name = "",
    breaks = df$Site,
    labels = df$Site.Legend,
    guide = "legend") + 
  guides(color = FALSE, fill=FALSE, 
         label=guide_legend(nrow=4, byrow=TRUE))
  

ggsave("./data/graphs/Fig6A_23globalmap.tiff", 
       units="cm", 
       height=10, 
       width = 17, 
       dpi=1000)
