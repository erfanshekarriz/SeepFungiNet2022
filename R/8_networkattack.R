# libraries
library(ggplot2)
library(igraph)
library(ggsignif)
library(tidyverse)
library(rlist)
library(ggpubr)



#### NOTE: This piece of code only works on the LARGEST CONNECTED COMPONENT
#          of your graphs. It will automatically remove all other vertices
#          and report the number of nodes removed as it processes. 


set.seed(324)
npermutations <- 100 # this number has been tested to be enough to capture variation
                     # in data while maintaing a clear plot
alphapoints <- 0.01 # sets the opacity for plots 


#FUNCTIONS####################################


# function to remove vertices with no edges
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
# function to cbind even if lenght of two columns aren't matching
cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
attack.robustness.random <- function(g, perm = 50){
  n <- length(V(g)) # get the number of vertices as an integer
  mat <- matrix(ncol=2,nrow=n, 0) # make an empty matrix with vertices as rows
  mat[,1] <- as_ids(V(g)) #add vertices names to matrix
  random.attack <- list(plotting_values = NULL,
                        permutation=NULL,
                        R= NULL,
                        V = NULL,
                        variance=NULL,
                        mean_xy =list(x=NULL, y=NULL),
                        raw_value_table=list(NULL),
                        method = "randomperm")
  for (j in 1:perm) {
    # j <- 1
    matri <- as.matrix(mat[sample(nrow(mat)),]) #randomly shuffle
    g2 <- g # make a copy of the graph
    clustersizes<-integer(n-1)
    for(i in 1:(n-1)){
      g2 <- delete.vertices(g2, v=which(as_ids(V(g2))==matri[i,1])) #delete by name
      maxcsize2 <- max(clusters(g2)$csize)
      clustersizes[i]<-maxcsize2
    }
    df<-as.data.frame(cbind(matri, c(clustersizes, NA))) # returns a data-frame
    df[,4] <- 1:nrow(df) / nrow(df) # Fraction of nodes removed
    df[,5] <- as.numeric(as.character(df[,3])) / n # Fraction of Largest Component 
    df[1,6] <- (sum(df[,5], na.rm = TRUE))/n #Robustness
    df[1,7] <- 0.5-(df[1,6]) # Vulnerability 
    options(warn=-1)
    toprow <- cbind(0, 0, (max(clusters(g)$csize))/n , 0, 1)
    df[,1:5] <- rbind(toprow, df[,1:5])
    options(warn=-0)
    
    names(df)<-c("Vertex Index Removed",
                 "Betweenness",
                 "Size of Largest Component ",
                 "Fraction of Nodes Removed",
                 "Fraction Size of Largest Component",
                 "R",
                 "V")
    df[is.na(df)] = 0
    
    random.attack$raw_permutation_values[[j]] <- df
  }
  
  random.attack$plotting_values <- list.rbind(lapply(random.attack$raw_permutation_values, `[`, 4:5))
  random.attack$variance <- var(random.attack$plotting_values[,2])
  
  xlist <- lapply(random.attack$raw_permutation_values, `[`, 4)
  xvalues <- matrix(ncol = 1, nrow = n)
  ylist <- lapply(random.attack$raw_permutation_values, `[`, 5)
  yvalues <- matrix(ncol = 1, nrow = n)
  
  for(w in 1:n){
    xvalues[w,] <- mean(unlist(lapply(xlist, `[`, w, 1)))
    yvalues[w,] <- mean(unlist(lapply(ylist, `[`, w, 1)))
  }
  random.attack$mean_xy$x <- xvalues
  random.attack$mean_xy$y <- yvalues
  
  
  random.attack$R <- mean(unlist(lapply(random.attack$raw_permutation_values, `[`, 1, 6)))
  random.attack$V <- mean(unlist(lapply(random.attack$raw_permutation_values, `[`, 1, 7)))
  random.attack$V.list <-  unlist(lapply(random.attack$raw_permutation_values, `[`, 1, 7))
  random.attack$permutation <- perm
  
  cat("\n#Permutations: ",  perm,
      "\nMean Normalized Robustness (R): ", random.attack$R,
      "\nMean Normalized Vulnerability (V): ",  random.attack$V,
      "\nRobustness Variance: ", random.attack$variance)
  
  return(random.attack)
}
attack.robustness.betweenness <- function(g){
  n <- length(V(g)) # get the number of vertices as an integer
  mat <- matrix(ncol=2,nrow=n, 0) # make an empty matrix with vertices as rows
  mat[,1] <- as_ids(V(g)) #add vertices names to matrix
  bet <- betweenness(g)  #calculate the betweenness of the vertex
  mat[,2] <- bet # add the degree to the graph
  matri <- mat[order(as.numeric(mat[,2]), decreasing =TRUE ),] #order in terms of decreasing between-ness
  betweenness.attack <- list(plotting_values = NULL,
                             permutation=NULL,
                             R= NULL,
                             V = NULL,
                             variance=NULL,
                             mean_xy =list(x=NULL, y=NULL),
                             raw_value_table=list(NULL),
                             method = "betweenness")
  g2 <- g # make a copy of the graph
  clustersizes<-integer(n-1)
  for(i in 1:(n-1)){
    g2 <- delete.vertices(g2, v=which(V(g2)$names==matri[i,1])) #index by name
    maxcsize2 <- max(clusters(g2)$csize)
    clustersizes[i]<-maxcsize2
  }
  df<-as.data.frame(cbind(mat, c(clustersizes, NA))) # returns a data-frame
  df[,4] <- 1:nrow(df) / nrow(df)
  df[,5] <- as.numeric(as.character(df[,3])) / n
  df[1,6] <- (sum(df[,5], na.rm = TRUE))/n
  df[1,7] <- 0.5-(df[1,6])
  options(warn=-1)
  toprow <- cbind(0, 0, (max(clusters(g)$csize))/n , 0, 1)
  df[,1:5] <- rbind(toprow, df[,1:5])
  options(warn=-0)
  
  names(df)<-c("Vertex Index Removed",
               "Betweenness",
               "Size of Largest Component ",
               "Fraction of Nodes Removed",
               "Fraction Size of Largest Component",
               "R",
               "V")
  df[is.na(df)] = 0
  
  
  betweenness.attack$plotting_values <- df[,4:5]
  betweenness.attack$permutation <- "only available with the (random.perm) method"
  betweenness.attack$R <- df[1,6]
  betweenness.attack$V <- df[1,7]
  betweenness.attack$variance <- "only available with the (random.perm) method"
  betweenness.attack$mean_xy <- "only available with the (random.perm) method"
  betweenness.attack$raw_value_table <- df
  
  cat("Robustness (R): ", betweenness.attack$R,
      "\nVulnerability (V): ",  betweenness.attack$V)
  
  return(betweenness.attack)
}
attack.robustness.degree <- function(g){
  n <- length(V(g)) # get the number of vertices as an integer
  mat <- matrix(ncol=2,nrow=n, 0) # make an empty matrix with vertices as rows
  mat[,1] <- as_ids(V(g)) #add vertices names to matrix
  deg <- degree(g) #calculate the degree of the vertex
  mat[,2] <- deg # add the degree to the graph
  matri <- mat[order(as.numeric(mat[,2]), decreasing =TRUE ),] #order in terms of decreasing between-ness
  g2 <- g # make a copy of the graph
  degree.attack <- list(plotting_values = NULL,
                        permutation=NULL,
                        R= NULL,
                        V = NULL,
                        variance=NULL,
                        mean_xy =list(x=NULL, y=NULL),
                        raw_value_table=list(NULL),
                        method = "degree")
  clustersizes<-integer(n-1)
  for(i in 1:(n-1)){
    g2 <- delete.vertices(g2, v=which(V(g2)$names==matri[i,1])) #index by name
    maxcsize2 <- max(components(g2)$csize)
    clustersizes[i]<-maxcsize2
  }
  df<-as.data.frame(cbind(mat, c(clustersizes, NA))) # returns a data-frame
  df[,4] <- 1:nrow(df) / nrow(df)
  df[,5] <- as.numeric(as.character(df[,3])) / n
  df[1,6] <- (sum(df[,5], na.rm = TRUE))/n
  df[1,7] <- 0.5-(df[1,6])
  options(warn=-1)
  toprow <- cbind(0, 0, (max(clusters(g)$csize))/n , 0, 1)
  df[,1:5] <- rbind(toprow, df[,1:5])
  options(warn=-0)
  
  names(df)<-c("Vertex Index Removed",
               "Betweenness",
               "Size of Largest Component ",
               "Fraction of Nodes Removed",
               "Fraction Size of Largest Component",
               "R",
               "V")
  df[is.na(df)] = 0
  degree.attack$plotting_values <- df[,4:5]
  degree.attack$permutation <- "only available with the (random.perm) method"
  degree.attack$R <- df[1,6]
  degree.attack$V <- df[1,7]
  degree.attack$variance <- "only available with the (random.perm) method"
  degree.attack$mean_xy <- "only available with the (random.perm) method"
  degree.attack$raw_value_table <- df
  
  cat("Robustness (R): ", degree.attack$R,
      "\nVulnerability (V): ",  degree.attack$V)
  
  return(degree.attack)
}
attack.robustness <- function(graphlist, method = "randomperm", permutation_n = 50) {
  if(missing(graphlist)){stop("Please input the correct concatenated LIST of your igraph networks to analyse...")}
  if(missing(method)){stop("Method missing...\nTry ?attack_robustness to view available methods...")}
  permutationz <- permutation_n
  attack.network.list <- list()
  for (i in 1:length(graphlist)) {
    cat("\n\nNow Processing: ", names(graphlist)[i])
    filtered.graph <- remove.0v(graphlist[[i]])
    if (method == "randomperm"){
      #if(missing(permutationz) | is.integer(permutationz)==FALSE){warning("\nPermutation parameter either left blank or not an integer: \n 50 Permutations was used by default if (randomperm) method is chosen")}
      cat("\nMethod: RANDOMPERM")
      attack.network.list[[i]] <- attack.robustness.random(filtered.graph, perm = permutationz)
    } else if (method == "betweenness"){
      attack.network.list[[i]] <- attack.robustness.betweenness(filtered.graph)
      cat("\nMethod: BETWEENNESS")
    } else if (method == "degree"){
      attack.network.list[[i]] <- attack.robustness.degree(filtered.graph)
      cat("\nMethod: DEGREE")
    } else {
      stop("Invalid method. Try ?attack_robustness to view available methods.")
    }
  }
  names(attack.network.list) <- names(graphlist)
  return(attack.network.list)
  print("done")
}
plot.robustness <- function(x, vulnerability=TRUE){
  require(ggplot2)
  vulnerability <- vulnerability
  
  # 1) initialize lists and matrices
  plot.matrix <- matrix()
  barplot.matrix <- matrix()
  plot.matrix.list <- list()
  barplot.matrix.list <- list()
  
  # 2) loop through the graphs and create a LIST containing extracted plotting values
  for (i in 1:length(x)) {
    plot.matrix.list[[i]] <- cbind(rbind(x[[i]]$plotting_values), Network = names(x)[i])
    barplot.matrix.list[[i]] <- cbind(rbind(x[[i]]$R), rbind(x[[i]]$V) , Network = names(x)[i])
  }
  
  # 3) convert list into a clean plottable data frame
  library(rlist)
  plot.matrix <- list.rbind(lapply(plot.matrix.list, `[`, 1:3))
  plot.df <- as.data.frame(plot.matrix)
  barplot.matrix <- list.rbind(lapply(barplot.matrix.list, `[`, 1:3))
  barplot.df <- as.data.frame(barplot.matrix)
  colnames(barplot.df) <- c("R", "V", "Network")
  barplot.df <- barplot.df
  
  if (x[[1]]$method == "randomperm"){
    random.matrix.plot.list <- list()
    for (j in 1:length(x)) {
      random.matrix.plot.list[[j]] <- cbind(as.numeric(attack.robustness.ran[[1]]$mean_xy$x), 
                                            as.numeric(attack.robustness.ran[[1]]$mean_xy$y), 
                                            Network = names(x)[j])
    }
    random.matrix.plot.list <- random.matrix.plot.list
    random.matrix.plot.matrix <- list.rbind(random.matrix.plot.list)
    random.matrix.plot.matrix <- random.matrix.plot.matrix
    random.matrix.plot.df <- as.data.frame(random.matrix.plot.matrix)
    colnames(random.matrix.plot.df) <- c("Fraction of Nodes Removed",
                                         "Fraction Size of Largest Component",
                                         "Network")
    random.matrix.plot.df <- random.matrix.plot.df
    random.matrix.plot.df <- random.matrix.plot.df %>%
      mutate(Network =  factor(Network, levels = names(graph.list))) %>%
      arrange(Network) 
    
    legendord <- names(sort(unlist(lapply(x, function(x) x$V))))
    attack.robustness.ran
    
    plot.df  <- plot.df %>%
      mutate(Network = factor(Network, levels = legendord)) %>%
      group_by(Network) %>%
      mutate(Fractionbin = cut(`Fraction of Nodes Removed`, breaks = 100), 
             Fractionbinnum = `Fraction of Nodes Removed`[as.numeric(Fractionbin)]) %>%
      group_by(Fractionbin) %>%
      mutate(meanFraction = mean(`Fraction Size of Largest Component`))
    
    plot.main <- ggplot(plot.df, 
                        aes(x= as.numeric(`Fraction of Nodes Removed`), 
                            y = as.numeric(`Fraction Size of Largest Component`), 
                            color = Network))+
      geom_point(size = 0.2, alpha=alphapoints) +
      geom_line(aes(x = `Fraction of Nodes Removed`, y = meanFraction, group=Network)) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme(
        # panel.grid.major = element_blank(),
        #     panel.grid.minor = element_blank(),
        #     panel.background = element_blank(),
        axis.title.x=element_text(margin = margin(t =10), size = 12),
        axis.title.y =element_text(margin = margin(r =10), size=12), 
        axis.text=element_text(size=10)) +
      labs(x = "Fraction of Nodes Removed",
           y = "Fraction Size of Largest Component",
           color = "Network\n") + 
      theme_bw() + 
      scale_x_continuous(breaks=seq(0, 1.0, by = 0.05), expand = c(0, 0)) + 
      scale_y_continuous(breaks=seq(0, 1.0, by = 0.05), expand = c(0, 0)) + 
      scale_color_brewer(palette = "Set1")
    
  }
  else {
    plot.main <- ggplot(plot.df, aes(x= `Fraction of Nodes Removed`, y = `Fraction Size of Largest Component`, color = Network))+
      geom_line(size = 0.8) +
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      theme(
        # panel.grid.major = element_blank(),
        #     panel.grid.minor = element_blank(),
        #     panel.background = element_blank(),
        axis.title.x=element_text(margin = margin(t =10), size = 12),
        axis.title.y =element_text(margin = margin(r =10), size=12),
        axis.text=element_text(size=10)) +
      labs(x = "Fraction of Nodes Removed",
           y = "Fraction Size of Largest Component",
           color = "Network\n") +
      expand_limits(x = 0, y = 1) + 
      theme_bw() + 
      scale_x_continuous(breaks=seq(0, 1.0, by = 0.05), expand = c(0, 0)) + 
      scale_y_continuous(breaks=seq(0, 1.0, by = 0.05), expand = c(0, 0)) + 
      scale_color_brewer(palette = "Set1")
    
  }
  if (vulnerability){
    barplot.inset <- ggplot(barplot.df, aes(x= Network, y = as.numeric(V), fill = Network)) +
      geom_bar(stat = "identity", width = 0.7, alpha=0.8, color = "black") +
      geom_text(aes(label= round(as.numeric(V), digits = 2)), vjust=-0.5, size = 3) + 
      scale_y_continuous(breaks=c(0, 0.25, 0.5)) +
      coord_cartesian(ylim=c(0,0.5)) +
      xlab("Networks") +
      ylab("V") +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, size = 7),
            axis.title.x=element_blank(),
            axis.title.y =element_text(size=8.5),
            axis.text.y =element_text(size=6)) +
      theme(legend.position= "none")  #+
    #geom_signif(comparisons = list(c("Prokaryotes.S", "CrossD.S")),
    #map_signif_level = TRUE, textsize = 3)
    
    plot.final <- plot.main + annotation_custom(ggplotGrob(barplot.inset), xmin = 0.72, ymax = 0.99,
                                                ymin = 0.70, xmax = 0.99)
  }
  else if (!vulnerability){
    plot.final <- plot.main
  }
  return(plot.final)
}


#LOAD_DATA####################################
files <- list.files("./data/networks/filtered/nonweighted/", 
                    recursive=FALSE, 
                    full.names = TRUE)

filesshort <- list.files("./data/networks/filtered/nonweighted/", 
                         recursive=FALSE, 
                         full.names = FALSE)

graph.list <- list()
for (i in 1:length(files)){
  graph.list <- append(graph.list, list(readRDS(files[i])))
  names(graph.list)[i] <- filesshort[i]
}


# add custom names to the graphs
names(graph.list) <- c("B", "BA", 
                       "BAF", "BF")


#RUN####################################
# randomperm 
# betweennness
# degree
attack.robustness.ran <- attack.robustness(graph.list, 
                                           method="randomperm", 
                                           permutation_n = npermutations)


plot.robustness(attack.robustness.ran, vulnerability = F) + 
  theme(axis.text = element_text(size = 8), 
        axis.title.x = element_text(size = 9, vjust = -4), 
        axis.title.y = element_text(size = 9, vjust= 4), 
        plot.margin = margin(0.5,0,0.8,0.8, "cm"))

ggsave("./data/graphs/Fig4A1_8networkattack.png" ,
       width = 18,
       height = 15,
       units = "cm",
       dpi = 1000 )



#PVALUES####################################

# choose which networks you'd like to compare
my_comparisons <- list(c("B", "BF"),
                       c("BF", "BAF"), 
                       c("BA", "BF"))

tofind <- paste(c("BAF","BA", "BF", "B"), collapse="|")

legendord <- names(sort(unlist(lapply(attack.robustness.ran, function(x) x$V))))

stack(lapply(attack.robustness.ran, function(l) l[['V.list']])) %>%
  rename(Network = ind, Vulnerability = values) %>%
  mutate(Network = factor(Network, levels=legendord)) %>%
  group_by(Network) %>%
  mutate(meanV = mean(Vulnerability), 
         method = unlist(str_extract_all(Network, tofind))) %>%
  ggplot(aes(y = Vulnerability, 
             x = reorder(Network, desc(meanV)), 
             color = factor(method))) +
  geom_jitter(size = 0.001, alpha = 1, height=0.1, width=0.2, color = "grey98") + 
  geom_boxplot(fill = "NA", outlier.shape = NA, width=0.2) +
  stat_compare_means(aes(label = ..p.signif..), 
                     comparisons = my_comparisons,
                     size = 4, 
                     method = "wilcox.test", 
                     color = "gray20", 
                     vjust=0.5,
                     label.y = c(0.43, 0.3, 0.38), 
  ) + 
  stat_compare_means(method= "anova", 
                     label.y = 0.49, 
                     label.x = length(attack.robustness.ran)-1.35,  
                     size = 2.2, 
                     color = "black", 
                     fontface = "bold") + 
  expand_limits(y = 0.50) +  # adjust y axis expansion
  theme_bw() +
  theme(
    # panel.grid.major = element_blank(),
    axis.title.x=element_text(margin = margin(t =10), size = 8),
    axis.title.y =element_text(margin = margin(r =10), size=8), 
    axis.text=element_text(size=8), 
    plot.background = element_blank(),
    panel.background = element_rect(fill = "grey90", color = "black"),
    panel.grid = element_line(color="white"),
    legend.position = "none") +
  labs(x = "Network",
       y = "Vulnerability (V)",
       color = "Network")  + 
  scale_color_manual(values=c("#4DAF4A" , "#984EA3", "#E41A1C", "#377EB8"))

ggsave("./data/graphs/Fig4A2_8networkattack.png" ,
       width = 5.1,
       height = 5,
       units = "cm",
       dpi = 1000 )


######
