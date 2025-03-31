# Script morphospace 

library("morphr")
library("chroma")
library("castr")
library("patchwork")
library("lubridate")
library("tidyverse")
library(dplyr)
library(imager)
library("FactoMineR")
library("factoextra")


# compute euclidean distance between each individual coordinates and cluster center coordinates 
library(vegan)

# Function to pre-process UVP images
preprocess <- function(x) {
  x |>
    # remove 31 pixels from the bottom (=the scale bar)
    img_chop(bottom=31) |>
    # change the gamma to see the light objects better
    img_adjust_gamma(gamma=0.7)
}

#write.csv(fordistCENTER1,'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/fordistCENTER1.csv', row.names = FALSE)
#write.csv(fordistCENTER2,'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/fordistCENTER2.csv', row.names = FALSE)
#write.csv(fordistCENTER3,'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/fordistCENTER3.csv', row.names = FALSE)
#write.csv(fordistCENTER4,'C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/fordistCENTER4.csv', row.names = FALSE)

#fordistCENTER3 <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/fordistCENTER3.csv", sep = ",", header = TRUE)
#distances3F <- dist(fordistCENTER3, method = "euclidean")[2:(nrow(indiv_clust4) + 1)]

qQR=c(0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95)
#qQR=c(0.75)

for (QR in qQR) 
{
  # cluster 1 
  
  q_1 <- quantile(distances1, QR)
  indiv_clust01 <- indiv_clust1 %>% mutate(distance_to_center = distances1,
                                          mean_distance_to_center = mean(distances1), 
                                          qrtile = q_1)
  # cluster 2
  
  q_2 <- quantile(distances2, QR)
  indiv_clust02 <- indiv_clust2 %>% mutate(distance_to_center = distances2,
                                          mean_distance_to_center = mean(distances2), 
                                          qrtile = q_2)
  
  # cluster 3
  
  q_3 <- quantile(distances3, QR)
  indiv_clust03 <- indiv_clust3 %>% mutate(distance_to_center = distances3,
                                          mean_distance_to_center = mean(distances3),
                                          qrtile = q_3)
  
  # cluster 4
  
  q_4 <- quantile(distances4, QR)
  indiv_clust04 <- indiv_clust4 %>% mutate(distance_to_center = distances4,
                                          mean_distance_to_center = mean(distances4), 
                                          qrtile = q_4)
  
  # compile all data set 
  indiv_clust <- rbind(indiv_clust01, indiv_clust02) %>%
    rbind(., indiv_clust03)%>%
    rbind(., indiv_clust04)%>%
    select(Dim.1, distance_to_center, mean_distance_to_center, qrtile)
  
  indiv_clust <- rbind(indiv_clust01, indiv_clust02,indiv_clust03, indiv_clust04)%>%
           select(Dim.1, distance_to_center, mean_distance_to_center, qrtile)
  
  # add the distance_to_center lesser than qrtile 
  
  objs0 <- indiv_clust %>% mutate(keep = distance_to_center <= qrtile,Dim1=Dim.1)
  
  # select objs 
  objs0 <- objs0 %>% select(Dim1, distance_to_center, qrtile, keep)
  angola_detritus0 <- merge(angola_detritus, objs0, by = 'Dim1')
  # unique sample id 
  uni<-unique(angola_detritus0$sample_id)
  cycleNum<-cbind(uni,c(1:118))
  cyclNumframe<-data.frame(cycleNum)
  cyclNumframe<-cyclNumframe %>% rename(sample_id=uni,CYCLE_NUMBER=V2)
  # merge this id and number given to id called cycle id
  angola_detritus0 <- merge(angola_detritus0, cyclNumframe, by = 'sample_id')
  
  objs_keep <- angola_detritus0 %>% filter(keep == TRUE)
  objs_delete <- angola_detritus0 %>% filter(keep == FALSE)
  
  # plot the station in the space
    ggplot() +
    geom_point(aes(x = Dim1, y = Dim2), data  = objs_keep, color = as.factor(objs_keep$cluster_kmeans), size = 4) +
    #geom_point(aes(x = Dim1, y = Dim2), data  = objs_delete, color = as.factor(objs_delete$cluster_kmeans), alpha = 0.05, size = 4) +
    #geom_label(aes(x = Dim.1, y = Dim.2, label = CYCLE_NUMBER), data = objs_delete, color = as.factor(objs_delete$cluster), alpha = 0.1)+
    geom_point(data=centers, aes(Dim.1, Dim.2), alpha=0.5, shape=16, size = 4, color = 'yellow') + #add to add cluster centers
    coord_fixed() + scale_colour_discrete(guide="none") +
    scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0)
  ggsave(paste0('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/FigClusterqrtile',as.character(QR*100),'.jpeg'), width=16, height=6)
   
  ########################################################################################################################################
  ################################################################################
  #
  ################################################################################
  
  
  ggplot(angola_detritus0) +
    geom_point(aes(Dim1, Dim2, colour=as.factor(cluster_kmeans)), size = 2) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    scale_colour_discrete("Cluster", guide=guide_legend(override.aes=list(shape=16)))
  
  
  
  
  ## Display and analysis of clusters -----
  
  # get some example images from each reduced cluster
  clust_imgs <- angola_detritus0 |>
    group_by(cluster_kmeans) |>
    sample_n(size=min(n(), 100)) |>
    group_map(.f=function(x, y) {
      ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
        labs(title=str_c("Cluster ", y$cluster_kmeans)) +
        theme(plot.title=element_text(hjust=0.5))
    })
  wrap_plots(clust_imgs, ncol=2)
  
  # compare with only selected individuals (distance to center =< mean(distance to center)) 
  clust_imgs <- objs_keep |>
    group_by(cluster_kmeans) |>
    sample_n(size=min(n(), 100)) |>
    group_map(.f=function(x, y) {
      ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
        labs(title=str_c("Cluster ", y$cluster_kmeans)) +
        theme(plot.title=element_text(hjust=0.5))
    })
  wrap_plots(clust_imgs, ncol=2)
  
  
  # read csv of cluster and morpho arranged in python
  #Angola_Morphodetrituscluster_used <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/Angola_Morphodetrituscluster_used.csv", sep = ",", header = TRUE)
  #angola_detritus0<-angola_detritus
  #angola_detritus<-Angola_Morphodetrituscluster_used
  
  
  # Box plot of the size, shape, brightness and structure of each cluster 
  
  # size ----
  
  size <- objs_keep %>% select(perim., Cluster)
  size$perim. <- log10(size$perim.)
  
  cluster_colors <- c("cluster 1" = "#DE64B1",
                      "cluster 2" = "#57A7B3",
                      "cluster 3" = "#F5F28E",
                      "cluster 4" = "grey33",
                      "cluster 5" = "#A6DAE6")
  
  size$Cluster <- factor(size$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))
  
  p1 <- ggplot(size, aes(x = Cluster, y = perim., fill = Cluster))+
    geom_boxplot() +
    scale_fill_manual(values = cluster_colors) +
    ggtitle("Size") +
    xlab("Cluster Name") +
    ylab("log(Perimeter)") +
    theme_bw()+
    theme(legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))
  p1
  
  # shape ----
  
  
  shape <- objs_keep %>% select(circ., Cluster)
  shape$circ. <- shape$circ.
  
  shape$Cluster <- factor(shape$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))
  
  p2 <- ggplot(shape, aes(x = Cluster, y = circ., fill = Cluster))+
    geom_boxplot() +
    scale_fill_manual(values = cluster_colors) +
    ggtitle("Shape") +
    xlab("Cluster Name") +
    ylab("Circularity") +
    theme_bw()+
    theme(legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))
  
  p2
  
  # Brightness ----
  
  brightness <- objs_keep %>% select(mean, Cluster)
  
  brightness$Cluster <- factor(brightness$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))
  
  p3 <- ggplot(brightness, aes(x = Cluster, y = mean, fill = Cluster))+
    geom_boxplot() +
    scale_fill_manual(values = cluster_colors) +
    ggtitle("Brightness") +
    xlab("Cluster Name") +
    ylab("mean grey level") +
    theme_bw()+
    theme(legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))
  p3
  
  # Structure ----
  
  structure <- objs_keep %>% select(kurt, Cluster)
  
  structure$Cluster <- factor(structure$Cluster, levels = c("cluster 1", "cluster 2", "cluster 3", "cluster 4", "cluster 5"))
  
  p4 <- ggplot(structure, aes(x = Cluster, y = kurt, fill = Cluster))+
    geom_boxplot() +
    scale_fill_manual(values = cluster_colors) +
    ggtitle("Structure") +
    xlab("Cluster Name") +
    ylab("kurtosis") +
    theme_bw()+
    theme(legend.position = "none",
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 15))
  
  p4
  
  # put the four plot into one figure ----
  library(gridExtra)
  jpeg(paste0('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/FigClusterCharacteristicqrtile',as.character(QR*100),'.jpeg'), width = 10, height = 4, units = 'in', res = 600)
  grid.arrange(p1, p2, p3, p4, ncol = 4) # not so well 
  dev.off()
  # select only columns that I'm interested in 
  
  angola_detritus0$datetime <- paste(angola_detritus0$date, angola_detritus0$time, sep = "")
  #angola_detritus <- angola_detritus %>% select(-date, -time)
  #angola_detritus$datetime <- as.POSIXct(angola_detritus$datetime, format = "%Y%m%d%H%M%S")
  angola_detritus0$datetime <- as.POSIXct(angola_detritus0$datetime, format = "%Y%m%d")
  
  #ts_df <-  angola_detritus %>% select(id, datetime, depth_min, Dim1, Dim2, Dim3, Dim4, Cluster)
  ts_df <-  angola_detritus0 %>% select(id,sample_id, CYCLE_NUMBER,depth_min, Dim1, Dim2, Dim3, Dim4, Cluster,datetime,keep)
  
  
  # save it 
  
  write.csv(ts_df, paste0('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/k-means_pca_resultscENTER_Q', as.character(QR*100),'.csv'), row.names = FALSE)
  
}

#ggsave("cluster_distrib_season.pdf", width=8, height=6)
