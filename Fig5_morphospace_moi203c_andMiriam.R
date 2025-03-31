
# Script morphospace  in (C:\Users\syawo\Downloads)

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

# Function to pre-process UVP images
preprocess <- function(x) {
  x |>
    # remove 31 pixels from the bottom (=the scale bar)
    img_chop(bottom=31) |>
    # change the gamma to see the light objects better
    img_adjust_gamma(gamma=0.7)
}

# Here i put the csv i sent to Miriam
#angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Allemande_bureau/etape0/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
#angola_detritus <- read.table("D:/AngolaRSTUDIO_MOPGA2023/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
# remove the "object_" part of column name 
names(angola_detritus) <- gsub("object_", "", names(angola_detritus), fixed = TRUE)

# load ms and centers and objs here before continue ...........................
ms<-readRDS("C:/Users/syawo/Downloads/ms_PCA_wosize.RData")

centers_PCA_wosize_kmeans10 <- read_csv("C:/Disque/all téléchargements/Télé/centers_PCA_wosize_kmeans10.csv")
#
centers_PCA_wosize_kmeans10 <- read_csv("C:/Disque/all téléchargements/Télé/centers_PCA_wosize_kmeans5.csv")
#
centers_PCA_wosize_kmeans10 <- read_delim("C:/Disque/all téléchargements/Télé/centers_umap_wosize_kmeans5.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
 

objs_PCA_wosize <- read_csv("C:/Disque/all téléchargements/Télé/objs_PCA_wosize.csv")

objs_PCA_wosize <- read_csv("C:/Disque/all téléchargements/Télé/objs_umap_wosize.csv")
# visual representation of the morphospace
ggmorph_tile(ms, objs_PCA_wosize$img_path,  steps=18, n_imgs=3, fun=preprocess)
#ggmorph_tile(ms, angola_detritus$img_path,  steps=18, n_imgs=3, fun=preprocess)

PCs <- 1:2
# get coordinates of active variables/columns
var <- ms$var$coord[,PCs] |> as.data.frame() |> rownames_to_column()
ggmorph_tile(ms, angola_detritus$img_path, steps=18, n_imgs=3, fun=preprocess) +
  #geom_point(data=objs, aes(Dim.1, Dim.2), shape=".", alpha=0.5) +
  geom_segment(data=var, aes(x = 0, y = 0, xend = (Dim.1*5),
                             yend = (Dim.2*5)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black") + geom_point(size = 3) + annotate("text", x = (var$Dim.1*5), y = (var$Dim.2*5),
           label = var$rowname)
ggtitle("Title")



pca.vars <- ms$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)

ggmorph_tile(ms, angola_detritus$img_path, steps=18, n_imgs=3, fun=preprocess) +
  geom_segment(data=pca.vars, aes(x = 0, xend = Dim.1*6, y = 0, yend = Dim.2*6),
               arrow=arrow(length = unit(0.025, "npc"), type = "open"),
               lwd=1, color="seagreen") +
  geom_text(data=pca.vars, aes(x=Dim.1*6.15, y=Dim.2*6.15, label=vars),
            color="seagreen", size=4, hjust=0) #+ # labels


#objs <- ms$ind$coord |> as_tibble()

# hierarchisation of the morphs
clust_h <- hclust(dist(centers_PCA_wosize_kmeans10), method="ward.D2")
plot(clust_h, hang=-1)
#clust_reduced <- cutree(clust_h, k=10) |> factor()

clust_reduced <- cutree(clust_h, k=5) |> factor()
#clust_reduced <- cutree(clust_h, k=4) |> factor()

# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs_PCA_wosize$umap_wosize_kmeans5), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centers_PCA_wosize_kmeans10)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

ggplot(objs_PCA_wosize) +
  geom_point(aes(Dim.1, Dim.2, colour=angola_detritus$cluster_reduced), size = 2) +
  coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
  scale_colour_discrete("Cluster", guide=guide_legend(override.aes=list(shape=16)))

## Display and analysis of clusters -----

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
       group_by(cluster_reduced) |>
       sample_n(size=min(n(), 20)) |>
       group_map(.f=function(x, y) {
           ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
               labs(title=str_c("Cluster ", y$cluster_reduced)) +
               theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=5)
0/0

#UMAP WOIZE K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Here i put the csv i sent to Miriam
#angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Allemande_bureau/etape0/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
# remove the "object_" part of column name 
names(angola_detritus) <- gsub("object_", "", names(angola_detritus), fixed = TRUE)

centers_UMAP_wosize_kmeans4 <- read_delim("C:/Disque/all téléchargements/Télé/centers_umap_wosize_kmeans4.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
objs_UMAP_wosize <- read_csv("C:/Disque/all téléchargements/Télé/objs_UMAP_wosize.csv")

# hierarchisation of the morphs
clust_h <- hclust(dist(centers_UMAP_wosize_kmeans4), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()

# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs_UMAP_wosize$umap_wosize_kmeans4), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centers_UMAP_wosize_kmeans4)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters -----

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=4)


# PCA WOIZE K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Here i put the csv i sent to Miriam
#angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Allemande_bureau/etape0/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)

# remove the "object_" part of column name 
names(angola_detritus) <- gsub("object_", "", names(angola_detritus), fixed = TRUE)

centers_PCA_wosize_kmeans4 <- read_csv("C:/Disque/all téléchargements/Télé/centers_PCA_wosize_kmeans4.csv")
objs_PCA_wosize <- read_csv("C:/Disque/all téléchargements/Télé/objs_PCA_wosize.csv")
objs_PCA_wosize$cluster_wosize_kmeans4[objs_PCA_wosize$cluster_wosize_kmeans4==2]<-22
objs_PCA_wosize$cluster_wosize_kmeans4[objs_PCA_wosize$cluster_wosize_kmeans4==3]<-2
objs_PCA_wosize$cluster_wosize_kmeans4[objs_PCA_wosize$cluster_wosize_kmeans4==22]<-3

# hierarchisation of the morphs
clust_h <- hclust(dist(centers_PCA_wosize_kmeans4), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()

# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs_PCA_wosize$cluster_wosize_kmeans4), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centers_PCA_wosize_kmeans4)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters -----

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=4)


#UMAP K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Here i put the csv i sent to Miriam
#angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Allemande_bureau/etape0/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)

angola_detritus <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/Angola_detritus_forMiriam.csv", sep = ",", header = TRUE)
# remove the "object_" part of column name 
names(angola_detritus) <- gsub("object_", "", names(angola_detritus), fixed = TRUE)

centers_UMAP_kmeans4 <- read_delim("C:/Disque/all téléchargements/Télé/centers_umap_kmeans4.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
centers_UMAP_kmeans4 <- read_delim("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/centers_morpho/centers_umap_kmeans4.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
objs_UMAPP <- read_csv("C:/Disque/all téléchargements/Télé/objs_umap.csv")
objs <- read.table("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/Angola_detritus_allclustersMiriamSept.csv", sep = ",", header = TRUE)

# hierarchisation of the morphs
clust_h <- hclust(dist(centers_UMAP_kmeans4), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()

# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs$umap_kmeans4), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centers_UMAP_kmeans4)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters ----- 

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=4)

#PCA K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centersK4 <- read_delim("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/centers_morpho/centers_PCA_kmeans4.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
# hierarchisation of the morphs
clust_h <- hclust(dist(centersK4), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()
# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs$cluster_kmeans4), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centersK4)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters ----- 

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=4)

#PCA K5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centersK5 <- read_delim("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/centers_morpho/centers_PCA_kmeans5.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
# hierarchisation of the morphs
clust_h <- hclust(dist(centersK5), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()
# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs$cluster_kmeans5), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centersK5)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters ----- 

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=5)

#PCA K10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centersK10 <- read_delim("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/centers_morpho/centers_PCA_kmeans10.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
# hierarchisation of the morphs
clust_h <- hclust(dist(centersK10), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()
# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs$cluster_kmeans10), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centersK10)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters ----- 

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=4)

#PCA WOsIZE K4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centersWOK4 <- read_delim("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/aMiriamAout2023/centers_morpho/centers_PCA_wosize_kmeans4.csv", delim = ",", escape_double = FALSE, trim_ws = TRUE)
# hierarchisation of the morphs
clust_h <- hclust(dist(centersWOK4), method="ward.D2")
plot(clust_h, hang=-1)
clust_reduced <- cutree(clust_h, k=4) |> factor()
# assign each original point to a reduced cluster
angola_detritus <- mutate(angola_detritus, cluster_kmeans=factor(objs$cluster_wosize_kmeans4), cluster_reduced=NA)
match_clusters <- tibble(cluster_kmeans=factor(1:nrow(centersWOK4)),cluster_reduced=clust_reduced)
angola_detritus <- left_join(select(angola_detritus, -cluster_reduced), match_clusters)

## Display and analysis of clusters ----- 

# get some example images from each reduced cluster
clust_imgs <- angola_detritus |>
  group_by(cluster_reduced) |>
  sample_n(size=min(n(), 100)) |>
  group_map(.f=function(x, y) {
    ggimg_grid(x$img_path, fun=preprocess, scale=0.002) +
      labs(title=str_c("Cluster ", y$cluster_reduced)) +
      theme(plot.title=element_text(hjust=0.5)) })
library("patchwork")
wrap_plots(clust_imgs, ncol=4)

