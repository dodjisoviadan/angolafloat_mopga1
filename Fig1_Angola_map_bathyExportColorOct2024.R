#Prerequisites
#You are going to need to install a few extra packages to follow along with this lecture.

# these are packages you will need, but probably already have.
# Don't bother installing if you already have them

# install.packages(c("ggplot2", "devtools", "dplyr", "stringr"))

# some standard map packages.
#install.packages(c("maps", "mapdata"))
#install.packages("ggOceanMaps")

# the github version of ggmap, which recently pulled in a small fix I had
# for a bug 
# devtools::install_github("dkahle/ggmap")
#Load up a few of the libraries we will use
library(ggplot2)
library(ggpubr)
library(ggmap)
library(maps)
library(mapdata)
library(vegan)
#Plotting maps-package maps with ggplot

library(ggOceanMaps)
#install.packages("marmap")
library(marmap)

# install.packages('dplyr')
library('dplyr')
#install.packages("ggnewscale")
library(ggnewscale)

#devtools::install_github("eliocamp/metR")
library(metR)
#install.packages('geomtextpath')
library(geomtextpath)

#install.packages('devEMF')
library('devEMF')

datageo=read.table('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Remi_Laxenaire_dataAngola/PostionAngolaclosestdistanceKm-Rstudio.csv',sep=';',header=TRUE)

#datageo=read.table('C:/Users/syawo/Downloads/Stage L3 UK no to send - but selection for students/WP2_200µm_stageNini/EnvprovVFin_andmaps.csv',sep=',',header=TRUE)

### Load the global coastline using map_data, you will overlay it on the bathymetry data which are negative (remove altitude on land cells)
world <- map_data('world',xlim=c(10, 15), ylim=c(-15, -10))
# get bathymetry data
b = getNOAA.bathy(lon1 = 9, lon2 = 16, lat1 = -16, lat2 = -9, resolution = 1) # 15
bathy <- fortify(b) # convert to ddf
class(bathy)
summary(bathy) ; head(bathy)

#palette = "Greys"

#emf(file = "C:/MOPGA 2022-2023 Dodji/ARGO float Angola/Remi_Laxenaire_dataAngola/Angola_MAPplotEMF.emf", emfPlus = TRUE)

databaty=bathy[bathy$z <= 0,]

#scale_fill_manual
map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = z), data = bathy[bathy$z <= 0,]) + 
  scale_fill_distiller(name = "Depth (m)", limits=c(-7000,0)) + 
  geom_polygon(data=world, aes(x=long, y=lat, group = group), fill = "grey75", colour = "black", size = 0.3)+
  theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey75")) + #,linetype = "dashed"
  
  new_scale('fill') +
  #ggplot_add.new_aes()+
  geom_point(shape=21,size=8, data=datageo, aes(x= as.numeric(Longitude) , y= as.numeric(Latitude), color=factor(ExportEventNofill), fill=factor(ExportEvent))) +
  scale_color_manual(values=c("red", "black","orange","#006164","cyan","#2c7be9","#2c7be9"))+ # "#EE6677"
  scale_fill_manual (values=c("red", "black","orange","#006164","cyan","#2c7be9","gray60"))+
  
  geom_contour(data = bathy[bathy$z <= 0,], aes(x = x, y = y, z = z), breaks = c(-200, -500, -1000,-2000, -3000), color = "grey70", size = 0.2) +
  
  #geom_contour(data=bathy[bathy$z <= 0,], aes(x = x, y = y, z = z), binwidth = 500, color = "black", size = 0.2) +
  
  #geom_text_contour(data=bathy[bathy$z <= 0,],aes(x = x, y = y, z = z), binwidth = 200,check_overlap=TRUE) +
  geom_textcontour(data=bathy[bathy$z <= 0,],aes(x = x, y = y, z = z), breaks = c(-200, -500, -1000,-2000,-3000),  size=6,straight=TRUE,text_smoothing=1, remove_long =TRUE ,halign='center', text_only=TRUE, hjust=0.5,vjust=0.5) + # 
  
  #new_scale("fill") +
  xlab("Longitude (°)" ) + ylab("Latitude (°)") + 

  #guides(fill=guide_legend())+
  labs(fill="Export Events")  +
  coord_sf(xlim=c(10, 15), ylim=c(-15,-10)) +
  annotate('text', y=-10.70,x=12.73000,label=expression('t'[0]),size=8)+
  annotate('text', y=-11.35,x=11.21500,label=expression('t'[f]),size=8)+
  guides(color="none") 


map1 +
  theme(axis.text.x=element_text(size=20, angle=90),axis.text.y=element_text(size=20, angle=90),
          plot.title = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20),
          legend.text = element_text(size=12)) 

map1$labels$fill<-"Export Events"
map1$labels$colour<-'Export Events'

ggsave('C:/Users/syawo/Downloads/Figures_pdf_Angola_geosciences_submissionOct20204/map1bonbon.pdf', width = 72*15/72, height = 72*8/72, dpi = 1200, device='pdf')#, units = "in"

#dev.off()

mm<-map1+annotate('text', y=-10.70,x=12.73000,label=expression('t'[0]),size=8)
mmm<-mm + annotate('text', y=-11.35,x=11.21500,label=expression('t'[f]),size=8)
mmm +
  theme(axis.text.x=element_text(size=20, angle=90),axis.text.y=element_text(size=20, angle=90),
        plot.title = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=12)) 

map1$labels$fill<-"Export Events"
map1$labels$colour<-'Export Events'

stop<- STOP #'do not use below now' #############################################################################################################################################

#geom_contour(data = bathy, 
             #aes(x = x, y = y, z = z),
             #binwidth = 200, color = "black", size = 0.1) +


library(ggOceanMaps)
#install.packages('ggspatial')
library(ggspatial) # for data plotting
basemap(limits = 60) # A synonym: basemap(60) 
basemap(limits = c(10, 15, -15, -10), bathymetry = TRUE,bathy.style = "poly_blues")

BB=basemap(limits = c(10, 15, -15, -10), bathy.style = "rcb") # synonym to "raster_continuous_blues"

BB

scale_color_manual(values=c("#EE6677", "black","orange","green","cyan","blue","blue","#7b3294"))
                            
######################################################################################" 

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

ggplot_add.new_aes <- function(object, plot, object_name) {
    plot$layers <- lapply(plot$layers, bump_aes, new_aes = object)
    plot$scales$scales <- lapply(plot$scales$scales, bump_aes, new_aes = object)
    plot$labels <- bump_aes(plot$labels, new_aes = object)
    plot
  }

