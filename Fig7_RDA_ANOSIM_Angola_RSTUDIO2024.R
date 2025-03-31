
library("vegan")
library(ggplot2)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages('devEMF')
library('devEMF') 
#install.packages("ggpubr")
library(ggpubr)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

MDSSS_dble=read.table('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/MDS_dbleQ100rstudio.csv', sep = ",", header = TRUE)
MDSSS_dble=read.table('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/MDS_dbleQ100Rstudio_vfOCT2023.csv', sep = ",", header = TRUE)

MDS_Environment=read.table('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/MDS_Environment.csv', sep = ",", header = TRUE) 
#MDS_Environment$DEPTH<-NULL
DATEOK=read.table('C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6/MDS_dbleQ100EX.csv', sep = ",", header = TRUE)
dataok<-DATEOK$dateok
SHAP<-MDSSS_dble$date
SHAP[SHAP=='period1']<-'Export Event1'
SHAP[SHAP=='period2']<-'Export Event2'
SHAP[SHAP=='period3']<-'Export Event3'
SHAP[SHAP=='period4']<-'Export Event4'
SHAP[SHAP=='period5']<-'Export Event5'
SHAP[SHAP=='period6']<-'Export Event6'
SHAP[grepl('00:00',SHAP)]<-'No Export'

group<-MDSSS_dble$date
group[group=='period1']<-'Export Event1'
group[group=='period2']<-'Export Event2'
group[group=='period3']<-'Export Event3'
group[group=='period4']<-'Export Event4'
group[group=='period5']<-'Export Event5'
group[group=='period6']<-'Export Event6'
group[1:120]<-'Export Event1'
group[226:240]<-'Export Event3'
group[301:310]<-'Export Event4'
group[371:405]<-'Export Event5'
group[536:580]<-'Export Event7'

names(MDSSS_dble)[names(MDSSS_dble)=="cluster.1"]<-"Flakes"
names(MDSSS_dble)[names(MDSSS_dble)=="cluster.2"]<-"Agglomerates"
names(MDSSS_dble)[names(MDSSS_dble)=="cluster.3"]<-"Strings"
names(MDSSS_dble)[names(MDSSS_dble)=="cluster.4"]<-"Spheres"

profile<-unique(MDSSS_dble$profile)
ishap=0
LC=0
counter = 1
for (profi in profile)
{profileX<-grep(profi,MDSSS_dble$profile)
ishap[counter]<-profileX[1]
LC=cbind(LC, cbind(counter,counter, counter, counter, counter))
counter = counter + 1
}

LCx=LC[2:ncol(LC)]

qQR=c('0-100m','100-200m','200-300m','300-400m','400-500m')
qQR=c('00m')

for (QR in qQR) 
{
  se=MDSSS_dble$layer==QR
  se<-grep(QR,MDSSS_dble$layer)
  data=log10(MDSSS_dble[se,4:7]+1)
  here.mds <- metaMDS(comm = data, distance = "bray", trace = FALSE, autotransform = FALSE)
  #plot(here.mds$points)
  MDS_xy <- data.frame(here.mds$points)
  #MDS_xy$dateok<-MDSSS_dble[se,]$dateok
  MDS_xy$dateok<-1:nrow(MDS_xy)
  MDS_xy$dateok<- LCx[se]
  MDS_xy$layer<- MDSSS_dble[se,]$layer
  MDS_xy$SH<-SHAP[ishap]
  par(bg="white")
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, nrow(MDS_xy)))# , color=MDS_xy$SH
  ggplot() + geom_point(aes(x=MDS_xy$MDS1, y=MDS_xy$MDS2, shape=MDS_xy$layer, size=4)) + xlim(min(0,min(MDS_xy[,1])),min(quantile(MDS_xy[,1],0.95))) + sc + #+ theme_bw()
    labs(x="MDSx",y="MDSy", title="MDS all layers",shape="Layer") + #,color="EVENTS"
    theme(plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour="white") )
  #ggplot(MDS_xy, aes(x=MDS1, y=MDS2, color=dateok)) + geom_point(size=4) + xlim(min(0,min(MDS_xy[,1])),min(1/3,quantile(MDS_xy[,1],0.95))) + sc #+ theme_bw()
  #geom_point(aes(x=sitecoord[,1]*3,y=sitecoord[,2]*3,label=rownames(sitecoord),color=dataenv_all$ocean[1:120]),size=log10(S1$S1+1)*coeff*2/4) +
  
  
  #ggsave(paste0("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/MDS_Rstudio/MDS" , QR , ".jpeg"), width=8, height=6)
} 


qQR=c('0-100m','100-200m','200-300m','300-400m','400-500m')
ANOpvalues=0
ANOstat=0
count=0
for (QR in qQR) 
{
  count=count+1
  se=MDSSS_dble$layer==QR
  se<-grep(QR,MDSSS_dble$layer)
  data=log10(MDSSS_dble[se,4:7]+1)
  here.mds <- metaMDS(comm = data, distance = "bray", trace = FALSE, autotransform = FALSE)
  #plot(here.mds$points)
  MDS_xy <- data.frame(here.mds$points)
  #MDS_xy$dateok<-MDSSS_dble[se,]$dateok
  MDS_xy$dateok<-1:nrow(MDS_xy)
  MDS_xy$dateok<- LCx[se]
  MDS_xy$layer<- MDSSS_dble[se,]$layer
  MDS_xy$SH<-SHAP[ishap]
  par(bg="white")
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, nrow(MDS_xy)))# 
  ggplot() + geom_point(aes(x=MDS_xy$MDS1, y=MDS_xy$MDS2, shape=MDS_xy$layer, color=MDS_xy$SH, size=4)) + xlim(min(0,min(MDS_xy[,1])),min(quantile(MDS_xy[,1],0.95)))  + #sc + #+ theme_bw()
    scale_color_manual(values=c("red", "black","orange","green","cyan","blue","gray"))+
    labs(x="MDSx",y="MDSy", title="MDS all layers",color="EVENTS",shape="Layer") + #
    theme(plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white", colour="white") )
  #ggplot(MDS_xy, aes(x=MDS1, y=MDS2, color=dateok)) + geom_point(size=4) + xlim(min(0,min(MDS_xy[,1])),min(1/3,quantile(MDS_xy[,1],0.95))) + sc #+ theme_bw()
  #geom_point(aes(x=sitecoord[,1]*3,y=sitecoord[,2]*3,label=rownames(sitecoord),color=dataenv_all$ocean[1:120]),size=log10(S1$S1+1)*coeff*2/4) +
  ano = anosim(data, MDS_xy$SH, distance = "bray", permutations = 9999)
  ANOpvalues[count]=ano$signif
  ANOstat[count]=ano$statistic
  
  #ggsave(paste0("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/MDS_Rstudio/MDS" , QR , ".jpeg"), width=8, height=6)
} 


0/0
#STOP

'dataenv3_nona <- dataenv3
for(i in 1:ncol(dataenv2)){
  dataenv3_nona[is.na(dataenv3_nona[,i]), i] <- mean(dataenv3_nona[,i], na.rm = TRUE)
}'

qQR=c('0-100m','100-200m','200-300m','300-400m','400-500m','00m')
ANOpvalues=0
ANOstat=0
count=0
RDA_dble<-MDSSS_dble[c(1:564,566:580),]
for (QR in qQR) 
{
  #se=RDA_dble$layer==QR
  se<-grep(QR,RDA_dble$layer)
  ENV<-MDS_Environment[MDS_Environment$profile != '0020a_WMO6903096',]
  ENV<-ENV[ENV$profile != '0029a_WMO6903096',]
  databio_all<-RDA_dble[se,4:7]
  row.names(databio_all) <- 1:nrow(databio_all)#MDSSS_dble$profile
  # TRANFORMED DATA HELLINGER
  databio_all_t<- decostand(databio_all, "hellinger") 
  bio_all_t.rda <- rda(databio_all_t, ENV[se,c(11:16)], scale=TRUE) # c(5,11:16)
  #varsupp<-predict(bio_all_t.rda, type = "sp", newdata=ENV[se,c(17:18)], scaling = "none")
  
  bio_all_t.rda$CCA$eig
  bio_all_t.rda$CCA$u
  bio_all_t.rda$CCA$v
  eigenval <- bio_all_t.rda$CCA$eig   # Here you will get the eigenvalues
  sitecoord <- bio_all_t.rda$CCA$u   # The site coordinates 
  varcoord <-bio_all_t.rda$CCA$v    # The species coordinates 
  envcoord <-bio_all_t.rda$CCA$biplot  # The environmental variables coordinates 
  rownames(envcoord)
  envlab=c('Density','Temperature','Salinity','Oxygen','bbpPOC','Chla')
  
  eig <- data.frame(eigenval)
  eig$nb <- c(1:length(eigenval))
  eig$prop <- eig$eigenval/sum(eig$eigenval)
  eig
  
  h<-0.1/2   # 193:232
  h<-0.2/2   # 193:232
  hy<-0.1/2   # 193:232
  
  if (QR == '0-100m') {
    sitecoord[,1]=sitecoord[,1]
    varcoord[,1]=varcoord[,1]
    envcoord[,1]=envcoord[,1]
    }
    
  # Which are the significant eigenvalues according to the Kaiser-Guttman criteria?
  # 3 axes
  
  text1 <- paste("RDA1 (",floor(eig$prop[1]*10000)/100," %)", sep = "")
  text2 <- paste("RDA2 (",floor(eig$prop[2]*10000)/100," %)", sep = "")
  text3 <- paste("RDA3 (",floor(eig$prop[3]*10000)/100," %)", sep = "")
  text4 <- paste("RDA4 (",floor(eig$prop[4]*10000)/100," %)", sep = "")
  S1<-rowSums(databio_all)
  S1<-as.data.frame(S1)
  SHAPI=SHAP[se]
  groupi=group[se]
  # Plotting the results of the RDA 
  # By default:
  plot(bio_all_t.rda)
  
  
  library('ggrepel')
  
  # Nicer with ggplot:

  
  coeff<-1
  GGP<-ggplot() +
    #geom_point(aes(x=sitecoord[,1]*3,y=sitecoord[,2]*3,label=rownames(sitecoord),color=SHAPI,shape=ENV[se,]$layer)) +
    geom_point(shape=21, size = 5, aes(x=-sitecoord[,1]*3,y=sitecoord[,2]*3,color=groupi, shape=ENV[se,]$layer, fill=factor(SHAPI))) +
    scale_color_manual(values=c("red", "black","orange","#006164","cyan","#2c7be9","#2c7be9"))+
    #geom_point(aes(x=sitecoord[c(1:120,226:240,301:310,371:405,536:580),1]*3,y=sitecoord[c(1:120,226:240,301:310,371:405,536:580),2]*3,fill=factor(group))) +
    scale_fill_manual(values=c("red", "black","orange","#006164","cyan","#2c7be9","gray"))+
    
    #scale_color_manual(values=c("red", "black","orange","green","cyan","blue","gray"))+
    geom_segment(mapping=aes(x=0, y=0, xend=-varcoord[,1], yend=(varcoord[,2])), arrow=arrow(angle = 20,length = unit(0.1, "inches")),color="black")+
    geom_segment(mapping=aes(x=0, y=0, xend=-envcoord[,1], yend=(envcoord[,2])), arrow=arrow(angle = 20,length = unit(0.1, "inches")),color="blue")+
    geom_text_repel(aes(x=-(varcoord[,1]+h*sign(varcoord[,1])),y=(varcoord[,2]+(h+0.15*h)*sign(varcoord[,2]))),label=colnames(databio_all),color="black",size=6, force_pull = 0)+
    geom_text_repel(aes(x=-(envcoord[,1]+(h)*sign(envcoord[,1])),y=(envcoord[,2]+h*sign(envcoord[,2]))),label=envlab,color="blue",size=6, force_pull=0)+
    
    #geom_segment(mapping=aes(x=0, y=0, xend=varsupp[,1], yend=varsupp[,2]), arrow=arrow(angle = 20,length = unit(0.1, "inches")),color="blue")+
    #geom_text_repel(aes(x=varsupp[,1]+h*sign(varsupp[,1])+h,y=varsupp[,2]+h*sign(varsupp[,2])),label=rownames(varsupp),color="blue",size=3,segment.alpha=0.5)+
    
    theme_bw() +
    geom_vline(aes(xintercept=0), linetype="dashed")+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    
    labs(x=text1,y=text2, title=paste0("RDA ",QR, " using helliger-transformed abundances"),color="columns",fill="", shape="Layers") + 
    guides(color="none") # fill="Export Event"

  ggpar(GGP,
          legend = "right", legend.title = "Export Event") + 
      font("legend.title", size=12)
  
  GGP+theme(axis.text.x=element_text(size=20, angle=90),axis.text.y=element_text(size=20, angle=90),
            plot.title = element_text(size=20),
            axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            legend.text = element_text(size=12)) 
  GGP
  anova(bio_all_t.rda)
  RsquareAdj(bio_all_t.rda)
  
  #ggsave(paste0("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/MDS_Rstudio/RDAVFoct2023" , QR, ".jpeg"), width=8, height=6)
  ######     ggsave(paste0("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/MDS_Rstudio/Aout2024/RDAVFAout2024VS" , QR, ".jpeg"), dpi=1200, width=8, height=6)
  #ggsave(paste0("C:/Users/syawo/Dropbox/codes 2024/fig RDA Angola/RDAVFFEV2024 " , QR, ".jpeg"), width=8, height=6)
  ggsave(paste0("C:/Users/syawo/Downloads/Figures_pdf_Angola_geosciences_submissionOct20204/RDAVFAout2024VSZbon" , QR, ".jpeg"), dpi=1200, width=7.5, height=6.2)
} 


0/INF/0

############### SAME AS ABOVE BUT X-AXIS CHANGES SIGN but not used in the manuscript  ###########################################3

qQR=c('0-100m','100-200m','200-300m','300-400m','400-500m','00m') 
ANOpvalues=0
ANOstat=0
count=0
RDA_dble<-MDSSS_dble[c(1:564,566:580),]
for (QR in qQR) 
{
  #se=RDA_dble$layer==QR
  se<-grep(QR,RDA_dble$layer)
  ENV<-MDS_Environment[MDS_Environment$profile != '0020a_WMO6903096',]
  ENV<-ENV[ENV$profile != '0029a_WMO6903096',]
  databio_all<-RDA_dble[se,4:7]
  row.names(databio_all) <- 1:nrow(databio_all)#MDSSS_dble$profile
  # TRANFORMED DATA HELLINGER
  databio_all_t<- decostand(databio_all, "hellinger") 
  bio_all_t.rda <- rda(databio_all_t, ENV[se,c(11:16)], scale=TRUE) # c(5,11:16)
  #varsupp<-predict(bio_all_t.rda, type = "sp", newdata=ENV[se,c(17:18)], scaling = "none")
  
  bio_all_t.rda$CCA$eig
  bio_all_t.rda$CCA$u
  bio_all_t.rda$CCA$v
  eigenval <- bio_all_t.rda$CCA$eig   # Here you will get the eigenvalues
  sitecoord <- bio_all_t.rda$CCA$u   # The site coordinates 
  varcoord <-bio_all_t.rda$CCA$v    # The species coordinates 
  envcoord <-bio_all_t.rda$CCA$biplot  # The environmental variables coordinates 
  rownames(envcoord)
  envlab=c('Density','Temperature','Salinity','Oxygen','bbpPOC','Chla')
  
  eig <- data.frame(eigenval)
  eig$nb <- c(1:length(eigenval))
  eig$prop <- eig$eigenval/sum(eig$eigenval)
  eig
  
  h<-0.1/2   # 193:232
  h<-0.2/2   # 193:232
  hy<-0.15/2   # 193:232
  hx<-0/2   # 193:232
  
  if (QR == '0-100m') {
    sitecoord[,1]=sitecoord[,1]
    varcoord[,1]=varcoord[,1]
    envcoord[,1]=envcoord[,1]
  }
  
  # Which are the significant eigenvalues according to the Kaiser-Guttman criteria?
  # 3 axes
  
  text1 <- paste("RDA1 (",floor(eig$prop[1]*10000)/100," %)", sep = "")
  text2 <- paste("RDA2 (",floor(eig$prop[2]*10000)/100," %)", sep = "")
  text3 <- paste("RDA3 (",floor(eig$prop[3]*10000)/100," %)", sep = "")
  text4 <- paste("RDA4 (",floor(eig$prop[4]*10000)/100," %)", sep = "")
  S1<-rowSums(databio_all)
  S1<-as.data.frame(S1)
  SHAPI=SHAP[se]
  groupi=group[se]
  # Plotting the results of the RDA 
  # By default:
  plot(bio_all_t.rda)
  
  
  library('ggrepel')
  
  # Nicer with ggplot:
  
  
  coeff<-1
  GGP<-ggplot() +
    #geom_point(aes(x=sitecoord[,1]*3,y=sitecoord[,2]*3,label=rownames(sitecoord),color=SHAPI,shape=ENV[se,]$layer)) +
    geom_point(shape=21, size = 3, aes(x=sitecoord[,1]*3,y=sitecoord[,2]*3,color=groupi, shape=ENV[se,]$layer, fill=factor(SHAPI))) +
    scale_color_manual(values=c("#EE6677", "black","orange","green","cyan","blue","blue"))+
    #geom_point(aes(x=sitecoord[c(1:120,226:240,301:310,371:405,536:580),1]*3,y=sitecoord[c(1:120,226:240,301:310,371:405,536:580),2]*3,fill=factor(group))) +
    scale_fill_manual(values=c("#EE6677", "black","orange","green","cyan","blue","gray"))+
    
    
    
    #scale_color_manual(values=c("red", "black","orange","green","cyan","blue","gray"))+
    geom_segment(mapping=aes(x=0, y=0, xend=varcoord[,1], yend=varcoord[,2]), arrow=arrow(angle = 20,length = unit(0.1, "inches")),color="black")+
    geom_segment(mapping=aes(x=0, y=0, xend=envcoord[,1], yend=envcoord[,2]), arrow=arrow(angle = 20,length = unit(0.1, "inches")),color="blue")+
    geom_text_repel(aes(x=(varcoord[,1]+h*sign(varcoord[,1])),y=varcoord[,2]+h*sign(varcoord[,2]),label=colnames(databio_all)),color="black",size=6)+
    geom_text_repel(aes(x=(envcoord[,1]+h*sign(envcoord[,1])+hx),y=envcoord[,2]+hy*sign(envcoord[,2])),label=envlab,color="blue",size=6)+
    
    #geom_segment(mapping=aes(x=0, y=0, xend=varsupp[,1], yend=varsupp[,2]), arrow=arrow(angle = 20,length = unit(0.1, "inches")),color="blue")+
    #geom_text_repel(aes(x=varsupp[,1]+h*sign(varsupp[,1])+h,y=varsupp[,2]+h*sign(varsupp[,2])),label=rownames(varsupp),color="blue",size=3,segment.alpha=0.5)+
    
    theme_bw() +
    geom_vline(aes(xintercept=0), linetype="dashed")+
    geom_hline(aes(yintercept=0), linetype="dashed")+
    
    labs(x=text1,y=text2, title=paste0("RDA ",QR, " using helliger-transformed abundances"),color="columns",fill="", shape="Layers") + 
    guides(color="none") # fill="Export Event"
  
  ggpar(GGP,
        legend = "right", legend.title = "Export Event") + 
    font("legend.title", size=12)
  
  GGP+theme(axis.text.x=element_text(size=20, angle=90),axis.text.y=element_text(size=20, angle=90),
            plot.title = element_text(size=20),
            axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            legend.text = element_text(size=12)) 
  GGP
  
  count=count+1
  eval(str2expression(paste0(paste0('Rsquare', count), ' = RsquareAdj(bio_all_t.rda)')))
  eval(str2expression(paste0(paste0('Anova', count), ' = anova(bio_all_t.rda)')))

  #ggsave(paste0("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/MDS_Rstudio/RDAVFoct2023" , QR, ".jpeg"), width=8, height=6)
  ######ggsave(paste0("C:/MOPGA 2022-2023 Dodji/ARGO float Angola/EcotaxaPart_BCG_UVP6_Plot/MDS_Rstudio/Aout2024/Positifxaxis/RDAVFAout2024" , QR, ".jpeg"), dpi=1200, width=8, height=6)
  #ggsave(paste0("C:/Users/syawo/Dropbox/codes 2024/fig RDA Angola/RDAVFFEV2024 " , QR, ".jpeg"), width=8, height=6)
  
} 


