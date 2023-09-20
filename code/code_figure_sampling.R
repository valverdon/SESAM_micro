
#plot points

library(tidyverse)
library(terra)
hillshade <- rast("../spatial_data/envlayers_terra/hillshade.tif")
dataSoil<-readRDS("../spatial_data/MpAlps_soil_data_VVerdon1903+.Rds")
OTUdata_BA<-readRDS("data/OTUdata_BA.Rds")
OTUdata_FU<-readRDS("data/OTUdata_FU.Rds")
OTUdata_PR<-readRDS("data/OTUdata_PR.Rds")

edaph_data_BA2<-dataSoil[dataSoil$sampleNameBA %in% rownames(OTUdata_BA),]
edaph_data_FU2<-dataSoil[dataSoil$sampleNameFU %in% rownames(OTUdata_FU),]
edaph_data_PR2<-dataSoil[dataSoil$sampleNamePR %in% rownames(OTUdata_PR),]

#BA points that are in FU and in PR
BA_FU_PR<-edaph_data_BA2[edaph_data_BA2$x%in%edaph_data_FU2$x&edaph_data_BA2$y%in%edaph_data_FU2$y&edaph_data_BA2$x%in%edaph_data_PR2$x&edaph_data_BA2$y%in%edaph_data_PR2$y,c("x","y")]
#Ba points that are in FU but not in PR
BA_FU<-edaph_data_BA2[edaph_data_BA2$x%in%edaph_data_FU2$x&edaph_data_BA2$y%in%edaph_data_FU2$y&(!(edaph_data_BA2$x%in%edaph_data_PR2$x))&(!(edaph_data_BA2$y%in%edaph_data_PR2$y)),c("x","y")]
#BA points that are in Pr but not in FU
BA_PR<-edaph_data_BA2[edaph_data_BA2$x%in%edaph_data_PR2$x&edaph_data_BA2$y%in%edaph_data_PR2$y&(!(edaph_data_BA2$x%in%edaph_data_FU2$x))&(!(edaph_data_BA2$y%in%edaph_data_FU2$y)),c("x","y")]
#Ba points that are neither in FU or PR
BA_only<-edaph_data_BA2[(!(edaph_data_BA2$x%in%edaph_data_PR2$x))&(!(edaph_data_BA2$y%in%edaph_data_PR2$y))&(!(edaph_data_BA2$x%in%edaph_data_FU2$x))&(!(edaph_data_BA2$y%in%edaph_data_FU2$y)),c("x","y")]

library(viridis)
viridis(4)
plot(hillshade,col=gray.colors(1000),legend=FALSE, axes=FALSE)
points(BA_FU_PR,pch=21,lwd=1,col="black",bg="#FDE725FF",cex=1)
points(BA_FU,pch=22,lwd=1,col="black",bg="#35B779FF",cex=1)
points(BA_PR,pch=23,col="black",bg="#31688EFF",cex=1)
points(BA_only,pch=24,col="black",bg="#440154FF",cex=1)

# points(edaph_data_BA[,c("x","y")],pch=16,col="#aec800",cex=1.2)
# points(edaph_data_FU[,c("x","y")],pch=17,col=alpha("#6B4C62",0.8))
# points(edaph_data_PR[,c("x","y")],pch=18,col=alpha("#457EB0",0.7))
legend(2555000,1130000,c("BAFP","BAF", "BAP","BA"),pch=c(21,22,23,24), col=rep("black",4),pt.bg=c("#FDE725FF","#35B779FF","#31688EFF","#440154FF"),bty="n")


png(file=paste0("figures/PAAB_selection/Datasets_sampling.png"),res=300,width=1961,height=1500)
plot(hillshade,col=gray.colors(1000),legend=FALSE, axes=FALSE)
points(BA_FU_PR,pch=21,lwd=1,col="black",bg="#FDE725FF",cex=.8)
points(BA_FU,pch=22,lwd=1,col="black",bg="#35B779FF",cex=.8)
points(BA_PR,pch=23,col="black",bg="#31688EFF",cex=.8)
points(BA_only,pch=24,col="black",bg="#440154FF",cex=.8)
legend(2552000,1130000,c("BAFP","BAF", "BAP","BA"),cex=.8,pch=c(21,22,23,24), col=rep("black",4),pt.bg=c("#FDE725FF","#35B779FF","#31688EFF","#440154FF"),bty="n")

dev.off()
