###############################################################
# ENVdata
#Get topoclimatic data from different sources                            #
# extract values at sampling points                           #
# merge everything and save as Envdata_XX.Rds                          #
#     #
###############################################################
#rm(list=ls())
library(terra)
library(tidyverse)

ENVstack<-rast("../spatial_data/ENVstack.tif")
dataSoil<-readRDS("../spatial_data/MpAlps_soil_data_VVerdon1903+.Rds")
names(dataSoil)[names(dataSoil)=="pH"]<-"pH_onsite"

for (group in c("BA","FU","PR")){#group="BA"
  print(paste0("group ",group))
  readcounts<-readRDS(file=paste0("data\\readcounts_",group,"_Guex.Rds"))
  edaph_data<-dataSoil[dataSoil[,paste0("sampleName",group)]%in% colnames(readcounts),]
  from_raster <- terra::extract(ENVstack, edaph_data[,c("x","y")])
  ENVdatatemp <- data.frame(from_raster[,-1],
                               edaph_data[,13:(ncol(edaph_data)-4)],
                               AllSand=edaph_data$ThickSand+edaph_data$ThinSand,
                               AllSilt=edaph_data$ThickSilt+edaph_data$ThinSilt)
  ENVdatatemp <- apply(ENVdatatemp,2,function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}) #z-score normalization (applied to all, even non normal variables)
  #Phyllosilicates + Quartz Feldspath_K Plagioclase_Na    Calcite Goethite   Indoses : not very informative (except calcite but correlated), win more plots.
  ENVdatatemp <- ENVdatatemp[,-which(colnames(ENVdatatemp)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]
  #remove HI OI cause too many empty data
  ENVdatatemp <- ENVdatatemp[,-which(colnames(ENVdatatemp)%in%c("HI","OI"))]
  ENVdatatemp <- ENVdatatemp[,-which(colnames(ENVdatatemp)%in%c("Limestones_SandyLimestones_.Marl_Shale","AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]
  # summary(ENVdatatemp)
  # par(mfrow=c(1,2)) #check per covariate
  # for (i in 1:ncol(ENVdatatemp)){
  #   hist(ENVdatatemp[,i],main=colnames(ENVdatatemp)[i])
  #   plot(ENVdatatemp[order(ENVdatatemp[,i]),i],main=colnames(ENVdatatemp)[i])
  #   invisible(readline(prompt="Press [enter] to continue"))
  # }
  # par(mfrow=c(1,0))
  ##not enough variability (not enough hit points), glacier is meh
  ENVdatatemp <- ENVdatatemp[,-which(colnames(ENVdatatemp)%in%c("agriculture_25","open_forest_25","overgrown_shrubland_unproductive_vegetation_25","permanent_crops_25"))]
  # cor(ENVdatatemp_BA[,"forest_25"],ENVdatatemp_BA[,"forestaggr_25"])#same variable
  # cor(ENVdatatemp_BA[,"pH_modeled"],ENVdatatemp_BA[,"pH_onsite"])#remove mdelled, but might replace for predictions
  # #remove pH modeled, might replace on ste pH if
  ENVdatatemp <- ENVdatatemp[,-which(colnames(ENVdatatemp)%in%c("closed_forest_25","forestaggr_25","pH_modeled"))]
  
  rownames(ENVdatatemp) <- edaph_data[,paste0("sampleName",group)]
  #Remove sites with NA
  #nrow(ENVdatatemp_BA)
  print(paste0("nsite_total: ",nrow(ENVdatatemp)))
  ENVdatatemp <- ENVdatatemp[complete.cases(ENVdatatemp),] #264 -> 250
  print(paste0("nsite_complete.case: ",nrow(ENVdatatemp)))
    # factorization of remaining binary var, doing it elsewhere to not create problem in cor matrix
    # binvar <- c()
    # for (i in colnames(ENVdatatemp)){#i =names(ENVdata) [31]
    #   if(length(unique(ENVdatatemp[,i]))==2){
    #     ENVdatatemp[,i]<-as.factor(ENVdatatemp[,i])
    #     binvar <- c(binvar,i)
    #   }
    # }
  # summary(ENVdatatemp_BA)
  # nrow(ENVdatatemp_BA) 250
  print(paste0("n_cov: ",ncol(ENVdatatemp)))
  # print(paste0("n_bin_cov: ",ncol(ENVdatatemp)))
  # ncol(ENVdatatemp_BA) 79 dont 12 bin
  ENVdata<-ENVdatatemp
  
  saveRDS(ENVdata,file=paste0("data/ENVdata_",group,".Rds"))
}






#old version (paper 1)
################################################################################
# readcounts_BA<-readRDS(file="data\\readcounts_BA_Guex.Rds")
# # colnames(OTU_data_bac_265samples)
# readcounts_FU<-readRDS(file="data\\readcounts_FU_Guex.Rds")
# readcounts_PR<-readRDS(file="data\\readcounts_PR_Guex.Rds")
# 
# 
# # saveRDS(dataSoil,file="../spatial_data/MpAlps_soil_data_VVerdon1903+.Rds")# !PATCH careful if rewrite dataSoil
# 
# edaph_data_BA <- dataSoil[dataSoil$sampleNameBA %in% colnames(readcounts_BA), ] # pick only sites for which there is BA OTUdata
# 
# edaph_data_FU <- dataSoil[dataSoil$sampleNameFU %in% colnames(readcounts_FU), ] # pick only sites for which there is FU OTUdata
# 
# edaph_data_PR <- dataSoil[dataSoil$sampleNamePR %in% colnames(readcounts_PR), ] # pick only sites for which there is PR OTUdata
# 
# # edaph_data_BA$altitude
# 
# from_raster_BA <- terra::extract(ENVstack, edaph_data_BA[,c("x","y")])
# from_raster_FU <- terra::extract(ENVstack, edaph_data_FU[,c("x","y")])
# from_raster_PR <- terra::extract(ENVstack, edaph_data_PR[,c("x","y")])
# 
# # plot(ENVstack$bio1_t)
# # points(edaph_data_BA[which(edaph_data_BA$sampleNameBA%in%rownames(OTUdata_BA_PA)),c("x","y")],col="red")
# # points(edaph_data_FU[which(edaph_data_FU$sampleNameFU%in%rownames(OTUdata_BA_FU)),c("x","y")],col="red")
# # points(edaph_data_PR[which(edaph_data_PR$sampleNamePR%in%rownames(OTUdata_BA_PR)),c("x","y")],col="red")
# 
# # ENVdatatemp_BA <- data.frame(from_raster_BA[,colnames(from_raster_BA)!="d13C"], #for stack with d13C that double the var
# #                              edaph_data_BA[,c(13:ncol(edaph_data_BA))],
# #                              AllSand=edaph_data_BA$ThickSand+edaph_data_BA$ThinSand, 
# #                              AllSilt=edaph_data_BA$ThickSilt+edaph_data_BA$ThinSilt)
# # ENVdatatemp_FU <- data.frame(from_raster_FU[,colnames(from_raster_FU)!="d13C"], 
# #                              edaph_data_FU[,c(13:ncol(edaph_data_FU))],
# #                              AllSand=edaph_data_FU$ThickSand+edaph_data_FU$ThinSand, 
# #                              AllSilt=edaph_data_FU$ThickSilt+edaph_data_FU$ThinSilt)
# # ENVdatatemp_PR <- data.frame(from_raster_PR[,colnames(from_raster_PR)!="d13C"], 
# #                              edaph_data_PR[,c(13:ncol(edaph_data_PR))],
# #                              AllSand=edaph_data_PR$ThickSand+edaph_data_PR$ThinSand, 
#                              # AllSilt=edaph_data_PR$ThickSilt+edaph_data_PR$ThinSilt)
# ENVdatatemp_BA <- data.frame(from_raster_BA[,-1],
#                              edaph_data_BA[,13:(ncol(edaph_data_BA)-4)],
#                              AllSand=edaph_data_BA$ThickSand+edaph_data_BA$ThinSand,
#                              AllSilt=edaph_data_BA$ThickSilt+edaph_data_BA$ThinSilt)
# ENVdatatemp_BA <- apply(ENVdatatemp_BA,2,function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}) #z-score normalization (applied to all, even non normal variables)
# 
# 
# 
# # summary(ENVdatatemp_BA) ; nrow(ENVdatatemp_BA)
# #remove some var 
# #c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt")
# #dry and wet because they are just here to have water content (multicol)
# #EC because we have 2 ways of measuring the same thing
# #AllSand and Allsilt because they are sum of more precise ones already in the table
# #removed Ankerite, because "04_Clust_UnivMod.R" cant fit univariate gam model for OTU15660 for that variable. And its not an important variable.
# #removed Dolomite, because not enough values != from 0 to perform variable selection with cubic regression.
# #(later classified in "trash")
# # varcor <- cor(ENVdatatemp_BA,use = "pairwise.complete.obs")
# # plot(hclust(as.dist(1-abs(var_cor)), "average"))
# # abline(h=0.3,lty=2, col="red", lwd=2)
# # abline(h=0.2,lty=2, col="black", lwd=2)
# # abline(h=0.25,lty=2, col="blue", lwd=2)
# 
# #Tmax removed, not interpretable and missing data
# #C.N removed : Obvious multicolinearity (literally Carbon / Nitrogen) and missing data
# #TOC removed : missing data and highly correlated with Nitrogen/Carbo /BSWC group
# #MINC removed : missing data and highly correlated with CaO/pH group
# #Phyllosilicates + Quartz Feldspath_K Plagioclase_Na    Calcite Goethite   Indoses : not very informative (except calcite but correlated), win more plots.
# ENVdatatemp_BA <- ENVdatatemp_BA[,-which(colnames(ENVdatatemp_BA)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]
# 
# #remove HI OI cause too many empty data
# ENVdatatemp_BA <- ENVdatatemp_BA[,-which(colnames(ENVdatatemp_BA)%in%c("HI","OI"))]
# 
# ENVdatatemp_BA <- ENVdatatemp_BA[,-which(colnames(ENVdatatemp_BA)%in%c("Limestones_SandyLimestones_.Marl_Shale","AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]
# # summary(ENVdatatemp_BA)
# # par(mfrow=c(1,2))
# # for (i in 1:ncol(ENVdatatemp_BA)){
# #   hist(ENVdatatemp_BA[,i],main=colnames(ENVdatatemp_BA)[i])
# #   plot(ENVdatatemp_BA[order(ENVdatatemp_BA[,i]),i],main=colnames(ENVdatatemp_BA)[i])
# #   invisible(readline(prompt="Press [enter] to continue"))
# # }
# # par(mfrow=c(1,0))
# 
# ##not enough variability (not enough hit points), glacier is meh
# ENVdatatemp_BA <- ENVdatatemp_BA[,-which(colnames(ENVdatatemp_BA)%in%c("agriculture_25","open_forest_25","overgrown_shrubland_unproductive_vegetation_25","permanent_crops_25"))]
# # cor(ENVdatatemp_BA[,"forest_25"],ENVdatatemp_BA[,"forestaggr_25"])#same variable
# # cor(ENVdatatemp_BA[,"pH_modeled"],ENVdatatemp_BA[,"pH_onsite"])#remove mdelled, but might replace for predictions
# # #remove pH modeled, might replace on ste pH if
# ENVdatatemp_BA <- ENVdatatemp_BA[,-which(colnames(ENVdatatemp_BA)%in%c("closed_forest_25","forestaggr_25","pH_modeled"))]
# 
# rownames(ENVdatatemp_BA) <- edaph_data_BA$sampleNameBA
# #Remove sites with NA
# #nrow(ENVdatatemp_BA)
# ENVdatatemp_BA <- ENVdatatemp_BA[complete.cases(ENVdatatemp_BA),] #264 -> 250
# # summary(ENVdatatemp_BA) ; nrow(ENVdatatemp_BA)
# 
# #remove binary var with not enough points with the var
# 
# #factorization of remaining binary var
# binvar <- c()
# for (i in colnames(ENVdatatemp_BA)){#i =names(ENVdata) [31]
#   if(length(unique(ENVdatatemp_BA[,i]))==2){
#     ENVdatatemp_BA[,i]<-as.factor(ENVdatatemp_BA[,i])
#     binvar <- c(binvar,i)
#   }
# }
# # summary(ENVdatatemp_BA)
# # nrow(ENVdatatemp_BA) 250
# # ncol(ENVdatatemp_BA) 79 dont 12 bin
# ENVdata_BA <-ENVdatatemp_BA
# 
# 
# 
# ENVdatatemp_FU <- data.frame(from_raster_FU[,-1], 
#                              edaph_data_FU[,13:(ncol(edaph_data_FU)-4)], 
#                              AllSand=edaph_data_FU$ThickSand+edaph_data_FU$ThinSand, 
#                              AllSilt=edaph_data_FU$ThickSilt+edaph_data_FU$ThinSilt)
# ENVdatatemp_FU <- apply(ENVdatatemp_FU,2,function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)})
# 
# ENVdatatemp_FU <- ENVdatatemp_FU[,-which(colnames(ENVdatatemp_FU)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]
# 
# #remove HI OI cause too many empty data
# ENVdatatemp_FU <- ENVdatatemp_FU[,-which(colnames(ENVdatatemp_FU)%in%c("HI","OI"))]
# 
# ENVdatatemp_FU <- ENVdatatemp_FU[,-which(colnames(ENVdatatemp_FU)%in%c("Limestones_SandyLimestones_.Marl_Shale","AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]
# # summary(ENVdatatemp_FU)
# # par(mfrow=c(1,2))
# # for (i in 1:ncol(ENVdatatemp_FU)){
# #   hist(ENVdatatemp_FU[,i],main=colnames(ENVdatatemp_FU)[i])
# #   plot(ENVdatatemp_FU[order(ENVdatatemp_FU[,i]),i],main=colnames(ENVdatatemp_FU)[i])
# #   invisible(readline(prompt="Press [enter] to continue"))
# # }
# # par(mfrow=c(1,0))
# 
# ##not enough variability (not enough hit points), glacier is meh
# ENVdatatemp_FU <- ENVdatatemp_FU[,-which(colnames(ENVdatatemp_FU)%in%c("agriculture_25","open_forest_25","overgrown_shrubland_unproductive_vegetation_25","permanent_crops_25"))]
# # cor(ENVdatatemp_FU[,"forest_25"],ENVdatatemp_FU[,"forestaggr_25"])#same variable
# # cor(ENVdatatemp_FU[,"pH_modeled"],ENVdatatemp_FU[,"pH_onsite"])#remove mdelled, but might replace for predictions
# # #remove pH modeled, might replace on ste pH if
# ENVdatatemp_FU <- ENVdatatemp_FU[,-which(colnames(ENVdatatemp_FU)%in%c("closed_forest_25","forestaggr_25","pH_modeled"))]
# 
# rownames(ENVdatatemp_FU) <- edaph_data_FU$sampleNameFU
# #Remove sites with NA
# #nrow(ENVdatatemp_FU)
# ENVdatatemp_FU <- ENVdatatemp_FU[complete.cases(ENVdatatemp_FU),] #232 -> 217
# # summary(ENVdatatemp_FU) ; nrow(ENVdatatemp_FU)
# 
# #remove binary var with not enough points with the var
# 
# #factorization of remaining binary var
# binvar <- c()
# for (i in colnames(ENVdatatemp_FU)){#i =names(ENVdata) [31]
#   if(length(unique(ENVdatatemp_FU[,i]))==2){
#     ENVdatatemp_FU[,i]<-as.factor(ENVdatatemp_FU[,i])
#     binvar <- c(binvar,i)
#   }
# }
# # summary(ENVdatatemp_FU)
# # nrow(ENVdatatemp_FU) 217
# # ncol(ENVdatatemp_FU) 79 dont 12 bin
# ENVdata_FU <-ENVdatatemp_FU
# 
# 
# ENVdatatemp_PR <- data.frame(from_raster_PR[,-1],
#                              edaph_data_PR[,13:(ncol(edaph_data_PR)-4)],
#                              AllSand=edaph_data_PR$ThickSand+edaph_data_PR$ThinSand,
#                              AllSilt=edaph_data_PR$ThickSilt+edaph_data_PR$ThinSilt)
# ENVdatatemp_PR <- apply(ENVdatatemp_PR,2,function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}) #z-score normalization (applied to all, even non normal variables)
# 
# 
# 
# # summary(ENVdatatemp_PR) ; nrow(ENVdatatemp_PR)
# #remove some var 
# #c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt")
# #dry and wet because they are just here to have water content (multicol)
# #EC because we have 2 ways of measuring the same thing
# #AllSand and Allsilt because they are sum of more precise ones already in the table
# #removed Ankerite, because "04_Clust_UnivMod.R" cant fit univariate gam model for OTU15660 for that variable. And its not an important variable.
# #removed Dolomite, because not enough values != from 0 to perform variable selection with cubic regression.
# #(later classified in "trash")
# # varcor <- cor(ENVdatatemp_PR,use = "pairwise.complete.obs")
# # plot(hclust(as.dist(1-abs(var_cor)), "average"))
# # abline(h=0.3,lty=2, col="red", lwd=2)
# # abline(h=0.2,lty=2, col="black", lwd=2)
# # abline(h=0.25,lty=2, col="blue", lwd=2)
# 
# #Tmax removed, not interpretable and missing data
# #C.N removed : Obvious multicolinearity (literally Carbon / Nitrogen) and missing data
# #TOC removed : missing data and highly correlated with Nitrogen/Carbo /BSWC group
# #MINC removed : missing data and highly correlated with CaO/pH group
# #Phyllosilicates + Quartz Feldspath_K Plagioclase_Na    Calcite Goethite   Indoses : not very informative (except calcite but correlated), win more plots.
# ENVdatatemp_PR <- ENVdatatemp_PR[,-which(colnames(ENVdatatemp_PR)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]
# 
# #remove HI OI cause too many empty data
# ENVdatatemp_PR <- ENVdatatemp_PR[,-which(colnames(ENVdatatemp_PR)%in%c("HI","OI"))]
# 
# ENVdatatemp_PR <- ENVdatatemp_PR[,-which(colnames(ENVdatatemp_PR)%in%c("Limestones_SandyLimestones_.Marl_Shale","AngularGravel_Blocks_SlopeScree","Conglomerates_weakmed_solidified_with_sandstone_and_deposits","Dolomites_Cornelia","GravelSand_currentdeposit","Dolomite","Goethite","Ankerite"))]
# # summary(ENVdatatemp_PR)
# # par(mfrow=c(1,2))
# # for (i in 1:ncol(ENVdatatemp_PR)){
# #   hist(ENVdatatemp_PR[,i],main=colnames(ENVdatatemp_PR)[i])
# #   plot(ENVdatatemp_PR[order(ENVdatatemp_PR[,i]),i],main=colnames(ENVdatatemp_PR)[i])
# #   invisible(readline(prompt="Press [enter] to continue"))
# # }
# # par(mfrow=c(1,0))
# 
# ##not enough variability (not enough hit points), glacier is meh
# ENVdatatemp_PR <- ENVdatatemp_PR[,-which(colnames(ENVdatatemp_PR)%in%c("agriculture_25","open_forest_25","overgrown_shrubland_unproductive_vegetation_25","permanent_crops_25"))]
# # cor(ENVdatatemp_PR[,"forest_25"],ENVdatatemp_PR[,"forestaggr_25"])#same variable
# # cor(ENVdatatemp_PR[,"forest_25"],ENVdatatemp_PR[,"closed_forest_25"])
# # cor(ENVdatatemp_PR[,"pH_modeled"],ENVdatatemp_PR[,"pH_onsite"])#remove mdelled, but might replace for predictions
# # plot(ENVdatatemp_PR[,"forest_25"],main=colnames(ENVdatatemp_PR)[i])
# # #remove pH modeled, might replace on ste pH if
# ENVdatatemp_PR <- ENVdatatemp_PR[,-which(colnames(ENVdatatemp_PR)%in%c("closed_forest_25","forestaggr_25","pH_modeled"))]
# 
# rownames(ENVdatatemp_PR) <- edaph_data_PR$sampleNamePR
# #Remove sites with NA
# #nrow(ENVdatatemp_PR)
# ENVdatatemp_PR <- ENVdatatemp_PR[complete.cases(ENVdatatemp_PR),] #174 -> 166
# # summary(ENVdatatemp_PR) ; nrow(ENVdatatemp_PR)
# 
# #remove binary var with not enough points with the var
# 
# #factorization of remaining binary var
# binvar <- c()
# for (i in colnames(ENVdatatemp_PR)){#i =names(ENVdata) [31]
#   if(length(unique(ENVdatatemp_PR[,i]))==2){
#     ENVdatatemp_PR[,i]<-as.factor(ENVdatatemp_PR[,i])
#     binvar <- c(binvar,i)
#   }
# }
# # summary(ENVdatatemp_PR)
# # nrow(ENVdatatemp_PR) 166
# # ncol(ENVdatatemp_PR) 79 dont 12 bin
# ENVdata_PR <-ENVdatatemp_PR
# 
# 
# 
# saveRDS(ENVdata_BA,file="data/ENVdata_BA.Rds")
# saveRDS(ENVdata_FU,file="data/ENVdata_FU.Rds")
# saveRDS(ENVdata_PR,file="data/ENVdata_PR.Rds")
