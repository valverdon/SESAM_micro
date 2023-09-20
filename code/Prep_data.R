###############################################################
# ENVdata
#Get topoclimatic data from different sources                            #
# extract values at sampling points                           #
# merge everything in Envdata_full                            #
#     #
###############################################################
#rm(list=ls())
library(raster)
library(tidyverse)
##################################### climatic and topographic variables ####################################################

load("../../spatial_data/ENVstack.Rda")
# names(ENVstack)
# plot(ENVstack[["Altitude"]])
getwd()
#problem : 5 points out of Altitude layer

########################################## edaphic variables #####################################################################
#point data
# load OTU data from 265 sites in addition to soil data to find matching sites
load("../../ASV_data\\BA_OTU_265samples.Rdata")
# colnames(OTU_data_bac_265samples)
load("../../ASV_data\\FU_OTU_232samples.Rdata")
load("../../ASV_data\\PR_OTU_179samples.Rdata")
PRtaxofull<-read.csv("../../../30_data/ASVtables_taxo/Tax_table_Protists.csv",sep=",")
Onlypro <- OTU_data_pro_179samples[rownames(OTU_data_pro_179samples)%in%PRtaxofull$X,]

load("../../spatial_data/MpAlps_soil_data_VVerdon1903+.Rda")
edaph_data_BA <- dataSoil[dataSoil$sampleNameBA %in% colnames(OTU_data_bac_265samples), ] # pick only sites for which there is BA OTUdata
edaph_data_FU <- dataSoil[dataSoil$sampleNameFU %in% colnames(OTU_data_fun_232samples), ] # pick only sites for which there is FU OTUdata
edaph_data_PR <- dataSoil[dataSoil$sampleNamePR %in% colnames(Onlypro), ] # pick only sites for which there is PR OTUdata

from_raster_BA <- terra::extract(ENVstack, edaph_data_BA[,c("x","y")])
from_raster_FU <- terra::extract(ENVstack, edaph_data_FU[,c("x","y")])
from_raster_PR <- terra::extract(ENVstack, edaph_data_PR[,c("x","y")])

ENVdatatemp_BA <- data.frame(from_raster_BA, 
                             edaph_data_BA[,13:ncol(edaph_data_BA)],
                             AllSand=edaph_data_BA$ThickSand+edaph_data_BA$ThinSand, 
                             AllSilt=edaph_data_BA$ThickSilt+edaph_data_BA$ThinSilt)
ENVdatatemp_FU <- data.frame(from_raster_FU, 
                             edaph_data_FU[,13:ncol(edaph_data_FU)], 
                             AllSand=edaph_data_FU$ThickSand+edaph_data_FU$ThinSand, 
                             AllSilt=edaph_data_FU$ThickSilt+edaph_data_FU$ThinSilt)
ENVdatatemp_PR <- data.frame(from_raster_PR, 
                             edaph_data_PR[,13:ncol(edaph_data_PR)], 
                             AllSand=edaph_data_PR$ThickSand+edaph_data_PR$ThinSand, 
                             AllSilt=edaph_data_PR$ThickSilt+edaph_data_PR$ThinSilt)


#Binarization of mostly 0 data.
#plot(ENVdatatemp_BA$MINC)
#list=c("Ankerite","Calcite","Dolomite","Feldspath_K","Goethite")

ENVdatatemp_BA[c("Ankerite","Calcite","Dolomite","Feldspath_K","Goethite")] <- as.factor(ENVdatatemp_BA[c("Ankerite","Calcite","Dolomite","Feldspath_K","Goethite")]!=0)

#summary(ENVdatatemp_BA )



#remove some var 
#c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt")
#dry and wet because they are just here to have water content (multicol)
#EC because we have 2 ways of measuring the same thing
#AllSand and Allsilt because they are sum of more precise ones already in the table
#removed Ankerite, because "04_Clust_UnivMod.R" cant fit univariate gam model for OTU15660 for that variable. And its not an important variable.
#removed Dolomite, because not enough values != from 0 to perform variable selection with cubic regression.
#(later classified in "trash")
# varcor <- cor(ENVdatatemp_BA,use = "pairwise.complete.obs")
# plot(hclust(as.dist(1-abs(varcor)), "average"))
# abline(h=0.2,lty=2, col="red", lwd=2)
# 

#Tmax removed, not interpretable and missing data
#C.N removed : Obvious multicolinearity (literally Carbon / Nitrogen) and missing data
#TOC removed : missing data and highly correlated with Nitrogen/Carbo /BSWC group
#MINC removed : missing data and highly correlated with CaO/pH group
#Phyllosilicates + Quartz Feldspath_K Plagioclase_Na    Calcite Goethite   Indoses : not very informative (except calcite but correlated), win more plots.
ENVdatatemp_BA <- ENVdatatemp_BA[,-which(names(ENVdatatemp_BA)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]

#remove HI OI cause too many empty data
ENVdatatemp_BA <- ENVdatatemp_BA[,-which(names(ENVdatatemp_BA)%in%c("HI","OI"))]

# ## test remove SITES with rock types NA
# coordinates(edaph_data_BA) <- edaph_data_BA[c("x","y")]
# napoints<-edaph_data_BA[rownames(ENVdatatemp_BA)%in%rownames(ENVdatatemp_BA)[which(is.na(ENVdatatemp_BA[,"Calcite"]))],]
# plot(ENVstack[["Altitude"]])
# plot(napoints,add=TRUE,col="red",pch=1,cex=0.5)
# coordinates(edaph_data_FU) <- edaph_data_FU[c("x","y")]
# napoints<-edaph_data_FU[rownames(ENVdatatemp_FU)%in%rownames(ENVdatatemp_FU)[which(is.na(ENVdatatemp_FU[,"Calcite"]))],]
# plot(napoints,add=TRUE,col="green",pch=1,cex=1)
# coordinates(edaph_data_PR) <- edaph_data_PR[c("x","y")]
# napoints<-edaph_data_PR[rownames(ENVdatatemp_PR)%in%rownames(ENVdatatemp_PR)[which(is.na(ENVdatatemp_PR[,"Calcite"]))],]
# plot(napoints,add=TRUE,col="blue",pch=1,cex=1.5)
# ##no real problem to remove them
# ## ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("Phyllosilicates","Quartz","Feldspath_K","Plagioclase_Na","Calcite","Goethite","Indoses"))]
# 
# ##idem soilTemp
# napoints<-edaph_data_BA[rownames(ENVdatatemp_BA)%in%rownames(ENVdatatemp_BA)[which(is.na(ENVdatatemp_BA[,"soilTemp"]))],]
# plot(ENVstack[["Altitude"]])
# plot(napoints,add=TRUE,col="red",pch=1,cex=0.5)
# napoints<-edaph_data_FU[rownames(ENVdatatemp_FU)%in%rownames(ENVdatatemp_FU)[which(is.na(ENVdatatemp_FU[,"soilTemp"]))],]
# plot(napoints,add=TRUE,col="green",pch=1,cex=1)
# napoints<-edaph_data_PR[rownames(ENVdatatemp_PR)%in%rownames(ENVdatatemp_PR)[which(is.na(ENVdatatemp_PR[,"soilTemp"]))],]
# plot(napoints,add=TRUE,col="blue",pch=1,cex=1.5)
# ##no real problem to remove them
# 
# ##idem chemistry
# napoints<-edaph_data_BA[rownames(ENVdatatemp_BA)%in%rownames(ENVdatatemp_BA)[which(is.na(ENVdatatemp_BA[,"SiO2"]))],]
# plot(ENVstack[["Altitude"]])
# plot(napoints,add=TRUE,col="red",pch=1,cex=0.5)
# napoints<-edaph_data_FU[rownames(ENVdatatemp_FU)%in%rownames(ENVdatatemp_FU)[which(is.na(ENVdatatemp_FU[,"SiO2"]))],]
# plot(napoints,add=TRUE,col="green",pch=1,cex=1)
# napoints<-edaph_data_PR[rownames(ENVdatatemp_PR)%in%rownames(ENVdatatemp_PR)[which(is.na(ENVdatatemp_PR[,"SiO2"]))],]
# plot(napoints,add=TRUE,col="blue",pch=1,cex=1.5)
# ##no real problem to remove them

# summary(ENVdatatemp_FU[!complete.cases(ENVdatatemp_FU),])
# summary(ENVdatatemp_PR[!complete.cases(ENVdatatemp_PR),])

ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]
ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("HI","OI"))]
# ENVdatatemp_FU <- ENVdatatemp_FU[,-which(names(ENVdatatemp_FU)%in%c("Phyllosilicates","Quartz","Feldspath_K","Plagioclase_Na","Calcite","Goethite","Indoses"))]
ENVdatatemp_PR <- ENVdatatemp_PR[,-which(names(ENVdatatemp_PR)%in%c("drySoilWeight","wetSoilWeight","EC_1_1","AllSand","AllSilt","Tmax","C.N","TOC","MINC"))]
ENVdatatemp_PR <- ENVdatatemp_PR[,-which(names(ENVdatatemp_PR)%in%c("HI","OI"))]
# ENVdatatemp_PR <- ENVdatatemp_PR[,-which(names(ENVdatatemp_PR)%in%c("Phyllosilicates","Quartz","Feldspath_K","Plagioclase_Na","Calcite","Goethite","Indoses"))]

rownames(ENVdatatemp_BA) <- edaph_data_BA$sampleNameBA
rownames(ENVdatatemp_FU) <- edaph_data_FU$sampleNameFU
rownames(ENVdatatemp_PR) <- edaph_data_PR$sampleNamePR


# summary(ENVdatatemp_BA) # chemistry -3 ; d13Cd15N -1 site
# summary(ENVdatatemp_FU) # PL layers - 7 sites ; hillshade -10 sites ; SCD -5 ; ; Altitude -2 ; chemistry -2 or -7
# summary(ENVdatatemp_PR) # PL layers - 11 sites ; hillshade -14 sites ; SCD -10 ; Altitude -5 ; chemistry -2 or -7

#Remove sites with NA
#nrow(ENVdatatemp_PR)
#nrow(ENVdata_PR)
ENVdata_BA <- ENVdatatemp_BA[complete.cases(ENVdatatemp_BA),] #264 -> 250
ENVdata_FU <- ENVdatatemp_FU[complete.cases(ENVdatatemp_FU),] #232 -> 217
ENVdata_PR <- ENVdatatemp_PR[complete.cases(ENVdatatemp_PR),] #174 -> 166
# summary(ENVdata_BA)


#ASV data preparation
#Bacteria
OTUdata_BA_AB <- OTU_data_bac_265samples[,which(colnames(OTU_data_bac_265samples)%in%dataSoil$sampleNameBA)]
OTUdata_BA_AB <- t(OTUdata_BA_AB)
# nrow(OTUdata_BA_AB) #264
# ncol(OTUdata_BA_AB) #60567 ASV

#remove plots that dont have measures for all variables
OTUdata_BA_AB <- OTUdata_BA_AB[complete.cases(ENVdatatemp_BA), ] 
# nrow(OTUdata_BA_AB) #250

OTUdata_BA_AB <- OTUdata_BA_AB[,colSums(OTUdata_BA_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them
#ncol(OTUdata_BA_AB) #58462

#Need treshold to remove very high and very low prevalence (algos wont be able to adjust to data)
#5 or 10% ?
numberpresence_BA <- apply(OTUdata_BA_AB,2,function(X){sum(X!=0)})
plot(numberpresence_BA,pch=20,cex=.5, main="number of plots in which the ASV \n is present")
abline(h=nrow(OTUdata_BA_AB)/10,col="blue")#10% cut
abline(h=nrow(OTUdata_BA_AB)/20,col="red") #5% cut
abline(h=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/10),col="blue")#10% cut
abline(h=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/20),col="red") #5% cut
plot(density(numberpresence_BA))
abline(v=nrow(OTUdata_BA_AB)/10,col="blue")#10% cut
abline(v=nrow(OTUdata_BA_AB)/20,col="red")#5% cut
abline(v=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/10),col="blue")#10% cut
abline(v=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/20),col="red") #5% cut



#remove less prevalent ones for AB
OTUdata_BA_AB <- OTUdata_BA_AB[,numberpresence_BA>=nrow(OTUdata_BA_AB)/20]


#transformation AB --> PA
OTUdata_BA_PA <- apply(OTUdata_BA_AB,c(1,2),FUN=function(X){X>0}) 


seqvec <- colnames(OTUdata_BA_AB)
names(seqvec) <- paste0("OTU",1:ncol(OTUdata_BA_AB))
save(seqvec,file="seqvecBA.Rda")
colnames(OTUdata_BA_AB) <- names(seqvec)



# ncol(OTUdata_BA_PA) #55649 ASV (removed if prevalence <5%)
# ncol(OTUdata_BA_AB) #55649 ASV (removed if prevalence <5%)
#253arrays of 220 ASV
TotSeqSum_BA <- rowSums(OTUdata_BA_AB)


#Fungi
OTUdata_FU_AB <- OTU_data_fun_232samples[,which(colnames(OTU_data_fun_232samples)%in%dataSoil$sampleNameFU)]
OTUdata_FU_AB <- t(OTUdata_FU_AB)
# nrow(OTUdata_FU_AB) #232
#ncol(OTUdata_FU_AB) #95387
OTUdata_FU_AB <- OTUdata_FU_AB[complete.cases(ENVdatatemp_FU), ] # #remove plots that dont have measures for all variables
# nrow(OTUdata_FU_AB) #217
OTUdata_FU_AB <- OTUdata_FU_AB[,colSums(OTUdata_FU_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them

# ncol(OTUdata_FU_AB) #92753

numberpresence_FU <- apply(OTUdata_FU_AB,2,function(X){sum(X!=0)})
plot(numberpresence_FU,pch=20,cex=.5, main="number of plots in which the ASV \n is present (only ASV with >100 reads) ")
abline(h=nrow(OTUdata_FU_AB)/10,col="blue")#10% cut 21
abline(h=nrow(OTUdata_FU_AB)/20,col="red") #5% cut 10
# abline(h=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/10),col="blue")#10% cut
# abline(h=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/20),col="red") #5% cut
plot(density(numberpresence_FU))
abline(v=nrow(OTUdata_FU_AB)/10,col="blue")#10% cut
abline(v=nrow(OTUdata_FU_AB)/20,col="red")#5% cut
# abline(v=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/10),col="blue")#10% cut
# abline(v=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/20),col="red") #5% cut

# OTUdata_FU_PA <- OTUdata_FU_PA[,numberpresence_FU_PA>=nrow(OTUdata_FU_PA)/20]

#remove less prevalent ones for AB
OTUdata_FU_AB <- OTUdata_FU_AB[,numberpresence_FU>=nrow(OTUdata_FU_AB)/20]


OTUdata_FU_PA <- apply(OTUdata_FU_AB,c(1,2),FUN=function(X){X>0})


seqvec <- colnames(OTUdata_FU_AB)
names(seqvec) <- paste0("OTU",1:ncol(OTUdata_FU_AB))
save(seqvec,file="seqvecFU.Rda")
colnames(OTUdata_FU_AB) <- names(seqvec)

# ncol(OTUdata_FU_AB) #51863
# ncol(OTUdata_FU_PA) #51863

TotSeqSum_FU <- rowSums(OTUdata_FU_AB)




#Protists
#nrow(Onlypro) # 3419
OTUdata_PR_AB <- Onlypro[,which(colnames(Onlypro)%in%dataSoil$sampleNamePR)]
OTUdata_PR_AB <- t(OTUdata_PR_AB)
# nrow(OTUdata_PR_AB) #174
# ncol(OTUdata_PR_AB) #3419
# TotSeqSum_PR_AB <- rowSums(OTUdata_PR) #sequencing depth for "precleaned" data
OTUdata_PR_AB <- OTUdata_PR_AB[complete.cases(ENVdatatemp_PR), ] # #remove plots that dont have measures for all variables

# nrow(OTUdata_PR_AB) #166
OTUdata_PR_AB <- OTUdata_PR_AB[,colSums(OTUdata_PR_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them

# ncol(OTUdata_PR_AB) #3323
numberpresence_PR <- apply(OTUdata_PR_AB,2,function(X){sum(X!=0)})
plot(numberpresence_PR,pch=20,cex=.5, main="number of plots in which the ASV \n is present (only ASV with >100 reads) ")
abline(h=nrow(OTUdata_PR_AB)/10,col="blue")#10% cut
abline(h=nrow(OTUdata_PR_AB)/20,col="red") #5% cut
# abline(h=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/10),col="blue")#10% cut
# abline(h=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20),col="red") #5% cut
plot(density(numberpresence_PR))
abline(v=nrow(OTUdata_PR_AB)/10,col="blue")#10% cut
abline(v=nrow(OTUdata_PR_AB)/20,col="red")#5% cut
# abline(v=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/10),col="blue")#10% cut
# abline(v=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20),col="red") #5% cut

#remove less prevalent ones
OTUdata_PR_AB <- OTUdata_PR_AB[,numberpresence_PR>=nrow(OTUdata_PR_AB)/20]

OTUdata_PR_PA <- apply(OTUdata_PR_AB,c(1,2),FUN=function(X){X>0})
# ncol(OTUdata_PR_AB) #3161
# ncol(OTUdata_PR_PA) #3161

seqvec <- colnames(OTUdata_PR_AB)
names(seqvec) <- paste0("OTU",1:ncol(OTUdata_PR_AB))
save(seqvec,file="../../seqvecPR.Rda")
colnames(OTUdata_PR_AB) <- names(seqvec)

#cut at 10% : BA = 26, FU = 23, PR = 17 
#cut at 5% : BA = 13 FU = 11 PR = 8
#Where to cut here? 
#10% removes most of fungi (most of them are present in a very few plots)
#5% leaves some protists with only 9 plots to fit the models
TotSeqSum_PR <- rowSums(OTUdata_PR_AB)





#save all as "OTUdata"
OTUdata <- OTUdata_BA_PA
ENVdata <- ENVdata_BA
save(OTUdata,file="BA/Outputs/OTUdata.Rda")
save(ENVdata,file="BA/Outputs/ENVdata.Rda")
OTUdata <- OTUdata_BA_AB
TotSeqSum <- TotSeqSum_BA
save(OTUdata,file="../Abundance/BA/Outputs/OTUdata.Rda")
save(ENVdata,file="../Abundance/BA/Outputs/ENVdata.Rda")
save(TotSeqSum,file="../Abundance/BA/Outputs/dataTotSeqSum.Rda")


OTUdata <- OTUdata_FU_PA
ENVdata <- ENVdata_FU
save(OTUdata,file="FU/Outputs/OTUdata.Rda")
save(ENVdata,file="FU/Outputs/ENVdata.Rda")
OTUdata <- OTUdata_FU_AB
TotSeqSum <- TotSeqSum_FU
save(OTUdata,file="../Abundance/FU/Outputs/OTUdata.Rda")
save(ENVdata,file="../Abundance/FU/Outputs/ENVdata.Rda")
save(TotSeqSum,file="../Abundance/FU/Outputs/dataTotSeqSum.Rda")

OTUdata <- OTUdata_PR_PA
ENVdata <- ENVdata_PR
save(OTUdata,file="PA/PR/Outputs/OTUdata.Rda")
save(ENVdata,file="PA/PR/Outputs/ENVdata.Rda")
OTUdata <- OTUdata_PR_AB
TotSeqSum <- TotSeqSum_PR
save(OTUdata,file="Abundance/PR/Outputs/OTUdata.Rda")
save(ENVdata,file="Abundance/PR/Outputs/ENVdata.Rda")
save(TotSeqSum,file="Abundance/PR/Outputs/dataTotSeqSum.Rda")

