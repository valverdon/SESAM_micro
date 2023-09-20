library(tidyverse)

for(group in c("BA","FU","PR")){#group="BA"
  print(paste0("group: ",group))
  readcounts<-readRDS(file=paste0("data\\readcounts_",group,"_Guex.Rds"))
  ENVdata<-readRDS(file=paste0("data\\ENVdata_",group,".Rds"))
  taxo<-readRDS(file=paste0("../../30_data/ASVtables_taxo/",group,"taxo2023_full.Rds"))
  OTUdata<- readcounts[,which(colnames(readcounts)%in%rownames(ENVdata))]
  # OTUdata[5600:5610,]
  TotSeqSum<-colSums(OTUdata) #only useful for abundance data to get rel. abundance of reads
  
  nOTU_prefilter<-nrow(OTUdata)
  print(paste0("nOTU_prefilter: ", nOTU_prefilter))

  OTUdata_taxo<-cbind(OTUdata ,taxo[match(rownames(OTUdata),taxo[,"Seq"]),])
  
  print(paste0("taxo filter"))

  OTUdata <- t(OTUdata)
  # OTUdata[,5600:5610]
  if (group=="BA"){
    print(paste0("n_unassigned(removed): ",sum(OTUdata_taxo$Kingdom=="unclassified_Root"))) #surement different FU PR #8186
    OTUdata <- OTUdata[,OTUdata_taxo$Kingdom!="unclassified_Root"] #surement different FU PR
    OTUdata_taxo<- OTUdata_taxo[OTUdata_taxo$Kingdom!="unclassified_Root",] #surement different FU PR
  }
  if(group=="FU"){
    print(paste0("n_unassigned(removed): ",sum(OTUdata_taxo$Kingdom!="Fungi"))) #surement different FU PR #8186
    OTUdata <- OTUdata[,OTUdata_taxo$Kingdom=="Fungi"] #surement different FU PR
    OTUdata_taxo<- OTUdata_taxo[OTUdata_taxo$Kingdom=="Fungi",] #surement different FU PR
  }
  if(group=="PR"){
    print(paste0("n_unassigned(removed): ",sum(OTUdata_taxo$Phylum!="Not_Protist"))) #surement different FU PR #8186
    OTUdata <- OTUdata[,OTUdata_taxo$Phylum!="Not_Protist"] #surement different FU PR
    OTUdata_taxo<- OTUdata_taxo[OTUdata_taxo$Phylum!="Not_Protist",] #surement different FU PR
  }
  print(paste0("no problem?: ",all(OTUdata_taxo$Seq==colnames(OTUdata)) & all(colnames(OTUdata)==rownames(OTUdata_taxo))))
  
  print(paste0("n_assigned(kept): ",ncol(OTUdata))) # 52381
  
  print(paste0("prevalence filter, enough ndata to model?"))
  numberpresence <- apply(OTUdata,2,function(X){sum(X!=0)})
  plot(density(numberpresence),main="OTU that can be modelled (enough/not too much presences")
  # abline(v=nrow(OTUdata_BA)/10,col="blue")#10% cut
  abline(v=nrow(OTUdata)/20,col="red")#5% cut
  # abline(v=nrow(OTUdata_BA)-(nrow(OTUdata_BA)/10),col="blue")#10% cut
  abline(v=nrow(OTUdata)-(nrow(OTUdata)/20),col="red") #5% cut
  densityx<-lapply(density(numberpresence)$x,function(X){ifelse(X>nrow(OTUdata)/20,nrow(OTUdata)/20,X)})
  polygon(densityx,density(numberpresence)$y,col="grey")
  densityx2<-lapply(density(numberpresence)$x,function(X){ifelse(X<nrow(OTUdata)-(nrow(OTUdata)/20),nrow(OTUdata)-(nrow(OTUdata)/20),X)})
  polygon(densityx2,density(numberpresence)$y,col="grey")
  
  removed_toorare<-OTUdata[,numberpresence<nrow(OTUdata)/20]
  removed_toorare_taxo<-OTUdata_taxo[numberpresence<nrow(OTUdata)/20,]
  removed_toogeneral<-OTUdata[,numberpresence>nrow(OTUdata)-nrow(OTUdata)/20]
  removed_toogeneral_taxo<-OTUdata_taxo[numberpresence>nrow(OTUdata)-nrow(OTUdata)/20,]
  
  # rownames(removed_toorare)
  # removed_toomuchrare <- sum(numberpresence<nrow(OTUdata)/20)
  # removed_toogeneral <- sum(numberpresence>nrow(OTUdata)-nrow(OTUdata)/20)

  # removed_toomuchrare #2404
  # removed_toogeneral #796  #they are only removed in PA models
  # removed_toomuchrare+removed_toogeneral #2857
  if(group=="BA"){
    print("removed_rare_Bacteria")
    print(sum(removed_toorare_taxo$Kingdom=="Bacteria",na.rm = TRUE))#2368
    print("removed_rare_Archaea")
    print(sum(removed_toorare_taxo$Kingdom=="Archaea",na.rm = TRUE))#36
    print("removed_ubiquist_Bacteria")
    print(sum(removed_toogeneral_taxo$Kingdom=="Bacteria",na.rm = TRUE))#796
    print("removed_ubiquist_Archaea")
    print(sum(removed_toogeneral_taxo$Kingdom=="Archaea",na.rm = TRUE))#0
  }
  if(group=="FU"){
    print("removed_rare_fungi")
    print(sum(removed_toorare_taxo$Kingdom=="Fungi",na.rm = TRUE))#2368
    print("removed_ubiquist_fungi")
    print(sum(removed_toogeneral_taxo$Kingdom=="Fungi",na.rm = TRUE))#796
  }
  if(group=="PR"){
    print("removed_rare_protist")
    print(sum(removed_toorare_taxo$Phylum!="Not_Protist",na.rm = TRUE))#2368
    print("removed_ubiquist_protist")
    print(sum(removed_toogeneral_taxo$Phylum!="Not_Protist",na.rm = TRUE))#796
  }

  OTUdata_AB <- OTUdata[,numberpresence>=nrow(OTUdata)/20]
  OTUdata <- OTUdata[,which((numberpresence>=nrow(OTUdata)/20)&(numberpresence<(nrow(OTUdata)-nrow(OTUdata)/20)))]
  # OTUdata[,5600:5610]
  OTUdata <- apply(OTUdata,2,FUN=function(X){ifelse(X>0,1,0)}) 

  # test2<-apply(OTUdata,2,FUN=function(X){ifelse(X>0,1,0)}) 
  #  test2[,5600:5610]

  print("nOTU postfilter: ")
  print(ncol(OTUdata))
  

  saveRDS(OTUdata,file=paste0("data/OTUdata_",group,".Rds"))
  print("check nsite AB ")
  print(length(TotSeqSum)==nrow(OTUdata_AB))
  saveRDS(OTUdata_AB,file=paste0("data/OTUdata_",group,"_AB.Rds"))
  saveRDS(TotSeqSum,file=paste0("data/TotSeqSum_",group,"_AB.Rds"))
  
}




#old stuff (article 1 version)
# ################################################################################
# readcounts_BA<-readRDS(file="data\\readcounts_BA_Guex.Rds")
# readcounts_FU<-readRDS(file="data\\readcounts_FU_Guex.Rds")
# readcounts_PR<-readRDS(file="data\\readcounts_PR_Guex.Rds")
# 
# ENVdata<-readRDS(file="data\\ENVdata_BA.Rds")
# saveRDS(ENVdata_BA,file="data\\ENVdata_BA.Rds")
# ENVdata_BA<-ENVdata
# load(file="data\\ENVdata_FU.Rds")
# saveRDS(ENVdata_BA,file="data\\ENVdata_BA.Rds")
# ENVdata_FU<-ENVdata
# load(file="data\\ENVdata_PR.Rds")
# ENVdata_PR<-ENVdata
# # load("../spatial_data/MpAlps_soil_data_VVerdon1903+.Rda")
# 
# 
# ### Only to have an idea about taxonomy in data filtering steps
# load("../../30_data/ASVtables_taxo/BAtaxo2023_full.Rda")
# # saveRDS(BAtaxo2019,paste0("../../30_data/ASVtables_taxo/BAtaxo2023_full.Rds"))
# load("../../30_data/ASVtables_taxo/PRtaxo2023_full.Rda")
# # saveRDS(PRtaxo2023,paste0("../../30_data/ASVtables_taxo/PRtaxo2023_full.Rds"))
# load("../../30_data/ASVtables_taxo/FUtaxo2023_full.Rda")
# # saveRDS(FUtaxo2021,paste0("../../30_data/ASVtables_taxo/FUtaxo2023_full.Rds"))
# 
# #ASV data preparation
# #Bacteria
# OTUdata_BA <- readcounts_BA[,which(colnames(readcounts_BA)%in%rownames(ENVdata_BA))]
# 
# # nrow(OTUdata_BA) #60567
# OTUdata_BA <- OTUdata_BA[rowSums(OTUdata_BA)>99,]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them
# 
# #corresponding taxo
# OTUdata_BA_taxo<-cbind(OTUdata_BA ,BAtaxo2019[match(rownames(OTUdata_BA),BAtaxo2019[,"Seq"]),])
# all(OTUdata_BA_taxo$Seq==rownames(OTUdata_BA)) #All good taxo?
# 
# 
# #############################
# OTUdata_BA <- t(OTUdata_BA)
# # nrow(OTUdata_BA) #250
# # ncol(OTUdata_BA) #60567 ASV
# # sum(OTUdata_BA_taxo$Kingdom=="Bacteria",na.rm = TRUE)# 50353 Bacteria
# # sum(OTUdata_BA_taxo$Kingdom=="Archaea",na.rm = TRUE)#  187 Archaea
# # sum(OTUdata_BA_taxo$Kingdom=="unclassified_Root",na.rm = TRUE)# 7922
# 
# 
# 
# #easier datset for future (remove unassigned taxa directly instead of post-modelling):
# OTUdata_BA <- OTUdata_BA[,OTUdata_BA_taxo$Kingdom!="unclassified_Root"]
# OTUdata_BA_taxo2<-BAtaxo2019[match(colnames(OTUdata_BA),BAtaxo2019[,"Seq"]),]
# # sum(OTUdata_BA_taxo2$Kingdom=="unclassified_Root",na.rm = TRUE)# 7922
# #ncol(OTUdata_BA) #50540
# 
# #Need treshold to remove very high and very low prevalence (algos wont be able to adjust to data)
# #5 or 10% ?
# numberpresence_BA <- apply(OTUdata_BA,2,function(X){sum(X!=0)})
# # plot(numberpresence_BA,pch=20,cex=.5, main="number of plots in which the ASV \n is present")
# # # abline(h=nrow(OTUdata_BA_AB)/10,col="blue")#10% cut
# # abline(h=nrow(OTUdata_BA)/20,col="red") #5% cut
# # # abline(h=nrow(OTUdata_BA_AB)-(nrow(OTUdata_BA_AB)/10),col="blue")#10% cut
# # abline(h=nrow(OTUdata_BA)-(nrow(OTUdata_BA)/20),col="red") #5% cut
# plot(density(numberpresence_BA),main="Bacteria OTU that can be modelled (enough/not too much presences")
# # abline(v=nrow(OTUdata_BA)/10,col="blue")#10% cut
# abline(v=nrow(OTUdata_BA)/20,col="red")#5% cut
# # abline(v=nrow(OTUdata_BA)-(nrow(OTUdata_BA)/10),col="blue")#10% cut
# abline(v=nrow(OTUdata_BA)-(nrow(OTUdata_BA)/20),col="red") #5% cut
# densityx<-lapply(density(numberpresence_BA)$x,function(X){ifelse(X>nrow(OTUdata_BA)/20,nrow(OTUdata_BA)/20,X)})
# polygon(densityx,density(numberpresence_BA)$y,col="grey")
# densityx2<-lapply(density(numberpresence_BA)$x,function(X){ifelse(X<nrow(OTUdata_BA)-(nrow(OTUdata_BA)/20),nrow(OTUdata_BA)-(nrow(OTUdata_BA)/20),X)})
# polygon(densityx2,density(numberpresence_BA)$y,col="grey")
# 
# removed_toomuchrare <- sum(numberpresence_BA<nrow(OTUdata_BA)/20)
# removed_toogeneral <- sum(numberpresence_BA>nrow(OTUdata_BA)-nrow(OTUdata_BA)/20)
# removed_toorare_taxo<-OTUdata_BA_taxo2[numberpresence_BA<nrow(OTUdata_BA)/20,]
# removed_toogeneral_taxo<-OTUdata_BA_taxo2[numberpresence_BA>nrow(OTUdata_BA)-nrow(OTUdata_BA)/20,]
# 
# # removed_toomuchrare #2061
# # removed_toogeneral #796  #they are only removed in PA models
# # removed_toomuchrare+removed_toogeneral #2857
# sum(removed_toorare_taxo$Kingdom=="Bacteria",na.rm = TRUE)#2037
# sum(removed_toorare_taxo$Kingdom=="Archaea",na.rm = TRUE)#24
# sum(removed_toorare_taxo$Kingdom=="unclassified_Root",na.rm = TRUE)# 752    if 0 -> removal ok
# sum(removed_toogeneral_taxo$Kingdom=="Bacteria",na.rm = TRUE)#796
# sum(removed_toogeneral_taxo$Kingdom=="Archaea",na.rm = TRUE)#0
# sum(removed_toogeneral_taxo$Kingdom=="unclassified_Root",na.rm = TRUE)#21
# 
# #remove extreme prevalence (not enough data to model)
# OTUdata_BA_AB <- OTUdata_BA[,numberpresence_BA>=nrow(OTUdata_BA)/20]
# OTUdata_BA <- OTUdata_BA[,(numberpresence_BA>=nrow(OTUdata_BA)/20)*(numberpresence_BA<nrow(OTUdata_BA)-nrow(OTUdata_BA)/20)]
# 
# OTUdata_BA <- apply(OTUdata_BA,c(1,2),FUN=function(X){X>0}) 
# 
# # ncol(OTUdata_BA) #47683 ASV (817 removed if prevalence <5%)
# # ncol(OTUdata_BA_AB) #48479 ASV
# #253arrays of 220 ASV
# TotSeqSum_BA <- rowSums(OTUdata_BA_AB)
# 
# OTUdata <- OTUdata_BA
# # save(OTUdata,file="data/OTUdata_BA.Rds")
# 
# OTUdata <- OTUdata_BA_AB
# TotSeqSum <- TotSeqSum_BA
# # save(OTUdata,file="data/OTUdata_AB_BA.Rds")
# # save(TotSeqSum,file="data/TotSeqSum_BA.Rds")
# 
# 
# 
# #################################Fungi
# 
# OTUdata_FU_AB <- OTU_data_fun_232samples[,which(colnames(OTU_data_fun_232samples)%in%dataSoil$sampleNameFU)]
# OTUdata_FU_AB_taxo<-cbind(OTUdata_FU_AB ,FUtaxo2021[match(rownames(OTUdata_FU_AB),FUtaxo2021[,"Seq"]),])
# all(OTUdata_FU_AB_taxo$Seq==rownames(OTUdata_FU_AB_taxo))
# 
# OTUdata_FU_AB <- t(OTUdata_FU_AB)
# # nrow(OTUdata_FU_AB) #232
# #ncol(OTUdata_FU_AB) #95387
# # sum(OTUdata_FU_AB_taxo$Kingdom=="Fungi",na.rm = TRUE) 31112 Fungi
# # sum(OTUdata_FU_AB_taxo$Kingdom!="Fungi",na.rm = TRUE) 64261 not Fungi (or unclassified)
# # sum(is.na(OTUdata_FU_AB_taxo$Kingdom)) 0
# 
# OTUdata_FU_AB <- OTUdata_FU_AB[rownames(ENVdata_FU), ] # #remove plots that dont have measures for all variables
# # nrow(OTUdata_FU_AB) #217
# OTUdata_FU_AB <- OTUdata_FU_AB[,colSums(OTUdata_FU_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them
# 
# # ncol(OTUdata_FU_AB) #92753
# OTUdata_FU_AB_taxo2<-FUtaxo2021[match(colnames(OTUdata_FU_AB),FUtaxo2021[,"Seq"]),]
# all(OTUdata_FU_AB_taxo2$Seq==colnames(OTUdata_FU_AB))
# 
# # sum(OTUdata_FU_AB_taxo2$Kingdom=="Fungi",na.rm = TRUE)# 30366 Fungi
# # sum(OTUdata_FU_AB_taxo2$Kingdom!="Fungi",na.rm = TRUE)# 62387 not Fungi (or unclassified)
# # sum(is.na(OTUdata_FU_AB_taxo2$Kingdom))# 0
# 
# #easier datset for future (remove unassigned taxa directly instead of post-modelling): #ignored for publication 2023
# OTUdata_FU_AB <- OTUdata_FU_AB[,OTUdata_FU_AB_taxo2$Kingdom=="Fungi"]
# OTUdata_FU_AB_taxo2<-FUtaxo2021[match(colnames(OTUdata_FU_AB),FUtaxo2021[,"Seq"]),]
# #
# 
# numberpresence_FU <- apply(OTUdata_FU_AB,2,function(X){sum(X!=0)})
# # plot(numberpresence_FU,pch=20,cex=.5, main="number of plots in which the ASV \n is present (only ASV with >100 reads) ")
# # abline(h=nrow(OTUdata_FU_AB)/10,col="blue")#10% cut 21
# # abline(h=nrow(OTUdata_FU_AB)/20,col="red") #5% cut 10
# # # abline(h=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/10),col="blue")#10% cut
# # # abline(h=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/20),col="red") #5% cut
# # plot(density(numberpresence_FU))
# # abline(v=nrow(OTUdata_FU_AB)/10,col="blue")#10% cut
# # abline(v=nrow(OTUdata_FU_AB)/20,col="red")#5% cut
# # #abline(v=nrow(OTUdata_FU)-(nrow(OTUdata_FU)/10),col="blue")#10% cut
# # #abline(v=nrow(OTUdata_FU_AB)-(nrow(OTUdata_FU_AB)/20),col="red") #5% cut
# 
# 
# 
# removed_toorare <- sum(numberpresence_FU<nrow(OTUdata_FU_AB)/20)
# removed_toogeneral <- sum(numberpresence_FU>nrow(OTUdata_FU_AB)-nrow(OTUdata_FU_AB)/20)
# removed_toorare_taxo<-OTUdata_FU_AB_taxo2[numberpresence_FU<nrow(OTUdata_FU_AB)/20,]
# removed_toogeneral_taxo<-OTUdata_FU_AB_taxo2[numberpresence_FU>nrow(OTUdata_FU_AB)-nrow(OTUdata_FU_AB)/20,]
# 
# # removed_toorare #40890
# # removed_toogeneral #34  #they are only removed in PA models
# # removed_toorare+removed_toogeneral #40924
# sum(removed_toorare_taxo$Kingdom=="Fungi",na.rm = TRUE)#13021
# sum(removed_toorare_taxo$Kingdom!="Fungi",na.rm = TRUE)#27869
# sum(is.na(removed_toorare_taxo$Kingdom!="Fungi"))
# 
# sum(removed_toogeneral_taxo$Kingdom=="Fungi",na.rm = TRUE)#27
# sum(removed_toogeneral_taxo$Kingdom!="Fungi",na.rm = TRUE)#7
# 
# # OTUdata_FU_PA <- OTUdata_FU_PA[,numberpresence_FU_PA>=nrow(OTUdata_FU_PA)/20]
# 
# #remove less prevalent ones for AB
# OTUdata_FU_AB <- OTUdata_FU_AB[,numberpresence_FU>=nrow(OTUdata_FU_AB)/20]
# OTUdata_FU_PA <- OTUdata_FU_AB[,(numberpresence_FU>=nrow(OTUdata_FU_AB)/20)*(numberpresence_FU<nrow(OTUdata_FU_AB)-nrow(OTUdata_FU_AB)/20)]
# 
# # ncol(OTUdata_FU_AB)# 51863
# # ncol(OTUdata_FU_PA)# 51829
# 
# OTUdata_FU_PA <- apply(OTUdata_FU_AB,c(1,2),FUN=function(X){X>0})
# 
# 
# seqvec <- colnames(OTUdata_FU_AB)
# names(seqvec) <- paste0("OTU",1:ncol(OTUdata_FU_AB))
# # save(seqvec,file="seqvecFU.Rda")
# 
# # ncol(OTUdata_FU_AB) #51863
# # ncol(OTUdata_FU_PA) #51863
# 
# TotSeqSum_FU <- rowSums(OTUdata_FU_AB)
# 
# 
# 
# 
# #Protists
# OTUdata_PR_AB <- OTU_data_pro_179samples[,which(colnames(OTU_data_pro_179samples)%in%dataSoil$sampleNamePR)]
# 
# 
# OTUdata_PR_AB_taxo<-cbind(OTUdata_PR_AB ,PRtaxo2023[match(rownames(OTUdata_PR_AB),PRtaxo2023[,"Seq"]),])
# all(OTUdata_PR_AB_taxo$Seq==rownames(OTUdata_PR_AB_taxo))
# OTUdata_PR_AB <- t(OTUdata_PR_AB)
# 
# # nrow(OTUdata_PR_AB) #174
# # ncol(OTUdata_PR_AB) #11239
# # sum(OTUdata_PR_AB_taxo$Phylum!="Not_Protist",na.rm = TRUE) 2331 Protists
# # sum(OTUdata_PR_AB_taxo$Phylum=="Not_Protist",na.rm = TRUE) 8908 not Fungi (or unclassified)
# # sum(is.na(OTUdata_PR_AB_taxo$Phylum)) 0
# 
# # TotSeqSum_PR_AB <- rowSums(OTUdata_PR) #sequencing depth for "precleaned" data
# OTUdata_PR_AB <- OTUdata_PR_AB[rownames(ENVdata_PR), ] # #remove plots that dont have measures for all variables
# 
# # nrow(OTUdata_PR_AB) #166
# OTUdata_PR_AB <- OTUdata_PR_AB[,colSums(OTUdata_PR_AB)>99 ]# the selection of plots made some OTU to go under the threshold of 100 total read counts, remove them
# # ncol(OTUdata_PR_AB) #10860
# OTUdata_PR_AB_taxo2<-PRtaxo2023[match(colnames(OTUdata_PR_AB),PRtaxo2023[,"Seq"]),]
# all(OTUdata_PR_AB_taxo2$Seq==colnames(OTUdata_PR_AB))
# 
# # sum(OTUdata_PR_AB_taxo2$Phylum!="Not_Protist",na.rm = TRUE) 2270 Protists
# # sum(OTUdata_PR_AB_taxo2$Phylum=="Not_Protist",na.rm = TRUE) 8590 not Fungi (or unclassified)
# # sum(is.na(OTUdata_PR_AB_taxo$Phylum)) 0
# 
# # GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT
# #Remove non protists
# # troph_euk<-read.csv("../../ASV_data/Troph_Prot.csv",sep=";")
# # troph_euk<-read.csv("Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/Tax_table_Protists.csv",sep=",")
# 
# # OTUdata_PR_AB <- OTUdata_PR_AB[,which(colnames(OTUdata_PR_AB)%in%troph_euk$X)]
# 
# 
# #easier datset for future (remove unassigned taxa directly instead of post-modelling): #ignored for publication 2023
# OTUdata_PR_AB <- OTUdata_PR_AB[,OTUdata_PR_AB_taxo2$Phylum!="Not_Protist"]
# OTUdata_PR_AB_taxo2<-PRtaxo2023[match(colnames(OTUdata_PR_AB),PRtaxo2023[,"Seq"]),]
# #
# # ncol(OTUdata_PR_AB ) #2270
# 
# numberpresence_PR <- apply(OTUdata_PR_AB,2,function(X){sum(X!=0)})
# # plot(numberpresence_PR,pch=20,cex=.5, main="number of plots in which the ASV \n is present (only ASV with >100 reads) ")
# # abline(h=nrow(OTUdata_PR_AB)/10,col="blue")#10% cut
# # abline(h=nrow(OTUdata_PR_AB)/20,col="red") #5% cut
# # # abline(h=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/10),col="blue")#10% cut
# # # abline(h=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20),col="red") #5% cut
# # plot(density(numberpresence_PR))
# # abline(v=nrow(OTUdata_PR_AB)/10,col="blue")#10% cut
# # abline(v=nrow(OTUdata_PR_AB)/20,col="red")#5% cut
# # # abline(v=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/10),col="blue")#10% cut
# # # abline(v=nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20),col="red") #5% cut
# 
# removed_toorare <- sum(numberpresence_PR<nrow(OTUdata_PR_AB)/20)
# removed_toogeneral <- sum(numberpresence_PR>nrow(OTUdata_PR_AB)-nrow(OTUdata_PR_AB)/20)
# removed_toorare_taxo<-OTUdata_PR_AB_taxo2[numberpresence_PR<nrow(OTUdata_PR_AB)/20,]
# removed_toogeneral_taxo<-OTUdata_PR_AB_taxo2[numberpresence_PR>nrow(OTUdata_PR_AB)-nrow(OTUdata_PR_AB)/20,]
# # removed_toorare 108
# # removed_toogeneral 8  #they are only removed in PA models
# # removed_toorare+removed_toogeneral #116
# sum(removed_toorare_taxo$Phylum!="Not_Protist",na.rm = TRUE)#108
# sum(removed_toorare_taxo$Phylum=="Not_Protist",na.rm = TRUE)#0
# sum(is.na(removed_toorare_taxo$Phylum))
# 
# sum(removed_toogeneral_taxo$Phylum!="Not_Protist",na.rm = TRUE)#27
# sum(removed_toogeneral_taxo$Phylum=="Not_Protist",na.rm = TRUE)#7
# 
# 
# #remove less prevalent ones for AB
# OTUdata_PR_AB <- OTUdata_PR_AB[,numberpresence_PR>=nrow(OTUdata_PR_AB)/20]
# OTUdata_PR_PA <- OTUdata_PR_AB[,(numberpresence_PR>=nrow(OTUdata_PR_AB)/20)*(numberpresence_PR<nrow(OTUdata_PR_AB)-nrow(OTUdata_PR_AB)/20)]
# 
# # ncol(OTUdata_PR_AB)# 2162
# # ncol(OTUdata_PR_PA)# 2154
# 
# OTUdata_PR_PA <- apply(OTUdata_PR_AB,c(1,2),FUN=function(X){X>0})
# # ncol(OTUdata_PR_AB) #2162
# # which(colnames(OTUdata_PR_AB)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
# #annoying  Champi = 3019
# 
# 
# 
# 
# # ncol(OTUdata_PR_PA) #3160
# 
# seqvec <- colnames(OTUdata_PR_AB)
# names(seqvec) <- paste0("OTU",1:ncol(OTUdata_PR_AB))
# # save(seqvec,file="seqvecPR.Rda")
# 
# #cut at 10% : BA = 26, FU = 23, PR = 17 
# #cut at 5% : BA = 13 FU = 11 PR = 8
# #Where to cut here? 
# #10% removes most of fungi (most of them are present in a very few plots)
# #5% leaves some protists with only 9 plots to fit the models
# TotSeqSum_PR <- rowSums(OTUdata_PR_AB)
# 
# # datapoints<-runif(100,10,40)
# # dporder<-datapoints[order(datapoints)]
# # 
# # data<-matrix(data=c(datapoints,))
# # plot(dporder,dnorm(dporder,25,3)*7.5,type="l",lwd=5,ylim=c(0,1))
# # abline(h=0)
# # abline(h=1)
# # points(dporder,rbinom(100,1,dnorm(dporder,25,3)*7.5),pch=20,col="blue",cex=2)
# 
# colnames(ENVdata_BA)
# #save all as "OTUdata"
# OTUdata <- OTUdata_BA_PA
# ENVdata <- ENVdata_BA
# # save(OTUdata,file="PA/BA/Outputs/OTUdata.Rda")
# 
# # save(OTUdata,file="PA/BA/data/OTUdata.Rda")
# # save(ENVdata,file="PA/BA/Outputs/ENVdata.Rda")
# 
# 
# OTUdata <- OTUdata_BA_AB
# TotSeqSum <- TotSeqSum_BA
# # save(OTUdata,file="AB/BA/Outputs/OTUdata.Rda")
# # save(ENVdata,file="AB/BA/Outputs/ENVdata.Rda")
# # save(TotSeqSum,file="AB/BA/Outputs/dataTotSeqSum.Rda")
# 
# 
# OTUdata <- OTUdata_FU_PA
# ENVdata <- ENVdata_FU
# # save(OTUdata,file="PA/FU/Outputs/OTUdata.Rda")
# # save(ENVdata,file="PA/FU/Outputs/ENVdata.Rda")
# # save(ENVdata,file="PA/FU/data/ENVdata.Rda")
# 
# OTUdata <- OTUdata_FU_AB
# TotSeqSum <- TotSeqSum_FU
# # save(OTUdata,file="AB//FU/Outputs/OTUdata.Rda")
# # save(ENVdata,file="AB//FU/Outputs/ENVdata.Rda")
# # save(TotSeqSum,file="AB//FU/Outputs/dataTotSeqSum.Rda")
# 
# OTUdata <- OTUdata_PR_PA
# ENVdata <- ENVdata_PR
# # save(OTUdata,file="PA/PR/Outputs/OTUdata.Rda")
# # save(ENVdata,file="PA/PR/Outputs/ENVdata.Rda")
# # save(ENVdata,file="PA/PR/data/ENVdata.Rda")
# OTUdata <- OTUdata_PR_AB
# TotSeqSum <- TotSeqSum_PR
# # save(OTUdata,file="AB//PR/Outputs/OTUdata.Rda")
# # save(ENVdata,file="AB//PR/Outputs/ENVdata.Rda")
# # save(TotSeqSum,file="AB//PR/Outputs/dataTotSeqSum.Rda")
# 
# # write.csv(ENVdata,file="ENVdata_FU")

