

PAAB="PA"
# PAAB="AB"
GtoM="PR"
algo="GBM"

load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo_grouped_phylum.Rda"))
Taxo<-eval(parse(text=paste0(GtoM,"taxo_grouped_phylum")))

# initialisations
Variable_importance$marker="16S_bacteria-archaea"
Variable_importance<-Variable_importance[,c(ncol(Variable_importance)-1,ncol(Variable_importance),1:(ncol(Variable_importance)-2))]
# nrow(Variable_importance)
#add taxo
# sum(Variable_importance$Phylum=="unknown")
Variable_importance<-cbind(Variable_importance[Variable_importance$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_importance$seq,"Phylum"])
colnames(Variable_importance)[ncol(Variable_importance)]<-"Phylum"
#remove not Bact-Arch
# Variable_importance<-Variable_importance[Variable_importance$Phylum!="unknown",]

# VarSel_PA_GLM <- cbind(Variable_importance[Variable_importance$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_importance$seq,"Phylum"])[FALSE,]

VarSel_PA_GLM <- VarSel_PA_RF <- VarSel_PA_GBM <- VarSel_PA_GAM <- VarSel_PA_GLM <- Variable_importance[FALSE,]



PAAB="AB"

load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo_grouped_phylum.Rda"))
Taxo<-eval(parse(text=paste0(GtoM,"taxo_grouped_phylum")))

# initialisations
Variable_importance$marker="16S_bacteria-archaea"
Variable_importance<-Variable_importance[,c(ncol(Variable_importance)-1,ncol(Variable_importance),1:(ncol(Variable_importance)-2))]
Variable_importance<-cbind(Variable_importance[Variable_importance$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_importance$seq,"Phylum"])
colnames(Variable_importance)[ncol(Variable_importance)]<-"Phylum"
VarSel_AB_RF <- VarSel_AB_GBM <- VarSel_AB_GAM <- VarSel_AB_GLM <- Variable_importance[FALSE,]

# filling
for(PAAB in c("PA","AB")){
  print(PAAB)
  for(algo in c("GLM","GAM","GBM","RF")){#Pb empty Variable_ranks.Rda GAM+GBM PA --> Bug Figure Selected Var
    print(algo)
    for(GtoM in c("BA","FU","PR")){
      print(GtoM)
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
      load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo_grouped_phylum.Rda"))
      Taxo<-eval(parse(text=paste0(GtoM,"taxo_grouped_phylum")))
      if(any(Variable_importance$seq!=Taxo$Seq)){ #will trigger warning but works
        Variable_importance<-Variable_importance[which(Variable_importance$seq%in%Taxo$Seq),]
        any(Variable_importance$seq!=Taxo$Seq)
      }
      if (GtoM=="PR"){
        Variable_importance <-data.frame(Variable_importance[Taxo$Phylum!="Not_Protist",],Phylum=Taxo$Phylum[Taxo$Phylum!="Not_Protist"])
        Variable_importance$marker="18S_protists"
        nrow(Variable_importance)
      }
      if (GtoM=="FU"){
        Variable_importance <-data.frame(Variable_importance[Taxo$Kingdom=="Fungi",],Phylum=Taxo$Phylum[Taxo$Kingdom=="Fungi"])
        Variable_importance$marker="ITS_fungi"
      }    
      if (GtoM=="BA"){
        Variable_importance <-data.frame(Variable_importance[Taxo$Kingdom!="unclassified_Root",],Phylum=Taxo$Phylum[Taxo$Kingdom!="unclassified_Root"])
        Variable_importance$marker="16S_bacteria-archaea"
      }
      
      Variable_importance <-Variable_importance[,c((ncol(Variable_importance)-2):ncol(Variable_importance),1:(ncol(Variable_importance)-3))]
      Variable_importance[,-(1:3)] <- t(apply(Variable_importance[,-(1:3)],1,function(X){X/max(X,na.rm=TRUE)})) #maxvalue = 1, then rescaled
      Variable_importance[,-(1:3)] <- t(apply(Variable_importance[,-(1:3)],1,function(X){ifelse(X==0,NA,X)})) #maxvalue = 1, then rescaled
      assign(paste0("VarSel_",PAAB,"_",algo),rbind(eval(parse(text=paste0("VarSel_",PAAB,"_",algo))),Variable_importance))
    }
  }
}

VarSel_PA_GLM[]<-lapply(VarSel_PA_GLM,function(X){ifelse(X==0,NA,X)})
VarSel_PA_GLM<-data.frame(VarSel_PA_GLM)
VarSel_PA_GLM<-VarSel_PA_GLM[c(1:2,ncol(VarSel_PA_GLM),3:(ncol(VarSel_PA_GLM)-1))]
VarSel_PA_GAM[]<-lapply(VarSel_PA_GAM,function(X){ifelse(X==0,NA,X)})
VarSel_PA_GAM<-data.frame(VarSel_PA_GAM)
VarSel_PA_GAM<-VarSel_PA_GAM[c(1:2,ncol(VarSel_PA_GAM),3:(ncol(VarSel_PA_GAM)-1))]
VarSel_PA_GBM[]<-lapply(VarSel_PA_GBM,function(X){ifelse(X==0,NA,X)})
VarSel_PA_GBM<-data.frame(VarSel_PA_GBM)
VarSel_PA_GBM<-VarSel_PA_GBM[c(1:2,ncol(VarSel_PA_GBM),3:(ncol(VarSel_PA_GBM)-1))]
VarSel_PA_RF[]<-lapply(VarSel_PA_RF,function(X){ifelse(X==0,NA,X)})
VarSel_PA_RF<-data.frame(VarSel_PA_RF)
VarSel_PA_RF<-VarSel_PA_RF[c(1:2,ncol(VarSel_PA_RF),3:(ncol(VarSel_PA_RF)-1))]
VarSel_AB_GLM[]<-lapply(VarSel_AB_GLM,function(X){ifelse(X==0,NA,X)})
VarSel_AB_GLM<-data.frame(VarSel_AB_GLM)
VarSel_AB_GLM<-VarSel_AB_GLM[c(1:2,ncol(VarSel_AB_GLM),3:(ncol(VarSel_AB_GLM)-1))]
VarSel_AB_GAM[]<-lapply(VarSel_AB_GAM,function(X){ifelse(X==0,NA,X)})
VarSel_AB_GAM<-data.frame(VarSel_AB_GAM)
VarSel_AB_GAM<-VarSel_AB_GAM[c(1:2,ncol(VarSel_AB_GAM),3:(ncol(VarSel_AB_GAM)-1))]
VarSel_AB_GBM[]<-lapply(VarSel_AB_GBM,function(X){ifelse(X==0,NA,X)})
VarSel_AB_GBM<-data.frame(VarSel_AB_GBM)
VarSel_AB_GBM<-VarSel_AB_GBM[c(1:2,ncol(VarSel_AB_GBM),3:(ncol(VarSel_AB_GBM)-1))]
VarSel_AB_RF[]<-lapply(VarSel_AB_RF,function(X){ifelse(X==0,NA,X)})
VarSel_AB_RF<-data.frame(VarSel_AB_RF)
VarSel_AB_RF<-VarSel_AB_RF[c(1:2,ncol(VarSel_AB_RF),3:(ncol(VarSel_AB_RF)-1))]

write.csv(VarSel_PA_GLM,file="figures/PAAB_selection/Figshare_tables/CovImp_GLM_PA.csv")
write.csv(VarSel_PA_GAM,file="figures/PAAB_selection/Figshare_tables/CovImp_GAM_PA.csv")
write.csv(VarSel_PA_GBM,file="figures/PAAB_selection/Figshare_tables/CovImp_GBM_PA.csv")
write.csv(VarSel_PA_RF,file="figures/PAAB_selection/Figshare_tables/CovImp_RF_PA.csv")
write.csv(VarSel_AB_GLM,file="figures/PAAB_selection/Figshare_tables/CovImp_GLM_AB.csv")
write.csv(VarSel_AB_GAM,file="figures/PAAB_selection/Figshare_tables/CovImp_GAM_AB.csv")
write.csv(VarSel_AB_GBM,file="figures/PAAB_selection/Figshare_tables/CovImp_GBM_AB.csv")
write.csv(VarSel_AB_RF,file="figures/PAAB_selection/Figshare_tables/CovImp_RF_AB.csv")

