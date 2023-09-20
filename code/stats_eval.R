library(tidyverse)
library(psych)

library(grid)
library(gridExtra)

dataset=c(AR='#5CB800',BA='#FFD700',FU='#6B4C62',PR='#457EB0')

###data gathering###
#PR

load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda"))
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.

#GLM
load("PA/PR/data/GLM/Eval.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGLM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(Eval)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGLM_PR2<-evalmetGLM_PR[1:18]
evalmetGLM_PR2$group<-"PR"
# summary(Eval)
#GAM
load("PA/PR/data/GAM/Eval.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGAM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(Eval)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGAM_PR2<-evalmetGAM_PR[1:18]
evalmetGAM_PR2$group<-"PR"

#GBM
load("PA/PR/data/GBM/Eval.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetGBM_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(Eval)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGBM_PR2<-evalmetGBM_PR[1:18]
evalmetGBM_PR2$group<-"PR"
#RF
load("PA/PR/data/RF/Eval.Rda")

if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
evalmetRF_PR<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(Eval)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetRF_PR2<-evalmetRF_PR[1:18]
evalmetRF_PR2$group<-"PR"

#FU
#GLM
load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.
load("PA/FU/data/GLM/Eval.Rda")
evalmetGLM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",-which(colnames(Eval)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGLM_FU2<-evalmetGLM_FU[1:18]
evalmetGLM_FU2$group<-"FU"
#GAM
load("PA/FU/data/GAM/Eval.Rda")
evalmetGAM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",-which(colnames(Eval)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGAM_FU2<-evalmetGAM_FU[1:18]
evalmetGAM_FU2$group<-"FU"
#GBM
load("PA/FU/data/GBM/Eval.Rda")
evalmetGBM_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",-which(colnames(Eval)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGBM_FU2<-evalmetGBM_FU[1:18]
evalmetGBM_FU2$group<-"FU"
#RF
load("PA/FU/data/RF/Eval.Rda")
evalmetRF_FU<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",-which(colnames(Eval)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetRF_FU2<-evalmetRF_FU[1:18]
evalmetRF_FU2$group<-"FU"

#BA
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum
#GLM
load("PA/BA/data/GLM/Eval.Rda")
evalmetGLM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGLM_BA2<-evalmetGLM_BA[1:18]
evalmetGLM_BA2$group<-"BA"
evalmetGLM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGLM_AR2<-evalmetGLM_AR[1:18]
evalmetGLM_AR2$group<-"AR"
#GAM
load("PA/BA/data/GAM/Eval.Rda")
evalmetGAM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGAM_BA2<-evalmetGAM_BA[1:18]
evalmetGAM_BA2$group<-"BA"
evalmetGAM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGAM_AR2<-evalmetGAM_AR[1:18]
evalmetGAM_AR2$group<-"AR"
#GBM
load("PA/BA/data/GBM/Eval.Rda")
evalmetGBM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGBM_BA2<-evalmetGBM_BA[1:18]
evalmetGBM_BA2$group<-"BA"
evalmetGBM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGBM_AR2<-evalmetGBM_AR[1:18]
evalmetGBM_AR2$group<-"AR"
#RF
load("PA/BA/data/RF/Eval.Rda")
evalmetRF_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetRF_BA2<-evalmetRF_BA[1:18]
evalmetRF_BA2$group<-"BA"
evalmetRF_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",-which(colnames(Eval)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetRF_AR2<-evalmetRF_AR[1:18]
evalmetRF_AR2$group<-"AR"

evalmetGLM_all<-rbind(evalmetGLM_PR2,evalmetGLM_FU2,evalmetGLM_BA2,evalmetGLM_AR2)
obj_to_plot_GLM<-evalmetGLM_all[!(is.na(evalmetGLM_all$TSS_adj)),c("TSS_adj","auc","kap","TSS","group")]
evalmetGAM_all<-rbind(evalmetGAM_PR2,evalmetGAM_FU2,evalmetGAM_BA2,evalmetGAM_AR2)
obj_to_plot_GAM<-evalmetGAM_all[!(is.na(evalmetGAM_all$TSS_adj)),c("TSS_adj","auc","kap","TSS","group")]
evalmetGBM_all<-rbind(evalmetGBM_PR2,evalmetGBM_FU2,evalmetGBM_BA2,evalmetGBM_AR2)
obj_to_plot_GBM<-evalmetGBM_all[!(is.na(evalmetGBM_all$TSS_adj)),c("TSS_adj","auc","kap","TSS","group")]
evalmetRF_all<-rbind(evalmetRF_PR2,evalmetRF_FU2,evalmetRF_BA2,evalmetRF_AR2)
obj_to_plot_RF<-evalmetRF_all[!(is.na(evalmetRF_all$TSS_adj)),c("TSS_adj","auc","kap","TSS","group")]
sum(evalmetGLM_all$TSS_sign,na.rm=TRUE)/length(evalmetGLM_all$TSS_adj)
sum(evalmetGLM_all$TSS_adj>0.4,na.rm=TRUE)/length(evalmetGLM_all$TSS_adj)
Summary_TSSadjGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_TSSadjGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_TSSadjGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_TSSadjRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS_adj, na.rm = TRUE),
    sd = sd(TSS_adj, na.rm = TRUE)
  )
Summary_aucGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(auc, na.rm = TRUE),
    sd = sd(auc, na.rm = TRUE)
  )
Summary_aucGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(auc, na.rm = TRUE),
    sd = sd(auc, na.rm = TRUE)
  )
Summary_aucGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(auc, na.rm = TRUE),
    sd = sd(auc, na.rm = TRUE)
  )
Summary_aucRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(auc, na.rm = TRUE),
    sd = sd(auc, na.rm = TRUE)
  )
Summary_kapGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(kap, na.rm = TRUE),
    sd = sd(kap, na.rm = TRUE)
  )
Summary_kapGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(kap, na.rm = TRUE),
    sd = sd(kap, na.rm = TRUE)
  )
Summary_kapGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(kap, na.rm = TRUE),
    sd = sd(kap, na.rm = TRUE)
  )
Summary_kapRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(kap, na.rm = TRUE),
    sd = sd(kap, na.rm = TRUE)
  )
Summary_TSSGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS, na.rm = TRUE),
    sd = sd(TSS, na.rm = TRUE)
  )
Summary_TSSGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS, na.rm = TRUE),
    sd = sd(TSS, na.rm = TRUE)
  )
Summary_TSSGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS, na.rm = TRUE),
    sd = sd(TSS, na.rm = TRUE)
  )
Summary_TSSRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(TSS, na.rm = TRUE),
    sd = sd(TSS, na.rm = TRUE)
  )

#taxo supplementary table
PRtaxo_supp<-PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
colnames(PRtaxo_supp)<-c("seq","level1","conf.1","level2","conf.2","level3","conf.3","level4","conf.4","level5","conf.5","level6","conf.6","level7","conf.7","level8","conf.8","level9","conf.9","Taxa")
# write.csv(PRtaxo_supp, file="figures/PAAB_selection/Figshare_tables/PRtaxo_supp.csv")
Taxo_BAARFU<-rbind(BAtaxo[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Archaea",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
unique(Taxo_BAARFU$Kingdom)
# write.csv(Taxo_BAARFU, file="figures/PAAB_selection/Figshare_tables/BAARFUtaxo_supp.csv")

#supp able model quality PA part
tobindPRGLM<-evalmetGLM_PR[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindPRGLM$marker<-"18S_protist"
tobindPRGLM<-tobindPRGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]
tobindFUGLM<-evalmetGLM_FU[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindFUGLM$marker<-"ITS_fungi"
tobindFUGLM<-tobindFUGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]
tobindBAGLM<-evalmetGLM_BA[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindBAGLM$marker<-"16S_bacteria-archaea"
tobindBAGLM<-tobindBAGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]
tobindARGLM<-evalmetGLM_AR[c("Seq","Phylum","auc","kap","TSS","TSS_adj")]
tobindARGLM$marker<-"16S_bacteria-archaea"
tobindARGLM<-tobindARGLM[c("Seq","marker","Phylum","auc","kap","TSS","TSS_adj")]

tobindPRGAM<-evalmetGAM_PR[c("auc","kap","TSS","TSS_adj")]
tobindFUGAM<-evalmetGAM_FU[c("auc","kap","TSS","TSS_adj")]
tobindBAGAM<-evalmetGAM_BA[c("auc","kap","TSS","TSS_adj")]
tobindARGAM<-evalmetGAM_AR[c("auc","kap","TSS","TSS_adj")]
tobindPRGBM<-evalmetGBM_PR[c("auc","kap","TSS","TSS_adj")]
tobindFUGBM<-evalmetGBM_FU[c("auc","kap","TSS","TSS_adj")]
tobindBAGBM<-evalmetGBM_BA[c("auc","kap","TSS","TSS_adj")]
tobindARGBM<-evalmetGBM_AR[c("auc","kap","TSS","TSS_adj")]
tobindPRRF<-evalmetRF_PR[c("auc","kap","TSS","TSS_adj")]
tobindFURF<-evalmetRF_FU[c("auc","kap","TSS","TSS_adj")]
tobindBARF<-evalmetRF_BA[c("auc","kap","TSS","TSS_adj")]
tobindARRF<-evalmetRF_AR[c("auc","kap","TSS","TSS_adj")]


evalmet_allPR<-cbind(tobindPRGLM,tobindPRGAM,tobindPRGBM,tobindPRRF)
colnames(evalmet_allPR)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allFU<-cbind(tobindFUGLM,tobindFUGAM,tobindFUGBM,tobindFURF)
colnames(evalmet_allFU)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allBA<-cbind(tobindBAGLM,tobindBAGAM,tobindBAGBM,tobindBARF)
colnames(evalmet_allBA)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allAR<-cbind(tobindARGLM,tobindARGAM,tobindARGBM,tobindARRF)
colnames(evalmet_allAR)<-c("Seq","marker","Phylum","GLM-auc","GLM-kap","GLM-TSS","GLM-TSS_adj","GAM-auc","GAM-kap","GAM-TSS","GAM-TSS_adj","GBM-auc","GBM-kap","GBM-TSS","GBM-TSS_adj","RF-auc","RF-kap","RF-TSS","RF-TSS_adj")
evalmet_allall<-rbind(evalmet_allPR,evalmet_allFU,evalmet_allBA,evalmet_allAR)
# write.csv(evalmet_allall, file="figures/PAAB_selection/Figshare_tables/Model_quality_phylotype_PA.csv")             

###Statstics###
anovaGLM<-aov(TSS_adj ~ group, data = evalmetGLM_all)
# summary(anovaGLM)
tuk_GLM<-TukeyHSD(anovaGLM,ordered=TRUE)$group
Cohens_D_GLM<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM)<-c("lower","effect","upper")
rownames(Cohens_D_GLM)<-rownames(tuk_GLM)
for (i in 1:nrow(tuk_GLM)){#i=1
  coupletolook <- word(rownames(tuk_GLM)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjGLM$mean[which(Summary_TSSadjGLM$group==coupletolook[1])] - Summary_TSSadjGLM$mean[which(Summary_TSSadjGLM$group==coupletolook[2])])/sd(evalmetGLM_all$TSS_adj, na.rm = TRUE)
  Cohens_D_GLM[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjGLM$count[which(Summary_TSSadjGLM$group==coupletolook[1])], n2=Summary_TSSadjGLM$count[which(Summary_TSSadjGLM$group==coupletolook[2])]) 
}
stats_res_GLM<-c()
stats_res_GLM$tukey<-c("a","b","c","d")
stats_res_GLM$cohen<-c("a","a","b","c")
stats_res_GLM<-as.data.frame(stats_res_GLM)

anovaGAM<-aov(TSS_adj ~ group, data = evalmetGAM_all)
# summary(anovaGAM)
tuk_GAM<-TukeyHSD(anovaGAM,ordered=TRUE)$group
Cohens_D_GAM<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM)<-c("lower","effect","upper")
rownames(Cohens_D_GAM)<-rownames(tuk_GAM)
for (i in 1:nrow(tuk_GAM)){#i=1
  coupletolook <- word(rownames(tuk_GAM)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjGAM$mean[which(Summary_TSSadjGAM$group==coupletolook[1])] - Summary_TSSadjGAM$mean[which(Summary_TSSadjGAM$group==coupletolook[2])])/sd(evalmetGAM_all$TSS_adj, na.rm = TRUE)
  Cohens_D_GAM[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjGAM$count[which(Summary_TSSadjGAM$group==coupletolook[1])], n2=Summary_TSSadjGAM$count[which(Summary_TSSadjGAM$group==coupletolook[2])]) 
}
stats_res_GAM<-c()
stats_res_GAM$tukey<-c("a","b","c","d")
stats_res_GAM$cohen<-c("a","a","b","c")
stats_res_GAM<-as.data.frame(stats_res_GAM)

anovaGBM<-aov(TSS_adj ~ group, data = evalmetGBM_all)
# summary(anovaGBM)
tuk_GBM<-TukeyHSD(anovaGBM,ordered=TRUE)$group
Cohens_D_GBM<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM)<-c("lower","effect","upper")
rownames(Cohens_D_GBM)<-rownames(tuk_GBM)
for (i in 1:nrow(tuk_GBM)){#i=1
  coupletolook <- word(rownames(tuk_GBM)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjGBM$mean[which(Summary_TSSadjGBM$group==coupletolook[1])] - Summary_TSSadjGBM$mean[which(Summary_TSSadjGBM$group==coupletolook[2])])/sd(evalmetGBM_all$TSS_adj, na.rm = TRUE)
  Cohens_D_GBM[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjGBM$count[which(Summary_TSSadjGBM$group==coupletolook[1])], n2=Summary_TSSadjGBM$count[which(Summary_TSSadjGBM$group==coupletolook[2])]) 
}
stats_res_GBM<-c()
stats_res_GBM$tukey<-c("a","b","c","d")
stats_res_GBM$cohen<-c("a","b","c","d")
stats_res_GBM<-as.data.frame(stats_res_GBM)

anovaRF<-aov(TSS_adj ~ group, data = evalmetRF_all)
# summary(anovaRF)
tuk_RF<-TukeyHSD(anovaRF,ordered=TRUE)$group
Cohens_D_RF<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF)<-c("lower","effect","upper")
rownames(Cohens_D_RF)<-rownames(tuk_RF)
for (i in 1:nrow(tuk_RF)){#i=1
  coupletolook <- word(rownames(tuk_RF)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSadjRF$mean[which(Summary_TSSadjRF$group==coupletolook[1])] - Summary_TSSadjRF$mean[which(Summary_TSSadjRF$group==coupletolook[2])])/sd(evalmetRF_all$TSS_adj, na.rm = TRUE)
  Cohens_D_RF[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSadjRF$count[which(Summary_TSSadjRF$group==coupletolook[1])], n2=Summary_TSSadjRF$count[which(Summary_TSSadjRF$group==coupletolook[2])]) 
}
stats_res_RF<-c()
stats_res_RF$tukey<-c("a","b","c","d")
stats_res_RF$cohen<-c("a","a","b","c")
stats_res_RF<-as.data.frame(stats_res_RF)


anovaGLM_auc<-aov(auc ~ group, data = evalmetGLM_all)
# summary(anovaGLM)
tuk_GLM_auc<-TukeyHSD(anovaGLM_auc,ordered=TRUE)$group # a b c d 
Cohens_D_GLM_auc<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM_auc)<-c("lower","effect","upper")
rownames(Cohens_D_GLM_auc)<-rownames(tuk_GLM_auc)
for (i in 1:nrow(tuk_GLM_auc)){#i=1
  coupletolook <- word(rownames(tuk_GLM_auc)[i],1:2,sep = "-")
  cohenD<-(Summary_aucGLM$mean[which(Summary_aucGLM$group==coupletolook[1])] - Summary_aucGLM$mean[which(Summary_aucGLM$group==coupletolook[2])])/sd(evalmetGLM_all$auc, na.rm = TRUE)
  Cohens_D_GLM_auc[i,]<-cohen.d.ci(d=cohenD,n1=Summary_aucGLM$count[which(Summary_aucGLM$group==coupletolook[1])], n2=Summary_aucGLM$count[which(Summary_aucGLM$group==coupletolook[2])]) 
}
stats_res_GLM_auc<-c()
stats_res_GLM_auc$tukey<-c("a","b","c","d")
stats_res_GLM_auc$cohen<-c("a","b","b","c")
stats_res_GLM_auc<-as.data.frame(stats_res_GLM_auc)

anovaGAM_auc<-aov(auc ~ group, data = evalmetGAM_all)
# summary(anovaGAM)
tuk_GAM_auc<-TukeyHSD(anovaGAM_auc,ordered=TRUE)$group # a b c d 
Cohens_D_GAM_auc<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM_auc)<-c("lower","effect","upper")
rownames(Cohens_D_GAM_auc)<-rownames(tuk_GAM_auc)
for (i in 1:nrow(tuk_GAM_auc)){#i=1
  coupletolook <- word(rownames(tuk_GAM_auc)[i],1:2,sep = "-")
  cohenD<-(Summary_aucGAM$mean[which(Summary_aucGAM$group==coupletolook[1])] - Summary_aucGAM$mean[which(Summary_aucGAM$group==coupletolook[2])])/sd(evalmetGAM_all$auc, na.rm = TRUE)
  Cohens_D_GAM_auc[i,]<-cohen.d.ci(d=cohenD,n1=Summary_aucGAM$count[which(Summary_aucGAM$group==coupletolook[1])], n2=Summary_aucGAM$count[which(Summary_aucGAM$group==coupletolook[2])]) 
}
stats_res_GAM_auc<-c()
stats_res_GAM_auc$tukey<-c("a","b","c","d")
stats_res_GAM_auc$cohen<-c("a","a","b","c")
stats_res_GAM_auc<-as.data.frame(stats_res_GAM_auc)

anovaGBM_auc<-aov(auc ~ group, data = evalmetGBM_all)
# summary(anovaGBM)
tuk_GBM_auc<-TukeyHSD(anovaGBM_auc,ordered=TRUE)$group # a b c d 
Cohens_D_GBM_auc<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM_auc)<-c("lower","effect","upper")
rownames(Cohens_D_GBM_auc)<-rownames(tuk_GBM_auc)
for (i in 1:nrow(tuk_GBM_auc)){#i=1
  coupletolook <- word(rownames(tuk_GBM_auc)[i],1:2,sep = "-")
  cohenD<-(Summary_aucGBM$mean[which(Summary_aucGBM$group==coupletolook[1])] - Summary_aucGBM$mean[which(Summary_aucGBM$group==coupletolook[2])])/sd(evalmetGBM_all$auc, na.rm = TRUE)
  Cohens_D_GBM_auc[i,]<-cohen.d.ci(d=cohenD,n1=Summary_aucGBM$count[which(Summary_aucGBM$group==coupletolook[1])], n2=Summary_aucGBM$count[which(Summary_aucGBM$group==coupletolook[2])]) 
}
stats_res_GBM_auc<-c()
stats_res_GBM_auc$tukey<-c("a","b","c","d")
stats_res_GBM_auc$cohen<-c("a","b","b","c")
stats_res_GBM_auc<-as.data.frame(stats_res_GBM_auc)

anovaRF_auc<-aov(auc ~ group, data = evalmetRF_all)
# summary(anovaRF)
tuk_RF_auc<-TukeyHSD(anovaRF_auc,ordered=TRUE)$group # a b c d 
Cohens_D_RF_auc<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF_auc)<-c("lower","effect","upper")
rownames(Cohens_D_RF_auc)<-rownames(tuk_RF_auc)
for (i in 1:nrow(tuk_RF_auc)){#i=1
  coupletolook <- word(rownames(tuk_RF_auc)[i],1:2,sep = "-")
  cohenD<-(Summary_aucRF$mean[which(Summary_aucRF$group==coupletolook[1])] - Summary_aucRF$mean[which(Summary_aucRF$group==coupletolook[2])])/sd(evalmetRF_all$auc, na.rm = TRUE)
  Cohens_D_RF_auc[i,]<-cohen.d.ci(d=cohenD,n1=Summary_aucRF$count[which(Summary_aucRF$group==coupletolook[1])], n2=Summary_aucRF$count[which(Summary_aucRF$group==coupletolook[2])]) 
}
stats_res_RF_auc<-c()
stats_res_RF_auc$tukey<-c("a","b","c","d")
stats_res_RF_auc$cohen<-c("a","a","b","c")
stats_res_RF_auc<-as.data.frame(stats_res_RF_auc)

anovaGLM_kap<-aov(kap ~ group, data = evalmetGLM_all)
# summary(anovaGLM)
tuk_GLM_kap<-TukeyHSD(anovaGLM_kap,ordered=TRUE)$group # a b c d 
Cohens_D_GLM_kap<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM_kap)<-c("lower","effect","upper")
rownames(Cohens_D_GLM_kap)<-rownames(tuk_GLM_kap)
for (i in 1:nrow(tuk_GLM_kap)){#i=1
  coupletolook <- word(rownames(tuk_GLM_kap)[i],1:2,sep = "-")
  cohenD<-(Summary_kapGLM$mean[which(Summary_kapGLM$group==coupletolook[1])] - Summary_kapGLM$mean[which(Summary_kapGLM$group==coupletolook[2])])/sd(evalmetGLM_all$kap, na.rm = TRUE)
  Cohens_D_GLM_kap[i,]<-cohen.d.ci(d=cohenD,n1=Summary_kapGLM$count[which(Summary_kapGLM$group==coupletolook[1])], n2=Summary_kapGLM$count[which(Summary_kapGLM$group==coupletolook[2])]) 
}
stats_res_GLM_kap<-c()
stats_res_GLM_kap$tukey<-c("a","a","b","c")
stats_res_GLM_kap$cohen<-c("a","a","b","c")
stats_res_GLM_kap<-as.data.frame(stats_res_GLM_kap)

anovaGAM_kap<-aov(kap ~ group, data = evalmetGAM_all)
# summary(anovaGAM)
tuk_GAM_kap<-TukeyHSD(anovaGAM_kap,ordered=TRUE)$group # a b c d 
Cohens_D_GAM_kap<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM_kap)<-c("lower","effect","upper")
rownames(Cohens_D_GAM_kap)<-rownames(tuk_GAM_kap)
for (i in 1:nrow(tuk_GAM_kap)){#i=1
  coupletolook <- word(rownames(tuk_GAM_kap)[i],1:2,sep = "-")
  cohenD<-(Summary_kapGAM$mean[which(Summary_kapGAM$group==coupletolook[1])] - Summary_kapGAM$mean[which(Summary_kapGAM$group==coupletolook[2])])/sd(evalmetGAM_all$kap, na.rm = TRUE)
  Cohens_D_GAM_kap[i,]<-cohen.d.ci(d=cohenD,n1=Summary_kapGAM$count[which(Summary_kapGAM$group==coupletolook[1])], n2=Summary_kapGAM$count[which(Summary_kapGAM$group==coupletolook[2])]) 
}
stats_res_GAM_kap<-c()
stats_res_GAM_kap$tukey<-c("a","a","b","c")
stats_res_GAM_kap$cohen<-c("a","a","b","c")
stats_res_GAM_kap<-as.data.frame(stats_res_GAM_kap)

anovaGBM_kap<-aov(kap ~ group, data = evalmetGBM_all)
# summary(anovaGBM)
tuk_GBM_kap<-TukeyHSD(anovaGBM_kap,ordered=TRUE)$group # a b c d 
Cohens_D_GBM_kap<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM_kap)<-c("lower","effect","upper")
rownames(Cohens_D_GBM_kap)<-rownames(tuk_GBM_kap)
for (i in 1:nrow(tuk_GBM_kap)){#i=1
  coupletolook <- word(rownames(tuk_GBM_kap)[i],1:2,sep = "-")
  cohenD<-(Summary_kapGBM$mean[which(Summary_kapGBM$group==coupletolook[1])] - Summary_kapGBM$mean[which(Summary_kapGBM$group==coupletolook[2])])/sd(evalmetGBM_all$kap, na.rm = TRUE)
  Cohens_D_GBM_kap[i,]<-cohen.d.ci(d=cohenD,n1=Summary_kapGBM$count[which(Summary_kapGBM$group==coupletolook[1])], n2=Summary_kapGBM$count[which(Summary_kapGBM$group==coupletolook[2])]) 
}
stats_res_GBM_kap<-c()
stats_res_GBM_kap$tukey<-c("a","a","b","c")
stats_res_GBM_kap$cohen<-c("a","a","b","c")
stats_res_GBM_kap<-as.data.frame(stats_res_GBM_kap)

anovaRF_kap<-aov(kap ~ group, data = evalmetRF_all)
# summary(anovaRF)
tuk_RF_kap<-TukeyHSD(anovaRF_kap,ordered=TRUE)$group # a b c d 
Cohens_D_RF_kap<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF_kap)<-c("lower","effect","upper")
rownames(Cohens_D_RF_kap)<-rownames(tuk_RF_kap)
for (i in 1:nrow(tuk_RF_kap)){#i=1
  coupletolook <- word(rownames(tuk_RF_kap)[i],1:2,sep = "-")
  cohenD<-(Summary_kapRF$mean[which(Summary_kapRF$group==coupletolook[1])] - Summary_kapRF$mean[which(Summary_kapRF$group==coupletolook[2])])/sd(evalmetRF_all$kap, na.rm = TRUE)
  Cohens_D_RF_kap[i,]<-cohen.d.ci(d=cohenD,n1=Summary_kapRF$count[which(Summary_kapRF$group==coupletolook[1])], n2=Summary_kapRF$count[which(Summary_kapRF$group==coupletolook[2])]) 
}
stats_res_RF_kap<-c()
stats_res_RF_kap$tukey<-c("a","a","b","c")
stats_res_RF_kap$cohen<-c("a","a","b","b")
stats_res_RF_kap<-as.data.frame(stats_res_RF_kap)

anovaGLM_TSS<-aov(TSS ~ group, data = evalmetGLM_all)
# summary(anovaGLM)
tuk_GLM_TSS<-TukeyHSD(anovaGLM_TSS,ordered=TRUE)$group # a b c d 
Cohens_D_GLM_TSS<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM_TSS)<-c("lower","effect","upper")
rownames(Cohens_D_GLM_TSS)<-rownames(tuk_GLM_TSS)
for (i in 1:nrow(tuk_GLM_TSS)){#i=1
  coupletolook <- word(rownames(tuk_GLM_TSS)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSGLM$mean[which(Summary_TSSGLM$group==coupletolook[1])] - Summary_TSSGLM$mean[which(Summary_TSSGLM$group==coupletolook[2])])/sd(evalmetGLM_all$TSS, na.rm = TRUE)
  Cohens_D_GLM_TSS[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSGLM$count[which(Summary_TSSGLM$group==coupletolook[1])], n2=Summary_TSSGLM$count[which(Summary_TSSGLM$group==coupletolook[2])]) 
}
stats_res_GLM_TSS<-c()
stats_res_GLM_TSS$tukey<-c("a","b","c","d")
stats_res_GLM_TSS$cohen<-c("a","ab","b","c")
stats_res_GLM_TSS<-as.data.frame(stats_res_GLM_TSS)

anovaGAM_TSS<-aov(TSS ~ group, data = evalmetGAM_all)
# summary(anovaGAM)
tuk_GAM_TSS<-TukeyHSD(anovaGAM_TSS,ordered=TRUE)$group # a b c d 
Cohens_D_GAM_TSS<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM_TSS)<-c("lower","effect","upper")
rownames(Cohens_D_GAM_TSS)<-rownames(tuk_GAM_TSS)
for (i in 1:nrow(tuk_GAM_TSS)){#i=1
  coupletolook <- word(rownames(tuk_GAM_TSS)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSGAM$mean[which(Summary_TSSGAM$group==coupletolook[1])] - Summary_TSSGAM$mean[which(Summary_TSSGAM$group==coupletolook[2])])/sd(evalmetGAM_all$TSS, na.rm = TRUE)
  Cohens_D_GAM_TSS[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSGAM$count[which(Summary_TSSGAM$group==coupletolook[1])], n2=Summary_TSSGAM$count[which(Summary_TSSGAM$group==coupletolook[2])]) 
}
stats_res_GAM_TSS<-c()
stats_res_GAM_TSS$tukey<-c("a","b","c","d")
stats_res_GAM_TSS$cohen<-c("a","a","b","c")
stats_res_GAM_TSS<-as.data.frame(stats_res_GAM_TSS)

anovaGBM_TSS<-aov(TSS ~ group, data = evalmetGBM_all)
# summary(anovaGBM)
tuk_GBM_TSS<-TukeyHSD(anovaGBM_TSS,ordered=TRUE)$group # a b c d 
Cohens_D_GBM_TSS<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM_TSS)<-c("lower","effect","upper")
rownames(Cohens_D_GBM_TSS)<-rownames(tuk_GBM_TSS)
for (i in 1:nrow(tuk_GBM_TSS)){#i=1
  coupletolook <- word(rownames(tuk_GBM_TSS)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSGBM$mean[which(Summary_TSSGBM$group==coupletolook[1])] - Summary_TSSGBM$mean[which(Summary_TSSGBM$group==coupletolook[2])])/sd(evalmetGBM_all$TSS, na.rm = TRUE)
  Cohens_D_GBM_TSS[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSGBM$count[which(Summary_TSSGBM$group==coupletolook[1])], n2=Summary_TSSGBM$count[which(Summary_TSSGBM$group==coupletolook[2])]) 
}
stats_res_GBM_TSS<-c()
stats_res_GBM_TSS$tukey<-c("a","b","c","d")
stats_res_GBM_TSS$cohen<-c("a","ab","b","c")
stats_res_GBM_TSS<-as.data.frame(stats_res_GBM_TSS)

anovaRF_TSS<-aov(TSS ~ group, data = evalmetRF_all)
# summary(anovaRF)
tuk_RF_TSS<-TukeyHSD(anovaRF_TSS,ordered=TRUE)$group # a b c d 
Cohens_D_RF_TSS<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF_TSS)<-c("lower","effect","upper")
rownames(Cohens_D_RF_TSS)<-rownames(tuk_RF_TSS)
for (i in 1:nrow(tuk_RF_TSS)){#i=1
  coupletolook <- word(rownames(tuk_RF_TSS)[i],1:2,sep = "-")
  cohenD<-(Summary_TSSRF$mean[which(Summary_TSSRF$group==coupletolook[1])] - Summary_TSSRF$mean[which(Summary_TSSRF$group==coupletolook[2])])/sd(evalmetRF_all$TSS, na.rm = TRUE)
  Cohens_D_RF_TSS[i,]<-cohen.d.ci(d=cohenD,n1=Summary_TSSRF$count[which(Summary_TSSRF$group==coupletolook[1])], n2=Summary_TSSRF$count[which(Summary_TSSRF$group==coupletolook[2])]) 
}
stats_res_RF_TSS<-c()
stats_res_RF_TSS$tukey<-c("a","a","b","c")
stats_res_RF_TSS$cohen<-c("a","a","b","c")
stats_res_RF_TSS<-as.data.frame(stats_res_RF_TSS)


###plots###
# 
# Cohens_D[rownames(Cohens_D)=="AR-BA",2]
# print("GLM")
# print(Cohens_D_GLM)
# print(tuk_GLM)
# c("a","ab","b","c")


plotopts<-c(y_title_size=10,y_text_size=8,x_text_size=8,stattext_size=4)
plotGLM<-ggplot(obj_to_plot_GLM, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_GLM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GLM, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GLM")
plotGLM_auc<-ggplot(obj_to_plot_GLM, aes(x=group, y=auc, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("AUC") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0.5,1) + 
  # geom_text(data=stats_res_GLM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GLM_auc, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GLM")

plotGLM_kap<-ggplot(obj_to_plot_GLM, aes(x=group, y=kap, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Kappa") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_GLM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GLM_kap, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GLM")
plotGLM_TSS<-ggplot(obj_to_plot_GLM, aes(x=group, y=TSS, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_GLM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GLM_TSS, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GLM")




plotGAM<-ggplot(obj_to_plot_GAM, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_GAM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GAM, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GAM")
plotGAM_auc<-ggplot(obj_to_plot_GAM, aes(x=group, y=auc, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("AUC") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0.5,1) + 
  # geom_text(data=stats_res_GAM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GAM_auc, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GAM")

plotGAM_kap<-ggplot(obj_to_plot_GAM, aes(x=group, y=kap, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Kappa") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_GAM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GAM_kap, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GAM")
plotGAM_TSS<-ggplot(obj_to_plot_GAM, aes(x=group, y=TSS, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_GAM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GAM_TSS, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GAM")



plotGBM<-ggplot(obj_to_plot_GBM, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_GBM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GBM, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GBM")
plotGBM_auc<-ggplot(obj_to_plot_GBM, aes(x=group, y=auc, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("AUC") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0.5,1) + 
  # geom_text(data=stats_res_GBM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GBM_auc, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GBM")

plotGBM_kap<-ggplot(obj_to_plot_GBM, aes(x=group, y=kap, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Kappa") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_GBM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GBM_kap, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GBM")
plotGBM_TSS<-ggplot(obj_to_plot_GBM, aes(x=group, y=TSS, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_GBM, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_GBM_TSS, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("GBM")



plotRF<-ggplot(obj_to_plot_RF, aes(x=group, y=TSS_adj, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS_adj") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed", col= "red") +
  # geom_text(data=stats_res_RF, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_RF, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("RF")
plotRF_auc<-ggplot(obj_to_plot_RF, aes(x=group, y=auc, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("AUC") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0.5,1) + 
  # geom_text(data=stats_res_RF, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_RF_auc, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("RF")

plotRF_kap<-ggplot(obj_to_plot_RF, aes(x=group, y=kap, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Kappa") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_RF, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_RF_kap, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("RF")
plotRF_TSS<-ggplot(obj_to_plot_RF, aes(x=group, y=TSS, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("TSS") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["x_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  ylim(0,1) + 
  # geom_text(data=stats_res_RF, aes(x=1:4,y=1,label=tukey),col="black") + 
  geom_text(data=stats_res_RF_TSS, aes(x=1:4,y=1,label=cohen),col="black",cex=plotopts["stattext_size"])+ ggtitle("RF")


pdf(file="figures/PAAB_selection/distribTSSadj/distribsTSS.pdf")
plotGLM
plotGAM
plotGBM
plotRF
dev.off()

png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsTSS_GLM.png"),res=300,width=1961,height=1500)
plotGLM
dev.off()
png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsTSS_GAM.png"),res=300,width=1961,height=1500)
plotGAM
dev.off()
png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsTSS_GBM.png"),res=300,width=1961,height=1500)
plotGBM
dev.off()
png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsTSS_RF.png"),res=300,width=1961,height=1500)
plotRF
dev.off()

png(file=paste0("figures/PAAB_selection/distribTSSadj/distribs_all.png"),res=300,width=1961,height=3000)
grid.arrange(plotGLM,plotGAM,plotGBM,plotRF,plotGLM_TSS,plotGAM_TSS,plotGBM_TSS,plotRF_TSS,plotGLM_kap,plotGAM_kap,plotGBM_kap,plotRF_kap,plotGLM_auc,plotGAM_auc,plotGBM_auc,plotRF_auc,nrow=4)
dev.off()



#BA166
load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum
#GLM
load("PA/BA_166/data/GLM/Eval.Rda")
evalmetGLM_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGLM_BA1662<-evalmetGLM_BA166[1:19]
evalmetGLM_BA1662$group<-"BA"
evalmetGLM_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGLM_AR1662<-evalmetGLM_AR166[1:19]
evalmetGLM_AR1662$group<-"AR"
#GAM
load("PA/BA_166/data/GAM/Eval.Rda")
evalmetGAM_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGAM_BA1662<-evalmetGAM_BA166[1:19]
evalmetGAM_BA1662$group<-"BA"
evalmetGAM_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGAM_AR1662<-evalmetGAM_AR166[1:19]
evalmetGAM_AR1662$group<-"AR"
#GBM
load("PA/BA_166/data/GBM/Eval.Rda")
evalmetGBM_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGBM_BA1662<-evalmetGBM_BA166[1:19]
evalmetGBM_BA1662$group<-"BA"
evalmetGBM_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGBM_AR1662<-evalmetGBM_AR166[1:19]
evalmetGBM_AR1662$group<-"AR"
#RF
load("PA/BA_166/data/RF/Eval.Rda")
evalmetRF_BA166<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetRF_BA1662<-evalmetRF_BA166[1:19]
evalmetRF_BA1662$group<-"BA"
evalmetRF_AR166<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetRF_AR1662<-evalmetRF_AR166[1:19]
evalmetRF_AR1662$group<-"AR"

boxplot(evalmetGLM_BA166$TSS_adj,evalmetGLM_BA$TSS_adj)
boxplot(evalmetGLM_AR166$TSS_adj,evalmetGLM_AR$TSS_adj)
t.test(evalmetGLM_BA166$TSS_adj,evalmetGLM_BA$TSS_adj,paired=TRUE)
t.test(evalmetGLM_AR166$TSS_adj,evalmetGLM_AR$TSS_adj,paired=TRUE)
cohenD<-(mean(evalmetGLM_BA$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_BA166$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_BA166$TSS_adj,evalmetGLM_BA$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_BA$TSS_adj))), n2=sum(!(is.na(evalmetGLM_BA166$TSS_adj))))
Cohens_D
cohenD<-(mean(evalmetGLM_AR$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_AR166$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_AR166$TSS_adj,evalmetGLM_AR$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_AR$TSS_adj))), n2=sum(!(is.na(evalmetGLM_AR166$TSS_adj))))
Cohens_D


t.test(evalmetGLM_BA166$TSS_adj,evalmetGLM_PR$TSS_adj)
t.test(evalmetGLM_AR166$TSS_adj,evalmetGLM_PR$TSS_adj)
cohenD<-(mean(evalmetGLM_BA166$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_PR$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_BA166$TSS_adj,evalmetGLM_PR$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_BA166$TSS_adj))), n2=sum(!(is.na(evalmetGLM_PR$TSS_adj))))
Cohens_D
cohenD<-(mean(evalmetGLM_AR166$TSS_adj,na.rm = TRUE) - mean(evalmetGLM_PR$TSS_adj,na.rm = TRUE))/sd(c(evalmetGLM_AR166$TSS_adj,evalmetGLM_PR$TSS_adj), na.rm = TRUE)
Cohens_D<-cohen.d.ci(d=cohenD,n1=sum(!(is.na(evalmetGLM_AR166$TSS_adj))), n2=sum(!(is.na(evalmetGLM_PR$TSS_adj))))
Cohens_D




















###################################################
#rel. ab

#GLM
load("AB/PR/data/GLM/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetGLM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(evalmet)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

evalmetGLM_PR2<-evalmetGLM_PR[1:16]
evalmetGLM_PR2$group<-"PR"
# summary(Eval)
#GAM
load("AB/PR/data/GAM/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetGAM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(evalmet)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGAM_PR2<-evalmetGAM_PR[1:16]
evalmetGAM_PR2$group<-"PR"

#GBM
load("AB/PR/data/GBM/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetGBM_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(evalmet)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetGBM_PR2<-evalmetGBM_PR[1:16]
evalmetGBM_PR2$group<-"PR"


#RF
load("AB/PR/data/RF/Eval_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
evalmetRF_PR<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",-which(colnames(evalmet)=="OTU")],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])
evalmetRF_PR2<-evalmetRF_PR[1:16]
evalmetRF_PR2$group<-"PR"


#FU
#GLM
load("AB/FU/data/GLM/Eval_Met.Rda")
evalmetGLM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",-which(colnames(evalmet)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGLM_FU2<-evalmetGLM_FU[1:16]
evalmetGLM_FU2$group<-"FU"
#GAM
load("AB/FU/data/GAM/Eval_Met.Rda")
evalmetGAM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",-which(colnames(evalmet)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGAM_FU2<-evalmetGAM_FU[1:16]
evalmetGAM_FU2$group<-"FU"
#GBM
load("AB/FU/data/GBM/Eval_Met.Rda")
evalmetGBM_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",-which(colnames(evalmet)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetGBM_FU2<-evalmetGBM_FU[1:16]
evalmetGBM_FU2$group<-"FU"
#RF
load("AB/FU/data/RF/Eval_Met.Rda")
evalmetRF_FU<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",-which(colnames(evalmet)=="OTU")],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
evalmetRF_FU2<-evalmetRF_FU[1:16]
evalmetRF_FU2$group<-"FU"

#BA
#GLM
load("AB/BA/data/GLM/Eval_Met.Rda")
evalmetGLM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGLM_BA2<-evalmetGLM_BA[1:16]
evalmetGLM_BA2$group<-"BA"
evalmetGLM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGLM_AR2<-evalmetGLM_AR[1:16]
evalmetGLM_AR2$group<-"AR"
#GAM
load("AB/BA/data/GAM/Eval_Met.Rda")
evalmetGAM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGAM_BA2<-evalmetGAM_BA[1:16]
evalmetGAM_BA2$group<-"BA"
evalmetGAM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGAM_AR2<-evalmetGAM_AR[1:16]
evalmetGAM_AR2$group<-"AR"
#GBM
load("AB/BA/data/GBM/Eval_Met.Rda")
evalmetGBM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetGBM_BA2<-evalmetGBM_BA[1:16]
evalmetGBM_BA2$group<-"BA"
evalmetGBM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetGBM_AR2<-evalmetGBM_AR[1:16]
evalmetGBM_AR2$group<-"AR"
#RF
load("AB/BA/data/RF/Eval_Met.Rda")
evalmetRF_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
evalmetRF_BA2<-evalmetRF_BA[1:16]
evalmetRF_BA2$group<-"BA"
evalmetRF_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",-which(colnames(evalmet)=="OTU")],BAtaxo[BAtaxo$Kingdom=="Archaea",])
evalmetRF_AR2<-evalmetRF_AR[1:16]
evalmetRF_AR2$group<-"AR"

evalmetGLM_all<-rbind(evalmetGLM_PR2,evalmetGLM_FU2,evalmetGLM_BA2,evalmetGLM_AR2)
obj_to_plot_GLM<-evalmetGLM_all[!(is.na(evalmetGLM_all$Dspear)),c("Dspear","RMSEs","group")]
obj_to_plot_GLM$RMSEs<-log10(obj_to_plot_GLM$RMSEs)
obj_to_plot_GLM$Dspear<-sapply(obj_to_plot_GLM$Dspear,function(X){ifelse(X<0,0,X)})
evalmetGAM_all<-rbind(evalmetGAM_PR2,evalmetGAM_FU2,evalmetGAM_BA2,evalmetGAM_AR2)
obj_to_plot_GAM<-evalmetGAM_all[!(is.na(evalmetGAM_all$Dspear)),c("Dspear","RMSEs","group")]
obj_to_plot_GAM$RMSEs<-sapply(obj_to_plot_GAM$RMSEs,function(X){ifelse(X>100000,100000,X)})
obj_to_plot_GAM$RMSEs<-log10(obj_to_plot_GAM$RMSEs)
obj_to_plot_GAM$Dspear<-sapply(obj_to_plot_GAM$Dspear,function(X){ifelse(X<0,0,X)})

evalmetGBM_all<-rbind(evalmetGBM_PR2,evalmetGBM_FU2,evalmetGBM_BA2,evalmetGBM_AR2)
obj_to_plot_GBM<-evalmetGBM_all[!(is.na(evalmetGBM_all$Dspear)),c("Dspear","RMSEs","group")]
obj_to_plot_GBM$RMSEs<-log10(obj_to_plot_GBM$RMSEs)
obj_to_plot_GBM$Dspear<-sapply(obj_to_plot_GBM$Dspear,function(X){ifelse(X<0,0,X)})
evalmetRF_all<-rbind(evalmetRF_PR2,evalmetRF_FU2,evalmetRF_BA2,evalmetRF_AR2)
obj_to_plot_RF<-evalmetRF_all[!(is.na(evalmetRF_all$Dspear)),c("Dspear","RMSEs","group")]
obj_to_plot_RF$RMSEs<-log10(obj_to_plot_RF$RMSEs)
obj_to_plot_RF$Dspear<-sapply(obj_to_plot_RF$Dspear,function(X){ifelse(X<0,0,X)})
sum(evalmetGLM_all$Dspear,na.rm=TRUE)/length(evalmetGLM_all$Dspear)
sum(evalmetGLM_all$Dspear>0.4,na.rm=TRUE)/length(evalmetGLM_all$Dspear)
sum(evalmetGLM_all$Dspear>0.2,na.rm=TRUE)/length(evalmetGLM_all$Dspear)
sum(evalmetGLM_all$RMSEs,na.rm=TRUE)/length(evalmetGLM_all$RMSEs)



Summary_DspearGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(Dspear, na.rm = TRUE),
    sd = sd(Dspear, na.rm = TRUE)
  )
Summary_DspearGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(Dspear, na.rm = TRUE),
    sd = sd(Dspear, na.rm = TRUE)
  )
Summary_DspearGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(Dspear, na.rm = TRUE),
    sd = sd(Dspear, na.rm = TRUE)
  )
Summary_DspearRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(Dspear, na.rm = TRUE),
    sd = sd(Dspear, na.rm = TRUE)
  )
Summary_CoVGLM<-group_by(evalmetGLM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(RMSEs, na.rm = TRUE),
    sd = sd(RMSEs, na.rm = TRUE)
  )
Summary_CoVGAM<-group_by(evalmetGAM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(RMSEs, na.rm = TRUE),
    sd = sd(RMSEs, na.rm = TRUE)
  )
Summary_CoVGBM<-group_by(evalmetGBM_all, group) %>%
  summarise(
    count = n(),
    mean = mean(RMSEs, na.rm = TRUE),
    sd = sd(RMSEs, na.rm = TRUE)
  )
Summary_CoVRF<-group_by(evalmetRF_all, group) %>%
  summarise(
    count = n(),
    mean = mean(RMSEs, na.rm = TRUE),
    sd = sd(RMSEs, na.rm = TRUE)
  )

#supp able model quality PA part
tobindPRGLM<-evalmetGLM_PR[c("Seq","Phylum","Dspear","RMSEs")]
tobindPRGLM$marker<-"18S_protist"
tobindPRGLM<-tobindPRGLM[c("Seq","marker","Phylum","Dspear","RMSEs")]
tobindFUGLM<-evalmetGLM_FU[c("Seq","Phylum","Dspear","RMSEs")]
tobindFUGLM$marker<-"ITS_fungi"
tobindFUGLM<-tobindFUGLM[c("Seq","marker","Phylum","Dspear","RMSEs")]
tobindBAGLM<-evalmetGLM_BA[c("Seq","Phylum","Dspear","RMSEs")]
tobindBAGLM$marker<-"16S_bacteria-archaea"
tobindBAGLM<-tobindBAGLM[c("Seq","marker","Phylum","Dspear","RMSEs")]
tobindARGLM<-evalmetGLM_AR[c("Seq","Phylum","Dspear","RMSEs")]
tobindARGLM$marker<-"16S_bacteria-archaea"
tobindARGLM<-tobindARGLM[c("Seq","marker","Phylum","Dspear","RMSEs")]

tobindPRGAM<-evalmetGAM_PR[c("Dspear","RMSEs")]
tobindFUGAM<-evalmetGAM_FU[c("Dspear","RMSEs")]
tobindBAGAM<-evalmetGAM_BA[c("Dspear","RMSEs")]
tobindARGAM<-evalmetGAM_AR[c("Dspear","RMSEs")]
tobindPRGBM<-evalmetGBM_PR[c("Dspear","RMSEs")]
tobindFUGBM<-evalmetGBM_FU[c("Dspear","RMSEs")]
tobindBAGBM<-evalmetGBM_BA[c("Dspear","RMSEs")]
tobindARGBM<-evalmetGBM_AR[c("Dspear","RMSEs")]
tobindPRRF<-evalmetRF_PR[c("Dspear","RMSEs")]
tobindFURF<-evalmetRF_FU[c("Dspear","RMSEs")]
tobindBARF<-evalmetRF_BA[c("Dspear","RMSEs")]
tobindARRF<-evalmetRF_AR[c("Dspear","RMSEs")]


evalmet_allPR<-cbind(tobindPRGLM,tobindPRGAM,tobindPRGBM,tobindPRRF)
colnames(evalmet_allPR)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-CoV","GAM-Dspear","GAM-CoV","GBM-Dspear","GBM-CoV","RF-Dspear","RF-CoV")
evalmet_allFU<-cbind(tobindFUGLM,tobindFUGAM,tobindFUGBM,tobindFURF)
colnames(evalmet_allFU)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-CoV","GAM-Dspear","GAM-CoV","GBM-Dspear","GBM-CoV","RF-Dspear","RF-CoV")
evalmet_allBA<-cbind(tobindBAGLM,tobindBAGAM,tobindBAGBM,tobindBARF)
colnames(evalmet_allBA)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-CoV","GAM-Dspear","GAM-CoV","GBM-Dspear","GBM-CoV","RF-Dspear","RF-CoV")
evalmet_allAR<-cbind(tobindARGLM,tobindARGAM,tobindARGBM,tobindARRF)
colnames(evalmet_allAR)<-c("Seq","marker","Phylum","GLM-Dspear","GLM-CoV","GAM-Dspear","GAM-CoV","GBM-Dspear","GBM-CoV","RF-Dspear","RF-CoV")
evalmet_allall<-rbind(evalmet_allPR,evalmet_allFU,evalmet_allBA,evalmet_allAR)
# write.csv(evalmet_allall, file="figures/PAAB_selection/Figshare_tables/Model_quality_phylotype_AB.csv")             

###Statstics###
anovaGLM_Dspear<-aov(Dspear ~ group, data = evalmetGLM_all)
anovaGLM_CoV<-aov(RMSEs ~ group, data = evalmetGLM_all)
# summary(anovaGLM)
tuk_GLM_Dspear<-TukeyHSD(anovaGLM_Dspear,ordered=TRUE)$group # a a b c 
stats_tuk_GLM_Dspear<-c("a","a","b","c")
tuk_GLM_CoV<-TukeyHSD(anovaGLM_CoV,ordered=TRUE)$group  # ab a b ab 
stats_tuk_GLM_CoV<-c("ab","a","b","ab")

Cohens_D_GLM_Dspear<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM_Dspear)<-c("lower","effect","upper")
rownames(Cohens_D_GLM_Dspear)<-rownames(tuk_GLM_Dspear)
for (i in 1:nrow(tuk_GLM_Dspear)){#i=1
  coupletolook <- word(rownames(tuk_GLM_Dspear)[i],1:2,sep = "-")
  cohenD<-(Summary_DspearGLM$mean[which(Summary_DspearGLM$group==coupletolook[1])] - Summary_DspearGLM$mean[which(Summary_DspearGLM$group==coupletolook[2])])/sd(evalmetGLM_all$Dspear, na.rm = TRUE)
  Cohens_D_GLM_Dspear[i,]<-cohen.d.ci(d=cohenD,n1=Summary_DspearGLM$count[which(Summary_DspearGLM$group==coupletolook[1])], n2=Summary_DspearGLM$count[which(Summary_DspearGLM$group==coupletolook[2])]) 
} # a a b c
stats_cohen_GLM_Dspear<-c("a", "a","b", "c")
Cohens_D_GLM_CoV<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GLM_CoV)<-c("lower","effect","upper")
rownames(Cohens_D_GLM_CoV)<-rownames(tuk_GLM_CoV)
for (i in 1:nrow(tuk_GLM_CoV)){#i=1
  coupletolook <- word(rownames(tuk_GLM_CoV)[i],1:2,sep = "-")
  cohenD<-(Summary_CoVGLM$mean[which(Summary_CoVGLM$group==coupletolook[1])] - Summary_CoVGLM$mean[which(Summary_CoVGLM$group==coupletolook[2])])/sd(evalmetGLM_all$RMSEs, na.rm = TRUE)
  Cohens_D_GLM_CoV[i,]<-cohen.d.ci(d=cohenD,n1=Summary_CoVGLM$count[which(Summary_CoVGLM$group==coupletolook[1])], n2=Summary_CoVGLM$count[which(Summary_CoVGLM$group==coupletolook[2])]) 
} # a a a a 
stats_cohen_GLM_CoV<-c("a", "a","a", "a")
stats_GLM_Dspear<-data.frame(cohen=stats_cohen_GLM_Dspear,tuckey=stats_tuk_GLM_Dspear)
stats_GLM_CoV<-data.frame(cohen=stats_cohen_GLM_CoV,tuckey=stats_tuk_GLM_CoV)

anovaGAM_Dspear<-aov(Dspear ~ group, data = evalmetGAM_all)
anovaGAM_CoV<-aov(RMSEs ~ group, data = evalmetGAM_all)
# summary(anovaGAM)
tuk_GAM_Dspear<-TukeyHSD(anovaGAM_Dspear,ordered=TRUE)$group # a a b c 
stats_tuk_GAM_Dspear<-c("a","a","b","c")
tuk_GAM_CoV<-TukeyHSD(anovaGAM_CoV,ordered=TRUE)$group  # a a a a  
stats_tuk_GAM_CoV<-c("a", "a","a", "a")

Cohens_D_GAM_Dspear<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM_Dspear)<-c("lower","effect","upper")
rownames(Cohens_D_GAM_Dspear)<-rownames(tuk_GAM_Dspear)
for (i in 1:nrow(tuk_GAM_Dspear)){#i=1
  coupletolook <- word(rownames(tuk_GAM_Dspear)[i],1:2,sep = "-")
  cohenD<-(Summary_DspearGAM$mean[which(Summary_DspearGAM$group==coupletolook[1])] - Summary_DspearGAM$mean[which(Summary_DspearGAM$group==coupletolook[2])])/sd(evalmetGAM_all$Dspear, na.rm = TRUE)
  Cohens_D_GAM_Dspear[i,]<-cohen.d.ci(d=cohenD,n1=Summary_DspearGAM$count[which(Summary_DspearGAM$group==coupletolook[1])], n2=Summary_DspearGAM$count[which(Summary_DspearGAM$group==coupletolook[2])]) 
} # a a b c
stats_cohen_GAM_Dspear<-c("a", "a","b", "c")
Cohens_D_GAM_CoV<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GAM_CoV)<-c("lower","effect","upper")
rownames(Cohens_D_GAM_CoV)<-rownames(tuk_GAM_CoV)
for (i in 1:nrow(tuk_GAM_CoV)){#i=1
  coupletolook <- word(rownames(tuk_GAM_CoV)[i],1:2,sep = "-")
  cohenD<-(Summary_CoVGAM$mean[which(Summary_CoVGAM$group==coupletolook[1])] - Summary_CoVGAM$mean[which(Summary_CoVGAM$group==coupletolook[2])])/sd(evalmetGAM_all$RMSEs, na.rm = TRUE)
  Cohens_D_GAM_CoV[i,]<-cohen.d.ci(d=cohenD,n1=Summary_CoVGAM$count[which(Summary_CoVGAM$group==coupletolook[1])], n2=Summary_CoVGAM$count[which(Summary_CoVGAM$group==coupletolook[2])]) 
} # a a a a $
stats_cohen_GAM_CoV<-c("a", "a","a", "a")
stats_GAM_Dspear<-data.frame(cohen=stats_cohen_GAM_Dspear,tuckey=stats_tuk_GAM_Dspear)
stats_GAM_CoV<-data.frame(cohen=stats_cohen_GAM_CoV,tuckey=stats_tuk_GAM_CoV)


anovaGBM_Dspear<-aov(Dspear ~ group, data = evalmetGBM_all)
anovaGBM_CoV<-aov(RMSEs ~ group, data = evalmetGBM_all)
# summary(anovaGBM)
tuk_GBM_Dspear<-TukeyHSD(anovaGBM_Dspear,ordered=TRUE)$group # ab a b b
stats_tuk_GBM_Dspear<-c("ab", "a","b", "b")
tuk_GBM_CoV<-TukeyHSD(anovaGBM_CoV,ordered=TRUE)$group  # a b a a 
stats_tuk_GBM_CoV<-c("a", "b","a", "a")

Cohens_D_GBM_Dspear<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM_Dspear)<-c("lower","effect","upper")
rownames(Cohens_D_GBM_Dspear)<-rownames(tuk_GBM_Dspear)
for (i in 1:nrow(tuk_GBM_Dspear)){#i=1
  coupletolook <- word(rownames(tuk_GBM_Dspear)[i],1:2,sep = "-")
  cohenD<-(Summary_DspearGBM$mean[which(Summary_DspearGBM$group==coupletolook[1])] - Summary_DspearGBM$mean[which(Summary_DspearGBM$group==coupletolook[2])])/sd(evalmetGBM_all$Dspear, na.rm = TRUE)
  Cohens_D_GBM_Dspear[i,]<-cohen.d.ci(d=cohenD,n1=Summary_DspearGBM$count[which(Summary_DspearGBM$group==coupletolook[1])], n2=Summary_DspearGBM$count[which(Summary_DspearGBM$group==coupletolook[2])]) 
} # a a a a
stats_cohen_GBM_Dspear<-c("a", "a","a", "a")
Cohens_D_GBM_CoV<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_GBM_CoV)<-c("lower","effect","upper")
rownames(Cohens_D_GBM_CoV)<-rownames(tuk_GBM_CoV)
for (i in 1:nrow(tuk_GBM_CoV)){#i=1
  coupletolook <- word(rownames(tuk_GBM_CoV)[i],1:2,sep = "-")
  cohenD<-(Summary_CoVGBM$mean[which(Summary_CoVGBM$group==coupletolook[1])] - Summary_CoVGBM$mean[which(Summary_CoVGBM$group==coupletolook[2])])/sd(evalmetGBM_all$RMSEs, na.rm = TRUE)
  Cohens_D_GBM_CoV[i,]<-cohen.d.ci(d=cohenD,n1=Summary_CoVGBM$count[which(Summary_CoVGBM$group==coupletolook[1])], n2=Summary_CoVGBM$count[which(Summary_CoVGBM$group==coupletolook[2])]) 
} # a b ab ab
stats_cohen_GBM_CoV<-c("a", "b","ab", "ab")
stats_GBM_Dspear<-data.frame(cohen=stats_cohen_GBM_Dspear,tuckey=stats_tuk_GBM_Dspear)
stats_GBM_CoV<-data.frame(cohen=stats_cohen_GBM_CoV,tuckey=stats_tuk_GBM_CoV)


anovaRF_Dspear<-aov(Dspear ~ group, data = evalmetRF_all)
anovaRF_CoV<-aov(RMSEs ~ group, data = evalmetRF_all)
# summary(anovaRF)
tuk_RF_Dspear<-TukeyHSD(anovaRF_Dspear,ordered=TRUE)$group # a a b c 
stats_tuk_RF_Dspear<-c("a", "a","b", "c")
tuk_RF_CoV<-TukeyHSD(anovaRF_CoV,ordered=TRUE)$group  # a b b c 
stats_tuk_RF_CoV<-c("a", "b","b", "c")

Cohens_D_RF_Dspear<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF_Dspear)<-c("lower","effect","upper")
rownames(Cohens_D_RF_Dspear)<-rownames(tuk_RF_Dspear)
for (i in 1:nrow(tuk_RF_Dspear)){#i=1
  coupletolook <- word(rownames(tuk_RF_Dspear)[i],1:2,sep = "-")
  cohenD<-(Summary_DspearRF$mean[which(Summary_DspearRF$group==coupletolook[1])] - Summary_DspearRF$mean[which(Summary_DspearRF$group==coupletolook[2])])/sd(evalmetRF_all$Dspear, na.rm = TRUE)
  Cohens_D_RF_Dspear[i,]<-cohen.d.ci(d=cohenD,n1=Summary_DspearRF$count[which(Summary_DspearRF$group==coupletolook[1])], n2=Summary_DspearRF$count[which(Summary_DspearRF$group==coupletolook[2])]) 
} # a a b c
stats_cohen_RF_Dspear<-c("a", "a","b", "c")
Cohens_D_RF_CoV<-matrix(NA,nrow=6,ncol=3)
colnames(Cohens_D_RF_CoV)<-c("lower","effect","upper")
rownames(Cohens_D_RF_CoV)<-rownames(tuk_RF_CoV)
for (i in 1:nrow(tuk_RF_CoV)){#i=1
  coupletolook <- word(rownames(tuk_RF_CoV)[i],1:2,sep = "-")
  cohenD<-(Summary_CoVRF$mean[which(Summary_CoVRF$group==coupletolook[1])] - Summary_CoVRF$mean[which(Summary_CoVRF$group==coupletolook[2])])/sd(evalmetRF_all$RMSEs, na.rm = TRUE)
  Cohens_D_RF_CoV[i,]<-cohen.d.ci(d=cohenD,n1=Summary_CoVRF$count[which(Summary_CoVRF$group==coupletolook[1])], n2=Summary_CoVRF$count[which(Summary_CoVRF$group==coupletolook[2])]) 
} # a a a a 
stats_cohen_RF_CoV<-c("a","a","a", "a")
stats_RF_Dspear<-data.frame(cohen=stats_cohen_RF_Dspear,tuckey=stats_tuk_RF_Dspear)
stats_RF_CoV<-data.frame(cohen=stats_cohen_RF_CoV,tuckey=stats_tuk_RF_CoV)


plotopts<-c(y_title_size=10,y_text_size=8,x_text_size=8,stattext_size=4)
plotGLM_Dspear<-ggplot(obj_to_plot_GLM, aes(x=group, y=Dspear, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Spearman correlation") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) + ggtitle("GLM") +
  geom_text(data=stats_GLM_Dspear, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])

plotGLM_CoV<-ggplot(obj_to_plot_GLM, aes(x=group, y=RMSEs, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("log10(Coefficint of Variation)") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  # geom_text(data=stats_res_GLM, aes(x=1:4,y=1,label=tukey),col="black") + 
  ggtitle("GLM")+
  geom_hline(yintercept=0, linetype="dashed", col= "red") + geom_hline(yintercept=-1, linetype="dashed", col= "black")+ geom_hline(yintercept=-2, linetype="dashed", col= "blue") +
  geom_text(data=stats_GLM_CoV, aes(x=1:4,y=4,label=cohen),col="black",cex=plotopts["stattext_size"])
#dashed line red : mean error is as high as mean of observation ; black : mean error is 10% of mean of observation ; blue : mean error is 1% of mean of obsevration

plotGAM_Dspear<-ggplot(obj_to_plot_GAM, aes(x=group, y=Dspear, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Spearman correlation") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none")  + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) + ggtitle("GAM") +
  geom_text(data=stats_GAM_Dspear, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])

plotGAM_CoV<-ggplot(obj_to_plot_GAM, aes(x=group, y=RMSEs, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("log10(Coefficint of Variation)") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  # geom_text(data=stats_res_GAM, aes(x=1:4,y=1,label=tukey),col="black") + 
  ggtitle("GAM")+
  geom_hline(yintercept=0, linetype="dashed", col= "red") + geom_hline(yintercept=-1, linetype="dashed", col= "black")+ geom_hline(yintercept=-2, linetype="dashed", col= "blue") +
  geom_text(data=stats_GAM_CoV, aes(x=1:4,y=4,label=cohen),col="black",cex=plotopts["stattext_size"])
  
#dashed line red : mean error is as high as mean of observation ; black : mean error is 10% of mean of observation ; blue : mean error is 1% of mean of obsevration

plotGBM_Dspear<-ggplot(obj_to_plot_GBM, aes(x=group, y=Dspear, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Spearman correlation") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) + ggtitle("GBM") +
  geom_text(data=stats_GBM_Dspear, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])

plotGBM_CoV<-ggplot(obj_to_plot_GBM, aes(x=group, y=RMSEs, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("log10(Coefficint of Variation)") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  # geom_text(data=stats_res_GBM, aes(x=1:4,y=1,label=tukey),col="black") + 
  ggtitle("GBM")+
  geom_hline(yintercept=0, linetype="dashed", col= "red") + geom_hline(yintercept=-1, linetype="dashed", col= "black")+ geom_hline(yintercept=-2, linetype="dashed", col= "blue") +
  geom_text(data=stats_GBM_CoV, aes(x=1:4,y=4,label=cohen),col="black",cex=plotopts["stattext_size"])
#dashed line red : mean error is as high as mean of observation ; black : mean error is 10% of mean of observation ; blue : mean error is 1% of mean of obsevration

plotRF_Dspear<-ggplot(obj_to_plot_RF, aes(x=group, y=Dspear, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("Spearman correlation") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) + ggtitle("RF")  +
  geom_text(data=stats_RF_Dspear, aes(x=1:4,y=0.95,label=cohen),col="black",cex=plotopts["stattext_size"])

plotRF_CoV<-ggplot(obj_to_plot_RF, aes(x=group, y=RMSEs, colour=group))+ 
  geom_violin(aes(fill=group), scale="width",alpha=0.2) +
  ylab("log10(Coefficint of Variation)") + geom_boxplot(width=0.2,lwd=0.2,pch=1, outlier.shape = NA) + labs(x="") + 
  theme_bw() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(margin=margin(r=10),size=plotopts["y_title_size"]), 
                     axis.text.y = element_text(size=plotopts["y_text_size"]),
                     axis.text.x = element_text(angle = 30, vjust=0.5,colour=c('#FFD700', '#5CB800','#6B4C62','#457EB0'),size=plotopts["y_text_size"]), 
                     plot.margin = unit(c(0,0,-.5,0.2),"cm"),
                     legend.position = "none") + 
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_fill_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  scale_x_discrete(labels=c("Archaea", "Bacteria", "Fungi","Protist")) +
  # geom_text(data=stats_res_RF, aes(x=1:4,y=1,label=tukey),col="black") + 
  ggtitle("RF")+
  geom_hline(yintercept=0, linetype="dashed", col= "red") + geom_hline(yintercept=-1, linetype="dashed", col= "black")+ geom_hline(yintercept=-2, linetype="dashed", col= "blue") +
  geom_text(data=stats_RF_CoV, aes(x=1:4,y=4,label=cohen),col="black",cex=plotopts["stattext_size"])
#dashed line red : mean error is as high as mean of observation ; black : mean error is 10% of mean of observation ; blue : mean error is 1% of mean of obsevration
png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsDspear_all.png"),res=300,width=1961,height=1500)
grid.arrange(plotGLM_Dspear,plotGAM_Dspear,plotGBM_Dspear,plotRF_Dspear,nrow=2)
dev.off()
png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsCoV_all.png"),res=300,width=1961,height=1500)
grid.arrange(plotGLM_CoV,plotGAM_CoV,plotGBM_CoV,plotRF_CoV,nrow=2)
dev.off()
png(file=paste0("figures/PAAB_selection/distribTSSadj/distribsDspearCoV_all.png"),res=300,width=1961,height=1500)
grid.arrange(plotGLM_Dspear,plotGAM_Dspear,plotGBM_Dspear,plotRF_Dspear,plotGLM_CoV,plotGAM_CoV,plotGBM_CoV,plotRF_CoV,nrow=2)
dev.off()