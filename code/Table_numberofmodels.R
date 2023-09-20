#create plots of number of models still passing step X
#best plots : Proportion_successful_models3 &

library(tidyverse)
library(ggplot2)
library(reshape2) 
library(grid)
library(gridExtra)
dataset=c(BA='#FFD700',AR='#5CB800',FU='#6B4C62',PR='#457EB0')

load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda"))
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.



load("PA/PR/data/GLM/Eval.Rda")
load("PA/PR/data/GLM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGLM<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGLM<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

# load(paste0("PA/PR/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/PR/data/GAM/Eval.Rda")
load("PA/PR/data/GAM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGAM<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGAM<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

load("PA/PR/data/GBM/Eval.Rda")
load("PA/PR/data/GBM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGBM<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGBM<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

load("PA/PR/data/RF/Eval.Rda")
load("PA/PR/data/RF/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(Eval)){
  Eval<- Eval[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetRF<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetRF<-cbind(Eval[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


# Result table, number of model per group
#number of ASV before modelling
nASV_PR<-c(length(evalmetGLM$TSS),
           length(evalmetGAM$TSS),
           length(evalmetGBM$TSS),
           length(evalmetRF$TSS))
#number of modelled (crashed models and prev>0.95 excluded)
nMod_PR<-c(sum(!(is.na(evalmetGLM$TSS))),
           sum(!(is.na(evalmetGAM$TSS))),
           sum(!(is.na(evalmetGBM$TSS))),
           sum(!(is.na(evalmetRF$TSS))))

#successfully modelled by all 4 algos
successAll_PR<-sum(apply(cbind(!(is.na(evalmetGLM$TSS)),!(is.na(evalmetGAM$TSS)),!(is.na(evalmetGBM$TSS)),!(is.na(evalmetRF$TSS))),1,sum)==4)
#successfully modelled by at least 1 algo
successOne_PR<-sum(apply(cbind(!(is.na(evalmetGLM$TSS)),!(is.na(evalmetGAM$TSS)),!(is.na(evalmetGBM$TSS)),!(is.na(evalmetRF$TSS))),1,sum)!=0)

#Number of nodel for which the fit is >0.2 (TSS)
nFit02_PR<-c(sum(fitmetGLM$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGAM$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGBM$TSS>=0.2,na.rm=TRUE),
             sum(fitmetRF$TSS>=0.2,na.rm=TRUE))
nEval02_PR<-c(sum(evalmetGLM$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGAM$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGBM$TSS>=0.2,na.rm=TRUE),
              sum(evalmetRF$TSS>=0.2,na.rm=TRUE))
nEval05_PR<-c(sum(evalmetGLM$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGAM$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGBM$TSS>=0.5,na.rm=TRUE),
              sum(evalmetRF$TSS>=0.5,na.rm=TRUE))
#TSS>05 for all 4 algos
TSS05_All_PR<-sum(apply(cbind(evalmetGLM$TSS>0.5,evalmetGAM$TSS>0.5,evalmetGBM$TSS>0.5,evalmetRF$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
TSS05_One_PR<-sum(apply(cbind(evalmetGLM$TSS>0.5,evalmetGAM$TSS>0.5,evalmetGBM$TSS>0.5,evalmetRF$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})!=0)

#number of models for which the fit is > 95th percentile of null models
nEval_null_PR<-c(sum(evalmetGLM$TSS_sign,na.rm=TRUE),
              sum(evalmetGAM$TSS_sign,na.rm=TRUE),
              sum(evalmetGBM$TSS_sign,na.rm=TRUE),
              sum(evalmetRF$TSS_sign,na.rm=TRUE))
#sum(evalmetGLM$TSS_sign,na.rm=TRUE)/sum(!(is.na(evalmetGLM$TSS)))
nEval_adj02_PR<-c(sum(evalmetGLM$TSS_adj>=0.2,na.rm=TRUE),
                 sum(evalmetGAM$TSS_adj>=0.2,na.rm=TRUE),
                 sum(evalmetGBM$TSS_adj>=0.2,na.rm=TRUE),
                 sum(evalmetRF$TSS_adj>=0.2,na.rm=TRUE))
nEval_adj04_PR<-c(sum(evalmetGLM$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.4,na.rm=TRUE))
nEval_adj06_PR<-c(sum(evalmetGLM$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.6,na.rm=TRUE))
nEval_adj08_PR<-c(sum(evalmetGLM$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.8,na.rm=TRUE))
nEval_adj05_PR<-c(sum(evalmetGLM$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.5,na.rm=TRUE))
#sum(evalmetGLM$TSS_adj>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM$TSS)))
#better than null for all 4 algos
Better_null_All_PR<-sum(apply(cbind(evalmetGLM$TSS_sign,evalmetGAM$TSS_sign,evalmetGBM$TSS_sign,evalmetRF$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
Better_null_One_PR<-sum(apply(cbind(evalmetGLM$TSS_sign,evalmetGAM$TSS_sign,evalmetGBM$TSS_sign,evalmetRF$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})!=0)


##############################################FU###########################

load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.

load("PA/FU/data/GLM/Eval.Rda")
load("PA/FU/data/GLM/Fit_Met.Rda")
fitmetGLM<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGLM<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])
# load(paste0("PA/FU/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/FU/data/GAM/Eval.Rda")
load("PA/FU/data/GAM/Fit_Met.Rda")
fitmetGAM<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGAM<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("PA/FU/data/GBM/Eval.Rda")
load("PA/FU/data/GBM/Fit_Met.Rda")
fitmetGBM<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGBM<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("PA/FU/data/RF/Eval.Rda")
load("PA/FU/data/RF/Fit_Met.Rda")
fitmetRF<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetRF<-cbind(Eval[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

# Result table, number of model per group
#number of ASV before modelling
nASV_FU<-c(length(fitmetGLM$TSS),
           length(fitmetGAM$TSS),
           length(fitmetGBM$TSS),
           length(fitmetRF$TSS))
#number of modelled (crashed models and prev>0.95 excluded)
nMod_FU<-c(sum(!(is.na(fitmetGLM$TSS))),
           sum(!(is.na(fitmetGAM$TSS))),
           sum(!(is.na(fitmetGBM$TSS))),
           sum(!(is.na(fitmetRF$TSS))))
#successfully modelled by all 4 algos
successAll_FU<-sum(apply(cbind(!(is.na(fitmetGLM$TSS)),!(is.na(fitmetGAM$TSS)),!(is.na(fitmetGBM$TSS)),!(is.na(fitmetRF$TSS))),1,sum)==4)
#successfully modelled by at least 1 algo
successOne_FU<-sum(apply(cbind(!(is.na(fitmetGLM$TSS)),!(is.na(fitmetGAM$TSS)),!(is.na(fitmetGBM$TSS)),!(is.na(fitmetRF$TSS))),1,sum)!=0)

#Number of nodel for which the fit is >0.2 (TSS)
nFit02_FU<-c(sum(fitmetGLM$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGAM$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGBM$TSS>=0.2,na.rm=TRUE),
             sum(fitmetRF$TSS>=0.2,na.rm=TRUE))
nEval02_FU<-c(sum(evalmetGLM$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGAM$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGBM$TSS>=0.2,na.rm=TRUE),
              sum(evalmetRF$TSS>=0.2,na.rm=TRUE))
nEval05_FU<-c(sum(evalmetGLM$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGAM$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGBM$TSS>=0.5,na.rm=TRUE),
              sum(evalmetRF$TSS>=0.5,na.rm=TRUE))
#better than null for all 4 algos
TSS05_All_FU<-sum(apply(cbind(evalmetGLM$TSS>0.5,evalmetGAM$TSS>0.5,evalmetGBM$TSS>0.5,evalmetRF$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
TSS05_One_FU<-sum(apply(cbind(evalmetGLM$TSS>0.5,evalmetGAM$TSS>0.5,evalmetGBM$TSS>0.5,evalmetRF$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})!=0)

nEval_null_FU<-c(sum(evalmetGLM$TSS_sign,na.rm=TRUE),
                 sum(evalmetGAM$TSS_sign,na.rm=TRUE),
                 sum(evalmetGBM$TSS_sign,na.rm=TRUE),
                 sum(evalmetRF$TSS_sign,na.rm=TRUE))
#sum(evalmetGLM$TSS_sign,na.rm=TRUE)/sum(!(is.na(evalmetGLM$TSS))) #81% #pourcentage better than null
nEval_adj02_FU<-c(sum(evalmetGLM$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.2,na.rm=TRUE))
nEval_adj04_FU<-c(sum(evalmetGLM$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.4,na.rm=TRUE))
nEval_adj06_FU<-c(sum(evalmetGLM$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.6,na.rm=TRUE))
nEval_adj08_FU<-c(sum(evalmetGLM$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.8,na.rm=TRUE))
nEval_adj05_FU<-c(sum(evalmetGLM$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGAM$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGBM$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetRF$TSS_adj>=0.5,na.rm=TRUE))
#sum(evalmetGLM$TSS_adj>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM$TSS)))
#better than null for all 4 algos
Better_null_All_FU<-sum(apply(cbind(evalmetGLM$TSS_sign,evalmetGAM$TSS_sign,evalmetGBM$TSS_sign,evalmetRF$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
Better_null_One_FU<-sum(apply(cbind(evalmetGLM$TSS_sign,evalmetGAM$TSS_sign,evalmetGBM$TSS_sign,evalmetRF$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})!=0)



##############################################BAAR###########################

load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum

load("PA/BA/data/GLM/Eval.Rda")
load("PA/BA/data/GLM/Fit_Met.Rda")

evalmetGLM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/BA/data/GAM/Eval.Rda")
load("PA/BA/data/GAM/Fit_Met.Rda")
evalmetGAM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("PA/BA/data/GBM/Eval.Rda")
load("PA/BA/data/GBM/Fit_Met.Rda")
evalmetGBM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


load("PA/BA/data/RF/Eval.Rda")
load("PA/BA/data/RF/Fit_Met.Rda")
evalmetRF_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


# Result table, number of model per group
#number of ASV before modelling
nASV_BA<-c(length(fitmetGLM_BA$TSS),
           length(fitmetGAM_BA$TSS),
           length(fitmetGBM_BA$TSS),
           length(fitmetRF_BA$TSS))
length(fitmetGLM_AR$TSS)
#number of modelled (crashed models and prev>0.95 excluded)
nMod_BA<-c(sum(!(is.na(fitmetGLM_BA$TSS))),
           sum(!(is.na(fitmetGAM_BA$TSS))),
           sum(!(is.na(fitmetGBM_BA$TSS))),
           sum(!(is.na(fitmetRF_BA$TSS))))
#successfully modelled by all 4 algos
successAll_BA<-sum(apply(cbind(!(is.na(fitmetGLM_BA$TSS)),!(is.na(fitmetGAM_BA$TSS)),!(is.na(fitmetGBM_BA$TSS)),!(is.na(fitmetRF_BA$TSS))),1,sum)==4)
successAll_AR<-sum(apply(cbind(!(is.na(fitmetGLM_AR$TSS)),!(is.na(fitmetGAM_AR$TSS)),!(is.na(fitmetGBM_AR$TSS)),!(is.na(fitmetRF_AR$TSS))),1,sum)==4)
#successfully modelled by at least 1 algo
successOne_BA<-sum(apply(cbind(!(is.na(fitmetGLM_BA$TSS)),!(is.na(fitmetGAM_BA$TSS)),!(is.na(fitmetGBM_BA$TSS)),!(is.na(fitmetRF_BA$TSS))),1,sum)!=0)
successOne_AR<-sum(apply(cbind(!(is.na(fitmetGLM_AR$TSS)),!(is.na(fitmetGAM_AR$TSS)),!(is.na(fitmetGBM_AR$TSS)),!(is.na(fitmetRF_AR$TSS))),1,sum)!=0)


#Number of nodel for which the fit is >0.2 (TSS)
nFit02_BA<-c(sum(fitmetGLM_BA$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGAM_BA$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGBM_BA$TSS>=0.2,na.rm=TRUE),
             sum(fitmetRF_BA$TSS>=0.2,na.rm=TRUE))
nEval02_BA<-c(sum(evalmetGLM_BA$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGAM_BA$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGBM_BA$TSS>=0.2,na.rm=TRUE),
              sum(evalmetRF_BA$TSS>=0.2,na.rm=TRUE))
nEval05_BA<-c(sum(evalmetGLM_BA$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGAM_BA$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGBM_BA$TSS>=0.5,na.rm=TRUE),
              sum(evalmetRF_BA$TSS>=0.5,na.rm=TRUE))
#better than null for all 4 algos
TSS05_All_BA<-sum(apply(cbind(evalmetGLM_BA$TSS>0.5,evalmetGAM_BA$TSS>0.5,evalmetGBM_BA$TSS>0.5,evalmetRF_BA$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
TSS05_One_BA<-sum(apply(cbind(evalmetGLM_BA$TSS>0.5,evalmetGAM_BA$TSS>0.5,evalmetGBM_BA$TSS>0.5,evalmetRF_BA$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})!=0)

nEval_null_BA<-c(sum(evalmetGLM_BA$TSS_sign,na.rm=TRUE),
                 sum(evalmetGAM_BA$TSS_sign,na.rm=TRUE),
                 sum(evalmetGBM_BA$TSS_sign,na.rm=TRUE),
                 sum(evalmetRF_BA$TSS_sign,na.rm=TRUE))
#sum(evalmetGLM_BA$TSS_sign,na.rm=TRUE)/sum(!(is.na(evalmetGLM_BA$TSS))) #91% #pourcentage better than null
#sum(evalmetGLM_AR$TSS_sign,na.rm=TRUE)/sum(!(is.na(evalmetGLM_AR$TSS))) #98% #pourcentage better than null
nEval_adj02_BA<-c(sum(evalmetGLM_BA$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetGAM_BA$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetGBM_BA$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetRF_BA$TSS_adj>=0.2,na.rm=TRUE))
nEval_adj04_BA<-c(sum(evalmetGLM_BA$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetGAM_BA$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetGBM_BA$TSS_adj>=0.4,na.rm=TRUE),
                  sum(evalmetRF_BA$TSS_adj>=0.4,na.rm=TRUE))
nEval_adj06_BA<-c(sum(evalmetGLM_BA$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetGAM_BA$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetGBM_BA$TSS_adj>=0.6,na.rm=TRUE),
                  sum(evalmetRF_BA$TSS_adj>=0.6,na.rm=TRUE))
nEval_adj08_BA<-c(sum(evalmetGLM_BA$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetGAM_BA$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetGBM_BA$TSS_adj>=0.8,na.rm=TRUE),
                  sum(evalmetRF_BA$TSS_adj>=0.8,na.rm=TRUE))
nEval_adj05_BA<-c(sum(evalmetGLM_BA$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGAM_BA$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGBM_BA$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetRF_BA$TSS_adj>=0.5,na.rm=TRUE))
#sum(evalmetGLM_BA$TSS_adj>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM_BA$TSS)))
# Result table, number of model per group
#number of ASV before modelling
nASV_AR<-c(length(fitmetGLM_AR$TSS),
           length(fitmetGAM_AR$TSS),
           length(fitmetGBM_AR$TSS),
           length(fitmetRF_AR$TSS))

#number of modelled (crashed models and prev>0.95 excluded)
nMod_AR<-c(sum(!(is.na(fitmetGLM_AR$TSS))),
           sum(!(is.na(fitmetGAM_AR$TSS))),
           sum(!(is.na(fitmetGBM_AR$TSS))),
           sum(!(is.na(fitmetRF_AR$TSS))))

#Number of nodel for which the fit is >0.2 (TSS)
nFit02_AR<-c(sum(fitmetGLM_AR$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGAM_AR$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGBM_AR$TSS>=0.2,na.rm=TRUE),
             sum(fitmetRF_AR$TSS>=0.2,na.rm=TRUE))
nEval02_AR<-c(sum(evalmetGLM_AR$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGAM_AR$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGBM_AR$TSS>=0.2,na.rm=TRUE),
              sum(evalmetRF_AR$TSS>=0.2,na.rm=TRUE))
nEval05_AR<-c(sum(evalmetGLM_AR$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGAM_AR$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGBM_AR$TSS>=0.5,na.rm=TRUE),
              sum(evalmetRF_AR$TSS>=0.5,na.rm=TRUE))
#better than null for all 4 algos
TSS05_All_AR<-sum(apply(cbind(evalmetGLM_AR$TSS>0.5,evalmetGAM_AR$TSS>0.5,evalmetGBM_AR$TSS>0.5,evalmetRF_AR$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
TSS05_One_AR<-sum(apply(cbind(evalmetGLM_AR$TSS>0.5,evalmetGAM_AR$TSS>0.5,evalmetGBM_AR$TSS>0.5,evalmetRF_AR$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})!=0)

nEval_null_AR<-c(sum(evalmetGLM_AR$TSS_sign,na.rm=TRUE),
                 sum(evalmetGAM_AR$TSS_sign,na.rm=TRUE),
                 sum(evalmetGBM_AR$TSS_sign,na.rm=TRUE),
                 sum(evalmetRF_AR$TSS_sign,na.rm=TRUE))
nEval_adj02_AR<-c(sum(evalmetGLM_AR$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetGAM_AR$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetGBM_AR$TSS_adj>=0.2,na.rm=TRUE),
                  sum(evalmetRF_AR$TSS_adj>=0.2,na.rm=TRUE))
nEval_adj04_AR<-c(sum(evalmetGLM_AR$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.4,na.rm=TRUE))
nEval_adj06_AR<-c(sum(evalmetGLM_AR$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.6,na.rm=TRUE))
nEval_adj08_AR<-c(sum(evalmetGLM_AR$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.8,na.rm=TRUE))
nEval_adj05_AR<-c(sum(evalmetGLM_AR$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGAM_AR$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetGBM_AR$TSS_adj>=0.5,na.rm=TRUE),
                  sum(evalmetRF_AR$TSS_adj>=0.5,na.rm=TRUE))
#sum(evalmetGLM_AR$TSS_adj>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM_AR$TSS)))
#better than null for all 4 algos
Better_null_All_BA<-sum(apply(cbind(evalmetGLM_BA$TSS_sign,evalmetGAM_BA$TSS_sign,evalmetGBM_BA$TSS_sign,evalmetRF_BA$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
Better_null_One_BA<-sum(apply(cbind(evalmetGLM_BA$TSS_sign,evalmetGAM_BA$TSS_sign,evalmetGBM_BA$TSS_sign,evalmetRF_BA$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})!=0)
#better than null for all 4 algos
Better_null_All_AR<-sum(apply(cbind(evalmetGLM_AR$TSS_sign,evalmetGAM_AR$TSS_sign,evalmetGBM_AR$TSS_sign,evalmetRF_AR$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
Better_null_One_AR<-sum(apply(cbind(evalmetGLM_AR$TSS_sign,evalmetGAM_AR$TSS_sign,evalmetGBM_AR$TSS_sign,evalmetRF_AR$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})!=0)

# sum(evalmetGLM_AR[which(evalmetGLM_AR$Phylum=="Thaumarchaeota"),]$TSS_adj>=0.4)

##############################################BA166###########################

load("PA/BA_166/data/GLM/Eval.Rda")
load("PA/BA_166/data/GLM/Fit_Met.Rda")
evalmetGLM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("PA/BA_166/data/GAM/Eval.Rda")
load("PA/BA_166/data/GAM/Fit_Met.Rda")
evalmetGAM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("PA/BA_166/data/GBM/Eval.Rda")
load("PA/BA_166/data/GBM/Fit_Met.Rda")
evalmetGBM_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


load("PA/BA_166/data/RF/Eval.Rda")
load("PA/BA_166/data/RF/Fit_Met.Rda")
evalmetRF_BA<-cbind(Eval[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR<-cbind(Eval[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]


# Result table, number of model per group
#number of ASV before modelling
nASV_BA_166<-c(length(fitmetGLM_BA$TSS),
           length(fitmetGAM_BA$TSS),
           length(fitmetGBM_BA$TSS),
           length(fitmetRF_BA$TSS))

#number of modelled (crashed models and prev>0.95 excluded)
nMod_BA_166<-c(sum(!(is.na(fitmetGLM_BA$TSS))),
           sum(!(is.na(fitmetGAM_BA$TSS))),
           sum(!(is.na(fitmetGBM_BA$TSS))),
           sum(!(is.na(fitmetRF_BA$TSS))))

#successfully modelled by all 4 algos
successAll_BA166<-sum(apply(cbind(!(is.na(fitmetGLM_BA$TSS)),!(is.na(fitmetGAM_BA$TSS)),!(is.na(fitmetGBM_BA$TSS)),!(is.na(fitmetRF_BA$TSS))),1,sum)==4)
successAll_AR166<-sum(apply(cbind(!(is.na(fitmetGLM_AR$TSS)),!(is.na(fitmetGAM_AR$TSS)),!(is.na(fitmetGBM_AR$TSS)),!(is.na(fitmetRF_AR$TSS))),1,sum)==4)
#successfully modelled by at least 1 algo
successOne_BA166<-sum(apply(cbind(!(is.na(fitmetGLM_BA$TSS)),!(is.na(fitmetGAM_BA$TSS)),!(is.na(fitmetGBM_BA$TSS)),!(is.na(fitmetRF_BA$TSS))),1,sum)!=0)
successOne_AR166<-sum(apply(cbind(!(is.na(fitmetGLM_AR$TSS)),!(is.na(fitmetGAM_AR$TSS)),!(is.na(fitmetGBM_AR$TSS)),!(is.na(fitmetRF_AR$TSS))),1,sum)!=0)


#Number of nodel for which the fit is >0.2 (TSS)
nFit02_BA_166<-c(sum(fitmetGLM_BA$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGAM_BA$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGBM_BA$TSS>=0.2,na.rm=TRUE),
             sum(fitmetRF_BA$TSS>=0.2,na.rm=TRUE))
nEval02_BA_166<-c(sum(evalmetGLM_BA$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGAM_BA$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGBM_BA$TSS>=0.2,na.rm=TRUE),
              sum(evalmetRF_BA$TSS>=0.2,na.rm=TRUE))
nEval05_BA_166<-c(sum(evalmetGLM_BA$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGAM_BA$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGBM_BA$TSS>=0.5,na.rm=TRUE),
              sum(evalmetRF_BA$TSS>=0.5,na.rm=TRUE))
#better than null for all 4 algos
TSS05_All_BA166<-sum(apply(cbind(evalmetGLM_BA$TSS>0.5,evalmetGAM_BA$TSS>0.5,evalmetGBM_BA$TSS>0.5,evalmetRF_BA$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
TSS05_One_BA166<-sum(apply(cbind(evalmetGLM_BA$TSS>0.5,evalmetGAM_BA$TSS>0.5,evalmetGBM_BA$TSS>0.5,evalmetRF_BA$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})!=0)


nEval_null_BA_166<-c(sum(evalmetGLM_BA$TSS_sign,na.rm=TRUE),
                 sum(evalmetGAM_BA$TSS_sign,na.rm=TRUE),
                 sum(evalmetGBM_BA$TSS_sign,na.rm=TRUE),
                 sum(evalmetRF_BA$TSS_sign,na.rm=TRUE))
nEval_adj02_BA_166<-c(sum(evalmetGLM_BA$TSS_adj>=0.2,na.rm=TRUE),
                     sum(evalmetGAM_BA$TSS_adj>=0.2,na.rm=TRUE),
                     sum(evalmetGBM_BA$TSS_adj>=0.2,na.rm=TRUE),
                     sum(evalmetRF_BA$TSS_adj>=0.2,na.rm=TRUE))
nEval_adj04_BA_166<-c(sum(evalmetGLM_BA$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetGAM_BA$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetGBM_BA$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetRF_BA$TSS_adj>=0.4,na.rm=TRUE))
nEval_adj06_BA_166<-c(sum(evalmetGLM_BA$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetGAM_BA$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetGBM_BA$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetRF_BA$TSS_adj>=0.6,na.rm=TRUE))
nEval_adj08_BA_166<-c(sum(evalmetGLM_BA$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetGAM_BA$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetGBM_BA$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetRF_BA$TSS_adj>=0.8,na.rm=TRUE))

# Result table, number of model per group
#number of ASV before modelling
nASV_AR_166<-c(length(fitmetGLM_AR$TSS),
           length(fitmetGAM_AR$TSS),
           length(fitmetGBM_AR$TSS),
           length(fitmetRF_AR$TSS))

#number of modelled (crashed models and prev>0.95 excluded)
nMod_AR_166<-c(sum(!(is.na(fitmetGLM_AR$TSS))),
           sum(!(is.na(fitmetGAM_AR$TSS))),
           sum(!(is.na(fitmetGBM_AR$TSS))),
           sum(!(is.na(fitmetRF_AR$TSS))))

#Number of nodel for which the fit is >0.2 (TSS)
nFit02_AR_166<-c(sum(fitmetGLM_AR$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGAM_AR$TSS>=0.2,na.rm=TRUE),
             sum(fitmetGBM_AR$TSS>=0.2,na.rm=TRUE),
             sum(fitmetRF_AR$TSS>=0.2,na.rm=TRUE))
nEval02_AR_166<-c(sum(evalmetGLM_AR$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGAM_AR$TSS>=0.2,na.rm=TRUE),
              sum(evalmetGBM_AR$TSS>=0.2,na.rm=TRUE),
              sum(evalmetRF_AR$TSS>=0.2,na.rm=TRUE))
nEval05_AR_166<-c(sum(evalmetGLM_AR$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGAM_AR$TSS>=0.5,na.rm=TRUE),
              sum(evalmetGBM_AR$TSS>=0.5,na.rm=TRUE),
              sum(evalmetRF_AR$TSS>=0.5,na.rm=TRUE))
#better than null for all 4 algos
TSS05_All_AR166<-sum(apply(cbind(evalmetGLM_AR$TSS>0.5,evalmetGAM_AR$TSS>0.5,evalmetGBM_AR$TSS>0.5,evalmetRF_AR$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
TSS05_One_AR166<-sum(apply(cbind(evalmetGLM_AR$TSS>0.5,evalmetGAM_AR$TSS>0.5,evalmetGBM_AR$TSS>0.5,evalmetRF_AR$TSS>0.5),1,function(X){sum(X,na.rm=TRUE)})!=0)

nEval_null_AR_166<-c(sum(evalmetGLM_AR$TSS_sign,na.rm=TRUE),
                 sum(evalmetGAM_AR$TSS_sign,na.rm=TRUE),
                 sum(evalmetGBM_AR$TSS_sign,na.rm=TRUE),
                 sum(evalmetRF_AR$TSS_sign,na.rm=TRUE))
nEval_adj02_AR_166<-c(sum(evalmetGLM_AR$TSS_adj>=0.2,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.2,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.2,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.2,na.rm=TRUE))
nEval_adj04_AR_166<-c(sum(evalmetGLM_AR$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.4,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.4,na.rm=TRUE))
nEval_adj06_AR_166<-c(sum(evalmetGLM_AR$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.6,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.6,na.rm=TRUE))
nEval_adj08_AR_166<-c(sum(evalmetGLM_AR$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetGAM_AR$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetGBM_AR$TSS_adj>=0.8,na.rm=TRUE),
                      sum(evalmetRF_AR$TSS_adj>=0.8,na.rm=TRUE))

#better than null for all 4 algos
Better_null_All_BA166<-sum(apply(cbind(evalmetGLM_BA$TSS_sign,evalmetGAM_BA$TSS_sign,evalmetGBM_BA$TSS_sign,evalmetRF_BA$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
Better_null_One_BA166<-sum(apply(cbind(evalmetGLM_BA$TSS_sign,evalmetGAM_BA$TSS_sign,evalmetGBM_BA$TSS_sign,evalmetRF_BA$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})!=0)
#better than null for all 4 algos
Better_null_All_AR166<-sum(apply(cbind(evalmetGLM_AR$TSS_sign,evalmetGAM_AR$TSS_sign,evalmetGBM_AR$TSS_sign,evalmetRF_AR$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})==4)
#better than nul for at least 1 algo
Better_null_One_AR166<-sum(apply(cbind(evalmetGLM_AR$TSS_sign,evalmetGAM_AR$TSS_sign,evalmetGBM_AR$TSS_sign,evalmetRF_AR$TSS_sign),1,function(X){sum(X,na.rm=TRUE)})!=0)




Tablenumbers<-data.frame(Number_of=c("ASV","model successfully fitted by at least 1 algorithm","model successfully fitted by all 4 algorithms","model with better predictive capacity than null models for at least 1 algo","model with better predictive capacity than null models for all 4 algo", "models with TSS>0.5 for at least 1 algo","model with TSS>0.5 for all 4 algos"),
                         Bacteria_250=c(max(nMod_BA),successOne_BA,successAll_BA,Better_null_One_BA,Better_null_All_BA,TSS05_One_BA,TSS05_All_BA),
                         Archaea_250=c(max(nMod_AR),successOne_AR,successAll_AR,Better_null_One_AR,Better_null_All_AR,TSS05_One_AR,TSS05_All_AR),
                         Fungi_217=c(max(nMod_FU),successOne_FU,successAll_FU,Better_null_One_FU,Better_null_All_FU,TSS05_One_FU,TSS05_All_FU),
                         Protist_166=c(max(nMod_PR),successOne_PR,successAll_PR,Better_null_One_PR,Better_null_All_PR,TSS05_One_PR,TSS05_All_PR),
                         Bacteria_166=c(max(nMod_BA_166),successOne_BA166,successAll_BA166,Better_null_One_BA166,Better_null_All_BA166,TSS05_One_BA166,TSS05_All_BA166),
                         Archaea_166=c(max(nMod_AR_166),successOne_AR166,successAll_AR166,Better_null_One_AR166,Better_null_All_AR166,TSS05_One_AR166,TSS05_All_AR166))
write.csv(Tablenumbers,file="figures/PAAB_selection/Table_numberofmodels_PA.csv")

Tableresult<- data.frame(Bacteria_250=c(nMod_BA,nFit02_BA,nEval02_BA,nEval05_BA,nEval_null_BA,nEval_adj02_BA,nEval_adj04_BA,nEval_adj06_BA,nEval_adj08_BA),
                         Bacteria_166=c(nMod_BA_166,nFit02_BA_166,nEval02_BA_166,nEval05_BA_166,nEval_null_BA_166,nEval_adj02_BA_166,nEval_adj04_BA_166,nEval_adj06_BA_166,nEval_adj08_BA_166),
                         Archaea_250=c(nMod_AR,nFit02_AR,nEval02_AR,nEval05_AR,nEval_null_AR,nEval_adj02_AR,nEval_adj04_AR,nEval_adj06_AR,nEval_adj08_AR),
                         Archaea_166=c(nMod_AR_166,nFit02_AR_166,nEval02_AR_166,nEval05_AR_166,nEval_null_AR_166,nEval_adj02_AR_166,nEval_adj04_AR_166,nEval_adj06_AR_166,nEval_adj08_AR_166),
                         Fungi_217=c(nMod_FU,nFit02_FU,nEval02_FU,nEval05_FU,nEval_null_FU,nEval_adj02_FU,nEval_adj04_FU,nEval_adj06_FU,nEval_adj08_FU),
                         Protist_166=c(nMod_PR,nFit02_PR,nEval02_PR,nEval05_PR,nEval_null_PR,nEval_adj02_PR,nEval_adj04_PR,nEval_adj06_PR,nEval_adj08_PR),
                         Model=rep(c("GLM","GAM","GBM","RF"),9),
                         Counted=c(rep("Fitted",4),rep("Fit02",4),rep("Eval02",4),rep("Eval05",4),rep("TSS > null",4),rep("TSSadj > 0.2",4),rep("TSSadj > 0.4",4),rep("TSSadj > 0.6",4),rep("TSSadj > 0.8",4)))
Tablesresult2<- rbind(c(c(length(fitmetGLM_BA$TSS),length(fitmetGLM_BA$TSS),length(fitmetGLM_AR$TSS),length(fitmetGLM_AR$TSS),51863,3160,"data","nASV_data")),
                     c(817,817,0,0,34,8,"data","nprev095"),
                     Tableresult)
Tablesresult2$Bacteria_250<-as.numeric(Tablesresult2$Bacteria_250)
Tablesresult2$Archaea_250<-as.numeric(Tablesresult2$Archaea_250)
Tablesresult2$Fungi_217<-as.numeric(Tablesresult2$Fungi_217)
Tablesresult2$Protist_166<-as.numeric(Tablesresult2$Protist_166)
Tablesresult2$Bacteria_166<-as.numeric(Tablesresult2$Bacteria_166)
Tablesresult2$Archaea_166<-as.numeric(Tablesresult2$Archaea_166)
Tablesresult2$Counted<-factor(Tablesresult2$Counted,levels=unique(Tablesresult2$Counted)[c(3,7,8,9,10,11)])


M_results<-melt(Tablesresult2[!(is.na(Tablesresult2$Counted)),], id= c("Model","Counted"))

levels(M_results$Counted)

pGLM<-ggplot(M_results[M_results$Model%in%c("data","GLM"),], aes(x=Counted, y=value,group=variable,color=variable)) +
  geom_line() + geom_point() + ggtitle("GLM") +scale_color_brewer(palette="Paired")+theme_minimal()+ scale_x_discrete(guide = guide_axis(n.dodge=2))
pGAM<-ggplot(M_results[M_results$Model%in%c("data","GAM"),], aes(x=Counted, y=value,group=variable,color=variable))+ 
  geom_line()+geom_point() + ggtitle("GAM") +scale_color_brewer(palette="Paired")+theme_minimal()+ scale_x_discrete(guide = guide_axis(n.dodge=2))

pGBM<-ggplot(M_results[M_results$Model%in%c("data","GBM"),], aes(x=Counted, y=value,group=variable,color=variable))+
  geom_line()+geom_point() + ggtitle("GBM") +scale_color_brewer(palette="Paired")+theme_minimal()+ scale_x_discrete(guide = guide_axis(n.dodge=2))
pRF<-ggplot(M_results[M_results$Model%in%c("data","RF"),], aes(x=Counted, y=value,group=variable,color=variable))+ 
  geom_line()+geom_point() + ggtitle("RF") +scale_color_brewer(palette="Paired")+theme_minimal()+ scale_x_discrete(guide = guide_axis(n.dodge=2))



plot_numb_passing_mods<-grid.arrange(pGLM,pGAM,pGBM,pRF)  
pdf(file="figures/PAAB_selection/testmodels/Number_successful_models2.pdf")
plot(plot_numb_passing_mods)
dev.off()


#same %
Tablesresult3<-Tablesresult2[!(is.na(Tablesresult2$Counted)),]
Tablesresult4<-Tablesresult2[!(is.na(Tablesresult2$Counted)),]
for (Mod in unique (Tablesresult3$Model)){#Mod="GLM"
Tablesresult3[Tablesresult3$Model==Mod,1:6]<-apply(Tablesresult3[Tablesresult3$Model==Mod,1:6],2,function(X){X/X[1]})
for (i in 1:(nrow(Tablesresult4[Tablesresult4$Model==Mod,1:6])-1)){#i =1
Tablesresult4[Tablesresult4$Model==Mod,1:6][i,]<-(Tablesresult4[Tablesresult4$Model==Mod,1:6][i,] - Tablesresult4[Tablesresult4$Model==Mod,1:6][i+1,])/Tablesresult2[Tablesresult2$Model==Mod,1:6][1,]
}
Tablesresult4[Tablesresult4$Model==Mod,1:6][nrow(Tablesresult4[Tablesresult4$Model==Mod,1:6]),] <- Tablesresult4[Tablesresult4$Model==Mod,1:6][nrow(Tablesresult4[Tablesresult4$Model==Mod,1:6]),]/Tablesresult2[Tablesresult2$Model==Mod,1:6][1,]
}

M_results3 <- melt(Tablesresult3, id= c("Model","Counted"))%>%
  mutate(nsite=ifelse(variable%in%c("Bacteria_166","Archea_166"),"166","250")) %>%
  mutate(variable2=substr(variable,1,8))

# pGLM3 <-ggplot(M_results3[M_results3$Model%in%c("data","GLM"),], aes(x=Counted, y=value )) +
#   geom_point(aes(shape=variable2, color=variable),size=5) +
#   scale_shape_manual(values = c(15,16,17,18)) +
#   scale_color_manual(
#     values = c("#47BF61", "#6EFF82", "#FFBD4A","#FFD155","#3B66FF","#FF4776"),labels = c("Bacteria", "Bacteria_166", "Archaea", "Archaea_166", "Fungi", "Protista_166")
#   ) +
#   geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                      axis.text.x = element_blank(),
#                      axis.text.y = element_text(size=20),
#                      legend.title = element_blank(),legend.text = element_blank()) + 
#   labs(x="") + labs(y="") +
#   guides(shape = "none",
#          color = guide_legend(
#     override.aes=list(shape = c(16,16,15,15,17,18))))
# 
# pGAM3<-ggplot(M_results3[M_results3$Model%in%c("data","GAM"),], aes(x=Counted, y=value )) +
#   geom_point(aes(shape=variable2, color=variable),size=5) +
#   scale_shape_manual(values = c(15,16,17,18)) +
#   scale_color_manual(
#     values = c("#47BF61", "#6EFF82", "#FFBD4A","#FFD155","#3B66FF","#FF4776"),labels = c("Bacteria", "Bacteria_166", "Archaea", "Archaea_166", "Fungi", "Protista_166")
#   ) +
#   geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                      axis.text.x = element_blank(),
#                      axis.text.y = element_text(size=20),
#                      legend.title = element_blank(),legend.text = element_text(size=20)) + 
#   labs(x="") + labs(y="") +
#   guides(shape = "none",
#          color = guide_legend(
#            override.aes=list(shape = c(16,16,15,15,17,18))))
# pGBM3<-ggplot(M_results3[M_results3$Model%in%c("data","GBM"),], aes(x=Counted, y=value )) +
#   geom_point(aes(shape=variable2, color=variable),size=5) +
#   scale_shape_manual(values = c(15,16,17,18)) +
#   scale_color_manual(
#     values = c("#47BF61", "#6EFF82", "#FFBD4A","#FFD155","#3B66FF","#FF4776"),labels = c("Bacteria", "Bacteria_166", "Archaea", "Archaea_166", "Fungi", "Protista_166")
#   ) +
#   geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                      axis.text.x = element_text(angle = 60, vjust=0.5, size=15),
#                      axis.text.y = element_text(size=20),
#                      legend.title = element_blank(),legend.text = element_blank()) + 
#   labs(x="") + labs(y="") +
#   guides(shape = "none",
#          color = guide_legend(
#            override.aes=list(shape = c(16,16,15,15,17,18))))
# pRF3<-ggplot(M_results3[M_results3$Model%in%c("data","RF"),], aes(x=Counted, y=value )) +
#   geom_point(aes(shape=variable2, color=variable),size=5) +
#   scale_shape_manual(values = c(15,16,17,18)) +
#   scale_color_manual(
#     values = c("#47BF61", "#6EFF82", "#FFBD4A","#FFD155","#3B66FF","#FF4776"),labels = c("Bacteria", "Bacteria_166", "Archaea", "Archaea_166", "Fungi", "Protista_166")
#   ) +
#   geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                      axis.text.x = element_text(angle = 60, vjust=0.5, size=15),
#                      axis.text.y = element_text(size=20),
#                      legend.title = element_blank(),legend.text = element_text(size=20)) + 
#   labs(x="") + labs(y="") +
#   guides(shape = "none",
#          color = guide_legend(
#            override.aes=list(shape = c(16,16,15,15,17,18))))
# 
# plot_prop_passing_mods<-grid.arrange(pGLM3,pGAM3,pGBM3,pRF3)  
# pdf(file="figures/testmodels/Proportion_successful_models2.pdf")
# plot(plot_prop_passing_mods)
# dev.off()



Fig_PA_modquality_per_dataset<-ggplot(M_results3, aes(x=Counted, y=value,group=Model)) +
  geom_point(aes(shape=variable2, color=variable),size=2) +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_color_manual(
    values = c('#FFD700', "#E0BF00", '#5CB800',"#66CC00",'#6B4C62','#457EB0'),labels = c("Bacteria", "Bacteria_166", "Archaea", "Archaea_166", "Fungi", "Protista_166")
  ) +
  facet_wrap(~Model) +
  geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(angle = 20, margin = margin(t = 5),hjust=1, size=10),
                     axis.text.y = element_text(size=10),
                     legend.position = c(1,1),
                     legend.justification = c(1, 1),
                     legend.text = element_text(size=5),
                     legend.key.size = unit(0.5,"cm"),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill = NA, colour = NA),
                     strip.text.x = element_text(size = 15)) + 
  labs(x="",y="Proportion of models") +
  guides(shape = "none",
         color = guide_legend(
           override.aes=list(shape = c(16,16,15,15,17,18))))
pdf(file="figures/PAAB_selection/testmodels/Proportion_successful_models_PA_6groups.pdf")
plot(Fig_PA_modquality_per_dataset)
dev.off()

png(file=paste0("figures/PAAB_selection/testmodels/Proportion_successful_models_PA_6groups.png"),res=300,width=1961,height=1500)
plot(Fig_PA_modquality_per_dataset)
dev.off()


Fig_PA_modquality_per_dataset_4group<-ggplot(M_results3[!(M_results3$variable%in%c("Bacteria_166","Archaea_166")),], aes(x=Counted, y=value,group=Model)) +
  geom_point(aes(shape=variable2, color=variable),size=2) +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  facet_wrap(~Model) +
  geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(angle = 20, margin = margin(t = 5),hjust=1, size=10),
                     axis.text.y = element_text(size=10),
                     legend.position = c(1,1),
                     legend.justification = c(1, 1),
                     legend.text = element_text(size=8),
                     legend.key.size = unit(0.4,"cm"),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill = NA, colour = NA),
                     strip.text.x = element_text(size = 15),
                     plot.margin = unit(c(0.2,0.2,0,0.2),"cm")) + 
  labs(x="",y="Proportion of models") +
  guides(shape = "none",
         color = guide_legend(
           override.aes=list(shape = c(16,15,17,18))))
pdf(file="figures/PAAB_selection/testmodels/Proportion_successful_models_PA_4groups.pdf")
plot(Fig_PA_modquality_per_dataset_4group)
dev.off()

png(file=paste0("figures/PAAB_selection/testmodels/Proportion_successful_models_PA_4groups.png"),res=300,width=1961,height=1500)
plot(Fig_PA_modquality_per_dataset_4group)
dev.off()

# library(ggplot2)
# library(hrbrthemes)
# M_results4 <- melt(Tablesresult4, id= c("Model","Counted"))
# ggplot(M_results3, aes(fill=variable, y=value,x=Model,group=variable)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_bar(position="dodge", stat="identity",fill="grey10",aes(alpha=Counted)) +
#   scale_alpha_manual(values=c(0,0.1,0.2,0.3,0.5,0.8),name="Counted",
#                      breaks=c("nMod_Fitted", "nMod_Eval>null", "nMod_TSSadj>0.2", "nMod_TSSadj>0.4","nMod_TSSadj>0.6","nMod_TSSadj>0.8")) 
# 
# ggplot(M_results3[M_results3$variable!="Bacteria_166"&M_results3$variable!="Archaea_166",], aes(fill=variable, y=value,x=Model,group=variable)) +
#   geom_bar(position="dodge", stat="identity") +
#   geom_bar(position="dodge", stat="identity",fill="grey10",aes(alpha=Counted)) +
#   scale_alpha_manual(values=c(0,0.1,0.2,0.3,0.5,0.8),name="Counted",
#                      breaks=c("nMod_Fitted", "nMod_Eval>null", "nMod_TSSadj>0.2", "nMod_TSSadj>0.4","nMod_TSSadj>0.6","nMod_TSSadj>0.8")) 
# 
# ggplot(M_results4[M_results4$variable!="Bacteria_166"&M_results4$variable!="Archaea_166",], aes(fill=Counted, y=value,x=variable)) +
#   geom_bar(position='stack', stat='identity') +
#   facet_wrap(~Model)
# 
# pdf(file="figures/testmodels/Proportion_successful_models3.pdf")
# ggplot(M_results4[M_results4$variable!="Bacteria_166"&M_results4$variable!="Archaea_166",], aes(fill=Counted, y=value,x=variable)) +
#   geom_bar(position='stack', stat='identity') +
#   facet_wrap(~Model)
# dev.off()





############################
##########SAME AB############


#normalised RMSE tests 
#(patching code 17/11/2022, real value should be computed in evalmetrics code).
#here im normalising the mean RMSE of all boots, each RMSE should be normalized before doing the mean to be correct)
#Maybe remove models that have RMSE error of double the amount of max value?
load("../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")
PRtaxo<-PRtaxo_grouped_phylum
load(paste0("PA/PR/data/OTUdata.Rda"))
all(PRtaxo$Seq==colnames(OTUdata)) #check same sequences same order.

# mins<-apply(OTUdata,2,min)
# maxs<-apply(OTUdata,2,max)

load("AB/PR/data/GLM/Eval_Met.Rda")
load("AB/PR/data/GLM/Fit_Met.Rda")
# evalmet<-evalmet[which(!(is.na(evalmet$Dspear))),]
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGLM<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGLM<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


load("AB/PR/data/GAM/Eval_Met.Rda")
load("AB/PR/data/GAM/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGAM<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGAM<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])

summary(evalmetGAM$Dspear)
# summary(evalmetGAM$Dspear)


load("AB/PR/data/GBM/Eval_Met.Rda")
load("AB/PR/data/GBM/Fit_Met.Rda")
# evalmet<-evalmet[which(!(is.na(evalmet$Dspear))),]
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetGBM<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetGBM<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


load("AB/PR/data/RF/Eval_Met.Rda")
load("AB/PR/data/RF/Fit_Met.Rda")
if(nrow(PRtaxo)!=nrow(evalmet)){
  evalmet<- evalmet[-3019,]
}
if(nrow(PRtaxo)!=nrow(fitmet)){
  fitmet<- fitmet[-3019,]
}
fitmetRF<-fitmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",]
evalmetRF<-cbind(evalmet[PRtaxo_grouped_phylum$Phylum!="Not_Protist",],PRtaxo[PRtaxo_grouped_phylum$Phylum!="Not_Protist",])


#remove non protists if not alredy done
# 
#   load("OTUdata10415PR.Rda")
#   load("seqvecPR.Rda")
#   evalmetRF <- evalmetRF[which(colnames(OTUdata_PR_AB)%in%seqvec),]
#   fitmetRF<- fitmetRF[which(colnames(OTUdata_PR_AB)%in%seqvec),]



# Result table, number of model per group
#number of ASV before modelling
nASV_PR<-c(length(fitmetGLM$Dspear),
           length(fitmetGAM$Dspear),
           length(fitmetGBM$Dspear),
           length(fitmetRF$Dspear))
#number of modelled (crashed models and Dpear<0 excluded)
nMod_PR<-c(sum(fitmetGLM$Dspear>0,na.rm=TRUE),
           sum(fitmetGAM$Dspear>0,na.rm=TRUE),
           sum(fitmetGBM$Dspear>0,na.rm=TRUE),
           sum(fitmetRF$Dspear>0,na.rm=TRUE))
success_All_PR<-sum(apply(cbind(!(is.na(evalmetGLM$Dspear)),!(is.na(evalmetGAM$Dspear)),!(is.na(evalmetGBM$Dspear)),!(is.na(evalmetRF$Dspear))),1,sum)==4)
#successfully modelled by at least 1 algo
success_One_PR<-sum(apply(cbind(!(is.na(evalmetGLM$Dspear)),!(is.na(evalmetGAM$Dspear)),!(is.na(evalmetGBM$Dspear)),!(is.na(evalmetRF$Dspear))),1,sum)!=0)

# prevdat<-apply(OTUdata,2,function(X){sum(X>0)/length(X)})
# ncol(OTUdata)
# nrow(fitmetRF)
# plot(prevdat~fitmetRF$Dpear[-3019])
# plot(prevdat~evalmetRF$Dpear[-3019])
# plot((fitmetGLM$Dpear-evalmetGLM$Dpear),prevdat)

# sum(fitmetRF$Dpear<evalmetRF$Dpear,na.rm=TRUE)
# sum(fitmetRF$Dpear>evalmetRF$Dpear,na.rm=TRUE)
# sum(fitmetRF$Dpear[fitmetRF$Dpear>0]<evalmetRF$Dpear[fitmetRF$Dpear>0],na.rm=TRUE)
# sum(fitmetRF$Dpear[fitmetRF$Dpear>0]>evalmetRF$Dpear[fitmetRF$Dpear>0],na.rm=TRUE)
#Number of nodel for which the fit is >0.2 (TSS)
nFit02_PR<-c(sum(fitmetGLM$Dspear>0.2,na.rm=TRUE),
             sum(fitmetGAM$Dspear>0.2,na.rm=TRUE),
             sum(fitmetGBM$Dspear>0.2,na.rm=TRUE),
             sum(fitmetRF$Dspear>0.2,na.rm=TRUE))
nEval02_PR<-c(sum(evalmetGLM$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.2,na.rm=TRUE))
#sum(evalmetGLM$Dspear>=0.2,na.rm=TRUE)/sum(!(is.na(evalmetGLM$Dspear))) #49%
nEval05_PR<-c(sum(evalmetGLM$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.5,na.rm=TRUE))
#sum(evalmetGLM$Dspear>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM$Dspear))) #3%
nEval04_PR<-c(sum(evalmetGLM$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.4,na.rm=TRUE))
#sum(evalmetGLM$Dspear>=0.4,na.rm=TRUE)/sum(!(is.na(evalmetGLM$Dspear))) #9%
nEval06_PR<-c(sum(evalmetGLM$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.6,na.rm=TRUE))
nEval08_PR<-c(sum(evalmetGLM$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.8,na.rm=TRUE))
#good for all 4 algos
GoodMod_All_PR<-sum(apply(cbind(evalmetGLM$Dspear>0.4,evalmetGAM$Dspear>0.4,evalmetGBM$Dspear>0.4,evalmetRF$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})==4)
#good for at least 1 algo
GoodMod_One_PR<-sum(apply(cbind(evalmetGLM$Dspear>0.4,evalmetGAM$Dspear>0.4,evalmetGBM$Dspear>0.4,evalmetRF$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})!=0)


UseMod_All_PR<-sum(apply(cbind(evalmetGLM$Dspear>0.2,evalmetGAM$Dspear>0.2,evalmetGBM$Dspear>0.2,evalmetRF$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})==4)
#successfully modelled by at least 1 algo
UseMod_One_PR<-sum(apply(cbind(evalmetGLM$Dspear>0.2,evalmetGAM$Dspear>0.2,evalmetGBM$Dspear>0.2,evalmetRF$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})!=0)
##############################################FU###########################

load("../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")
load(paste0("PA/FU/data/OTUdata.Rda"))
FUtaxo2021<-FUtaxo_grouped_phylum
all(FUtaxo2021$Seq==colnames(OTUdata)) #check same sequences same order.

load("AB/FU/data/GLM/Eval_Met.Rda")
load("AB/FU/data/GLM/Fit_Met.Rda")
fitmetGLM<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGLM<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("AB/FU/data/GAM/Eval_Met.Rda")
load("AB/FU/data/GAM/Fit_Met.Rda")
fitmetGAM<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGAM<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])

load("AB/FU/data/GBM/Eval_Met.Rda")
load("AB/FU/data/GBM/Fit_Met.Rda")
fitmetGBM<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetGBM<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])


load("AB/FU/data/RF/Eval_Met.Rda")
load("AB/FU/data/RF/Fit_Met.Rda")
fitmetRF<-fitmet[FUtaxo2021$Kingdom=="Fungi",]
evalmetRF<-cbind(evalmet[FUtaxo2021$Kingdom=="Fungi",],FUtaxo2021[FUtaxo2021$Kingdom=="Fungi",])


# Result table, number of model per group
#number of ASV before modelling
nASV_FU<-c(length(fitmetGLM$Dspear),
           length(fitmetGAM$Dspear),
           length(fitmetGBM$Dspear),
           length(fitmetRF$Dspear))
#number of modelled (crashed models and prev>0.95 excluded)
nMod_FU<-c(sum(!(is.na(fitmetGLM$Dspear))),
           sum(!(is.na(fitmetGAM$Dspear))),
           sum(!(is.na(fitmetGBM$Dspear))),
           sum(!(is.na(fitmetRF$Dspear))))
success_All_FU<-sum(apply(cbind(!(is.na(evalmetGLM$Dspear)),!(is.na(evalmetGAM$Dspear)),!(is.na(evalmetGBM$Dspear)),!(is.na(evalmetRF$Dspear))),1,sum)==4)
#successfully modelled by at least 1 algo
success_One_FU<-sum(apply(cbind(!(is.na(evalmetGLM$Dspear)),!(is.na(evalmetGAM$Dspear)),!(is.na(evalmetGBM$Dspear)),!(is.na(evalmetRF$Dspear))),1,sum)!=0)

#Number of nodel for which the fit is >0.2 (TSS)
nFit02_FU<-c(sum(fitmetGLM$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetGAM$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetGBM$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetRF$Dspear>=0.2,na.rm=TRUE))
nEval02_FU<-c(sum(evalmetGLM$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.2,na.rm=TRUE))
#sum(evalmetGLM$Dspear>=0.2,na.rm=TRUE)/sum(!(is.na(evalmetGLM$Dspear))) #63%
nEval05_FU<-c(sum(evalmetGLM$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.5,na.rm=TRUE))
#sum(evalmetGLM$Dspear>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM$Dspear))) #7%

nEval04_FU<-c(sum(evalmetGLM$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.4,na.rm=TRUE))
#sum(evalmetGLM$Dspear>=0.4,na.rm=TRUE)/sum(!(is.na(evalmetGLM$Dspear))) #20%
nEval06_FU<-c(sum(evalmetGLM$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.6,na.rm=TRUE))
nEval08_FU<-c(sum(evalmetGLM$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGAM$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGBM$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetRF$Dspear>=0.8,na.rm=TRUE))

#good for all 4 algos

GoodMod_All_FU<-sum(apply(cbind(evalmetGLM$Dspear>0.4,evalmetGAM$Dspear>0.4,evalmetGBM$Dspear>0.4,evalmetRF$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})==4)
#good for at least 1 algo
GoodMod_One_FU<-sum(apply(cbind(evalmetGLM$Dspear>0.4,evalmetGAM$Dspear>0.4,evalmetGBM$Dspear>0.4,evalmetRF$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})!=0)

UseMod_All_FU<-sum(apply(cbind(evalmetGLM$Dspear>0.2,evalmetGAM$Dspear>0.2,evalmetGBM$Dspear>0.2,evalmetRF$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})==4)
#successfully modelled by at least 1 algo
UseMod_One_FU<-sum(apply(cbind(evalmetGLM$Dspear>0.2,evalmetGAM$Dspear>0.2,evalmetGBM$Dspear>0.2,evalmetRF$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})!=0)
##############################################BA###########################

load(paste0("PA/BA/data/OTUdata.Rda"))
load("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")
all(BAtaxo_grouped_phylum$Seq==colnames(OTUdata)) #check same sequences same order.
BAtaxo<-BAtaxo_grouped_phylum

load("AB/BA/data/GLM/Eval_Met.Rda")
load("AB/BA/data/GLM/Fit_Met.Rda")
evalmetGLM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGLM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGLM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGLM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

# load(paste0("PA/BA/data/OTUdata.Rda"))
# plot(evalmet$TSS~apply(OTUdata,2,function(X){sum(X)/length(X)}))
#pas de lien prevalence - model quality
load("AB/BA/data/GAM/Eval_Met.Rda")
load("AB/BA/data/GAM/Fit_Met.Rda")
evalmetGAM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGAM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGAM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGAM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]



load("AB/BA/data/GBM/Eval_Met.Rda")
load("AB/BA/data/GBM/Fit_Met.Rda")
evalmetGBM_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetGBM_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetGBM_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetGBM_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

load("AB/BA/data/RF/Eval_Met.Rda")
load("AB/BA/data/RF/Fit_Met.Rda")
evalmetRF_BA<-cbind(evalmet[BAtaxo$Kingdom=="Bacteria",],BAtaxo[BAtaxo$Kingdom=="Bacteria",])
fitmetRF_BA<-fitmet[BAtaxo$Kingdom=="Bacteria",]
evalmetRF_AR<-cbind(evalmet[BAtaxo$Kingdom=="Archaea",],BAtaxo[BAtaxo$Kingdom=="Archaea",])
fitmetRF_AR<-fitmet[BAtaxo$Kingdom=="Archaea",]

# Result table, number of model per group
#number of ASV before modelling
nASV_BA<-c(length(fitmetGLM_BA$Dspear),
           length(fitmetGAM_BA$Dspear),
           length(fitmetGBM_BA$Dspear),
           length(fitmetRF_BA$Dspear))

#number of modelled (crashed models and prev>0.95 excluded)
nMod_BA<-c(sum(!(is.na(fitmetGLM_BA$Dspear))),
           sum(!(is.na(fitmetGAM_BA$Dspear))),
           sum(!(is.na(fitmetGBM_BA$Dspear))),
           sum(!(is.na(fitmetRF_BA$Dspear))))

success_All_BA<-sum(apply(cbind(!(is.na(evalmetGLM_BA$Dspear)),!(is.na(evalmetGAM_BA$Dspear)),!(is.na(evalmetGBM_BA$Dspear)),!(is.na(evalmetRF_BA$Dspear))),1,sum)==4)
#successfully modelled by at least 1 algo
success_One_BA<-sum(apply(cbind(!(is.na(evalmetGLM_BA$Dspear)),!(is.na(evalmetGAM_BA$Dspear)),!(is.na(evalmetGBM_BA$Dspear)),!(is.na(evalmetRF_BA$Dspear))),1,sum)!=0)

#Number of nodel for which the fit is >0.2 (TSS)
nFit02_BA<-c(sum(fitmetGLM_BA$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetGAM_BA$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetGBM_BA$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetRF_BA$Dspear>=0.2,na.rm=TRUE))
nEval02_BA<-c(sum(evalmetGLM_BA$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGAM_BA$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGBM_BA$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetRF_BA$Dspear>=0.2,na.rm=TRUE))
#sum(evalmetGLM_BA$Dspear>=0.2,na.rm=TRUE)/sum(!(is.na(evalmetGLM_BA$Dspear))) #85%
nEval05_BA<-c(sum(evalmetGLM_BA$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGAM_BA$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGBM_BA$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetRF_BA$Dspear>=0.5,na.rm=TRUE))
#sum(evalmetGLM_BA$Dspear>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM_BA$Dspear))) #22%
nEval04_BA<-c(sum(evalmetGLM_BA$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGAM_BA$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGBM_BA$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetRF_BA$Dspear>=0.4,na.rm=TRUE))
#sum(evalmetGLM_BA$Dspear>=0.4,na.rm=TRUE)/sum(!(is.na(evalmetGLM_BA$Dspear))) #40%
nEval06_BA<-c(sum(evalmetGLM_BA$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGAM_BA$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGBM_BA$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetRF_BA$Dspear>=0.6,na.rm=TRUE))
nEval08_BA<-c(sum(evalmetGLM_BA$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGAM_BA$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGBM_BA$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetRF_BA$Dspear>=0.8,na.rm=TRUE))

GoodMod_All_BA<-sum(apply(cbind(evalmetGLM_BA$Dspear>0.4,evalmetGAM_BA$Dspear>0.4,evalmetGBM_BA$Dspear>0.4,evalmetRF_BA$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})==4)
#successfully modelled by at least 1 algo
GoodMod_One_BA<-sum(apply(cbind(evalmetGLM_BA$Dspear>0.4,evalmetGAM_BA$Dspear>0.4,evalmetGBM_BA$Dspear>0.4,evalmetRF_BA$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})!=0)

UseMod_All_BA<-sum(apply(cbind(evalmetGLM_BA$Dspear>0.2,evalmetGAM_BA$Dspear>0.2,evalmetGBM_BA$Dspear>0.2,evalmetRF_BA$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})==4)
#successfully modelled by at least 1 algo
UseMod_One_BA<-sum(apply(cbind(evalmetGLM_BA$Dspear>0.2,evalmetGAM_BA$Dspear>0.2,evalmetGBM_BA$Dspear>0.2,evalmetRF_BA$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})!=0)


# Result table, number of model per group
#number of ASV before modelling
nASV_AR<-c(length(fitmetGLM_AR$Dspear),
           length(fitmetGAM_AR$Dspear),
           length(fitmetGBM_AR$Dspear),
           length(fitmetRF_AR$Dspear))

#number of modelled (crashed models and prev>0.95 excluded)
nMod_AR<-c(sum(!(is.na(fitmetGLM_AR$Dspear))),
           sum(!(is.na(fitmetGAM_AR$Dspear))),
           sum(!(is.na(fitmetGBM_AR$Dspear))),
           sum(!(is.na(fitmetRF_AR$Dspear))))
success_All_AR<-sum(apply(cbind(!(is.na(evalmetGLM_AR$Dspear)),!(is.na(evalmetGAM_AR$Dspear)),!(is.na(evalmetGBM_AR$Dspear)),!(is.na(evalmetRF_AR$Dspear))),1,sum)==4)
#successfully modelled by at least 1 algo
success_One_AR<-sum(apply(cbind(!(is.na(evalmetGLM_AR$Dspear)),!(is.na(evalmetGAM_AR$Dspear)),!(is.na(evalmetGBM_AR$Dspear)),!(is.na(evalmetRF_AR$Dspear))),1,sum)!=0)

#Number of nodel for which the fit is >0.2 (Dpear)
nFit02_AR<-c(sum(fitmetGLM_AR$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetGAM_AR$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetGBM_AR$Dspear>=0.2,na.rm=TRUE),
             sum(fitmetRF_AR$Dspear>=0.2,na.rm=TRUE))
nEval02_AR<-c(sum(evalmetGLM_AR$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGAM_AR$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetGBM_AR$Dspear>=0.2,na.rm=TRUE),
              sum(evalmetRF_AR$Dspear>=0.2,na.rm=TRUE))
#sum(evalmetGLM_AR$Dspear>=0.2,na.rm=TRUE)/sum(!(is.na(evalmetGLM_AR$Dspear))) #83%
nEval05_AR<-c(sum(evalmetGLM_AR$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGAM_AR$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetGBM_AR$Dspear>=0.5,na.rm=TRUE),
              sum(evalmetRF_AR$Dspear>=0.5,na.rm=TRUE))
#sum(evalmetGLM_AR$Dspear>=0.5,na.rm=TRUE)/sum(!(is.na(evalmetGLM_AR$Dspear))) #18%
nEval04_AR<-c(sum(evalmetGLM_AR$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGAM_AR$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetGBM_AR$Dspear>=0.4,na.rm=TRUE),
              sum(evalmetRF_AR$Dspear>=0.4,na.rm=TRUE))
#sum(evalmetGLM_AR$Dspear>=0.4,na.rm=TRUE)/sum(!(is.na(evalmetGLM_AR$Dspear))) #40%
nEval06_AR<-c(sum(evalmetGLM_AR$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGAM_AR$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetGBM_AR$Dspear>=0.6,na.rm=TRUE),
              sum(evalmetRF_AR$Dspear>=0.6,na.rm=TRUE))
nEval08_AR<-c(sum(evalmetGLM_AR$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGAM_AR$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetGBM_AR$Dspear>=0.8,na.rm=TRUE),
              sum(evalmetRF_AR$Dspear>=0.8,na.rm=TRUE))

UseMod_All_AR<-sum(apply(cbind(evalmetGLM_AR$Dspear>0.2,evalmetGAM_AR$Dspear>0.2,evalmetGBM_AR$Dspear>0.2,evalmetRF_AR$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})==4)
UseMod_One_AR<-sum(apply(cbind(evalmetGLM_AR$Dspear>0.2,evalmetGAM_AR$Dspear>0.2,evalmetGBM_AR$Dspear>0.2,evalmetRF_AR$Dspear>0.2),1,function(X){sum(X,na.rm=TRUE)})!=0)
GoodMod_All_AR<-sum(apply(cbind(evalmetGLM_AR$Dspear>0.4,evalmetGAM_AR$Dspear>0.4,evalmetGBM_AR$Dspear>0.4,evalmetRF_AR$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})==4)
GoodMod_One_AR<-sum(apply(cbind(evalmetGLM_AR$Dspear>0.4,evalmetGAM_AR$Dspear>0.4,evalmetGBM_AR$Dspear>0.4,evalmetRF_AR$Dspear>0.4),1,function(X){sum(X,na.rm=TRUE)})!=0)

Tablenumbers<-data.frame(Number_of=c("ASV","model successfully fitted by at least 1 algorithm","model successfully fitted by all 4 algorithms","model with rho>0.2 for at least 1 algo","model with rho>0.2 for all 4 algo", "models with rho>0.4 for at least 1 algo","model with rho>0.4 for all 4 algos"),
                         Bacteria_250=c(max(nMod_BA),success_One_BA,success_All_BA,UseMod_One_BA,UseMod_All_BA,GoodMod_One_BA,GoodMod_All_BA),
                         Archaea_250=c(max(nMod_AR),success_One_AR,success_All_AR,UseMod_One_AR,UseMod_All_AR,GoodMod_One_AR,GoodMod_All_AR),
                         Fungi_217=c(max(nMod_FU),success_One_FU,success_All_FU,UseMod_One_FU,UseMod_All_FU,GoodMod_One_FU,GoodMod_All_FU),
                         Protist_166=c(max(nMod_PR),success_One_PR,success_All_PR,UseMod_One_PR,UseMod_All_PR,GoodMod_One_PR,GoodMod_All_PR))
write.csv(Tablenumbers,file="figures/PAAB_selection/Table_numberofmodels_Ab.csv")



Tableresult<- data.frame(Bacteria=c(nMod_BA,nFit02_BA,nEval02_BA,nEval05_BA, nEval04_BA,nEval06_BA,nEval08_BA), 
                         Archaea=c(nMod_AR,nFit02_AR,nEval02_AR,nEval05_AR, nEval04_AR,nEval06_AR,nEval08_AR),
                         Fungi=c(nMod_FU,nFit02_FU,nEval02_FU,nEval05_FU, nEval04_FU,nEval06_FU,nEval08_FU),
                         Protist=c(nMod_PR,nFit02_PR,nEval02_PR,nEval05_PR, nEval04_PR,nEval06_PR,nEval08_PR),
                         Model=rep(c("GLM","GAM","GBM","RF"),7),
                         Counted=c(rep("Fitted",4),rep("Fit02",4),rep("rho>0.2",4),rep("rho>0.5",4),rep("rho>0.4",4),rep("rho>0.6",4),rep("rho>0.8",4)))
Tablesresult2<- rbind(c(c(max(nMod_BA),max(nMod_AR),max(nMod_FU),max(nMod_PR),"data","nASV_data")),
                      Tableresult)
Tablesresult2$Bacteria<-as.numeric(Tablesresult2$Bacteria)
Tablesresult2$Archaea<-as.numeric(Tablesresult2$Archaea)
Tablesresult2$Fungi<-as.numeric(Tablesresult2$Fungi)
Tablesresult2$Protist<-as.numeric(Tablesresult2$Protist)
Tablesresult2$Counted<-factor(Tablesresult2$Counted,levels=unique(Tablesresult2$Counted))

Tablesresult3<-Tablesresult2
Tablesresult3[,1:4]<-apply(Tablesresult3[,1:4],2,function(X){X/X[1]})

M_results3 <- melt(Tablesresult3, id= c("Model","Counted"))




M_results<-melt(Tablesresult2[-2,], id= c("Model","Counted"))

levels(M_results$Counted)
library(tidyverse)
library(ggplot2)
library(reshape2) 
pGLM<-ggplot(M_results[M_results$Model%in%c("data","GLM"),], aes(x=Counted, y=value,group=variable,color=variable)) +
  geom_line() + geom_point() + ggtitle("GLM") +scale_color_brewer(palette="Paired")+theme_minimal()
pGAM<-ggplot(M_results[M_results$Model%in%c("data","GAM"),], aes(x=Counted, y=value,group=variable,color=variable))+ 
  geom_line()+geom_point() + ggtitle("GAM") +scale_color_brewer(palette="Paired")+theme_minimal()

pGBM<-ggplot(M_results[M_results$Model%in%c("data","GBM"),], aes(x=Counted, y=value,group=variable,color=variable))+
  geom_line()+geom_point() + ggtitle("GBM") +scale_color_brewer(palette="Paired")+theme_minimal()
pRF<-ggplot(M_results[M_results$Model%in%c("data","RF"),], aes(x=Counted, y=value,group=variable,color=variable))+ 
  geom_line()+geom_point() + ggtitle("RF") +scale_color_brewer(palette="Paired")+theme_minimal()

library(grid)
library(gridExtra)

grid.arrange(pGLM,pGAM,pGBM,pRF)  


#same %
M_results3AB <- melt(Tablesresult3[Tablesresult3$Model!="data"&Tablesresult3$Counted!="rho>0.5"&Tablesresult3$Counted!="Fit02",], id= c("Model","Counted"))

# pGLM3<-ggplot(M_results3AB[M_results3AB$Model%in%c("data","GLM"),], aes(x=Counted, y=value,group=variable,color=variable)) +
#   geom_line() + geom_point() + ggtitle("GLM") +scale_color_brewer(palette="Paired")+theme_minimal()
# pGAM3<-ggplot(M_results3AB[M_results3AB$Model%in%c("data","GAM"),], aes(x=Counted, y=value,group=variable,color=variable))+ 
#   geom_line()+geom_point() + ggtitle("GAM")+scale_color_brewer(palette="Paired")+theme_minimal()
# 
# pGBM3<-ggplot(M_results3AB[M_results3AB$Model%in%c("data","GBM"),], aes(x=Counted, y=value,group=variable,color=variable))+
#   geom_line()+geom_point() + ggtitle("GBM")+scale_color_brewer(palette="Paired")+theme_minimal()
# pRF3<-ggplot(M_results3AB[M_results3AB$Model%in%c("data","RF"),], aes(x=Counted, y=value,group=variable,color=variable))+ 
#   geom_line()+geom_point() + ggtitle("RF")+scale_color_brewer(palette="Paired")+theme_minimal()
# grid.arrange(pGLM3,pGAM3,pGBM3,pRF3)  
# 
# 
# pdf(file="figures/testmodels/Proportion_successful_models_AB.pdf")
# grid.arrange(pGLM,pGAM,pGBM,pRF) 
# grid.arrange(pGLM3,pGAM3,pGBM3,pRF3)
# dev.off()

Fig_AB_modquality_per_dataset_4group<-ggplot(M_results3AB, aes(x=Counted, y=value,group=Model)) +
  geom_point(aes(shape=variable, color=variable),size=2) +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_color_manual(
    values = c('#FFD700', '#5CB800','#6B4C62','#457EB0'),labels = c("Bacteria", "Archaea", "Fungi", "Protista")
  ) +
  facet_wrap(~Model) +
  scale_x_discrete(labels=c("rho>0.2" = expression(rho~">0.2"), "rho>0.4" = expression(rho~">0.4"),"rho>0.6" = expression(rho~">0.6"),"rho>0.8" = expression(rho~">0.8"))) +
  geom_line(aes(x=Counted, y=value, group=variable, color= variable ))  +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(angle = 20, margin = margin(t = 5),hjust=1, size=10),
                     axis.text.y = element_text(size=10),
                     legend.position = c(1,1),
                     legend.justification = c(1, 1),
                     legend.text = element_text(size=8),
                     legend.key.size = unit(0.4,"cm"),
                     legend.title = element_blank(),
                     legend.background = element_rect(fill = NA, colour = NA),
                     strip.text.x = element_text(size = 15),
                     plot.margin = unit(c(0.2,0.2,0,0.2),"cm")) + 
  labs(x="",y="Proportion of models") +
  guides(shape = "none",
         color = guide_legend(
           override.aes=list(shape = c(15,16,17,18))))
pdf(file="figures/PAAB_selection/testmodels/Proportion_successful_modelsAB3.pdf")
plot(Fig_AB_modquality_per_dataset_4group)
dev.off()

png(file=paste0("figures/PAAB_selection/testmodels/Proportion_successful_modelsAB3.png"),res=300,width=1961,height=1500)
plot(Fig_AB_modquality_per_dataset_4group)
dev.off()

png(file=paste0("figures/PAAB_selection/testmodels/Proportion_successful_modelsPAAB.png"),res=300,width=1961,height=2000)
grid.arrange(Fig_PA_modquality_per_dataset_4group,Fig_AB_modquality_per_dataset_4group,nrow=2)
dev.off()
