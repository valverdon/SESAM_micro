
########### WORK IN PROGRESS

###################################################################################################
# Full pipe doing fit of models. Called from cluster using arguments    :                            # 
#   GtoM for microbial group ; 
#   Modeltype for Model algorythm used ;                                      #
#   PAAB for modelling binary response (P-A data) or count response (Richness, rel Abundance)       #                                    #
#     e.g. sbatch 01_Fit.txt PR GAM PA   doing Fit of PRotists modeled by GAM algo on Presence-Absence data   
# First argument is the batch job (.txt) that tells to run that R script, 
# arg 2 : group to model; tells the script in which folder to read data and find results (mine are called PR FU and BA)
# other arguments are passed to R as options to activate/inactivate parts of the script
# arg 3 : Modelling algorithm to use : either GAM GLM RF LGBM or All; 
# activates the corresponding parts of the code
# arg 4 : either PA AB or both, activates binary or count modelling
# 
# directory from which the bash code is read should be the "Both" folder, 
# with correct arborescence to reach files to read and files to write
# 
# The idea of the script is to run a set of ASV determined by the "arrayID" so that the cluster 
# will run the following script several time in parallel on the cluster (using socketting : each node in the cluster 
# works on one array ID, and each array ID correspond to a subset of all ASV I want to model 
###################################################################################################
.libPaths("/work/FAC/FBM/DEE/aguisan/sometalp/rlib/4.2") #where to find R libraries on the cluster

#load needed packages
library(plyr) #parallel
library(tidyverse) #data manip
library(doParallel) #parallel
#for evalmetrics # loaded in evalmetrics
# library(modEvA)
# library(ecospat)
# library(ROCR)

# set modelling options

#seed for reproducibility
set.seed(seed=08022022)#date of the day I added a seed
#get the modelling options directly from SLURM
args <- unlist(commandArgs(trailingOnly = TRUE) )
# args <- c(1,"BA","All","PA","test")
# args <- c(200,"PR","GLM","AB")

arrayID= as.numeric(args[1])  #will determine which subset of ASVs to model
GtoM <- args[2]               #will determine which group to model 
Modeltype <- args[3]          # activates the corresponding parts of the code
PAAB <- ifelse(args[4]=="PA" | args[4]=="AB", args[4], "PA") # activates the corresponding parts of the code, PA by default
test <- ifelse("test"%in%args,TRUE,FALSE) #activate the testing mode, not test by default

#load Environmental data
ENVdata<-as.data.frame(readRDS(paste0("data/ENVdata_",GtoM,".Rds")))
if (PAAB=="PA"){
  OTUdata<-readRDS(paste0("data/OTUdata_",GtoM,".Rds"))
}
if (PAAB=="AB"){
  OTUdata<-readRDS(paste0("data/OTUdata_",GtoM,"_AB.Rds"))
}

#end of data preparation, following are all codes that can be switched on and off by options
#later some commented things don't make sense in english.
#I use them to dig into the code for debugging, just ignore them
# source("code/Fit_GAM.R")
source("code/Fit_GLM.R")
source("code/Eval_GLM.R")
source("code/Fit_RF.R")
source("code/Eval_RF.R")
source("code/Fit_GBM.R")
source("code/Eval_GBM.R")

if(test){
  OTUstoRun <-sample(ncol(OTUdata),10)
  if(Modeltype=="All"|Modeltype=="GLM"){
      time1<-Sys.time()
      Modlist<-Fit_GLM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,test=TRUE)
      time2<-Sys.time()
      print(time2-time1)
      Eval_GLM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=TRUE)
      # Eval_GLM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,Model_list=Modlist,test=TRUE, validation.method="random.split-sample")
      # Eval_GLM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,Model_list=Modlist,test=TRUE, validation.method="LOO")
      # Eval_GLM_Varsel(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=TRUE)
      time3<-Sys.time()
      print(time3-time1)
      print(time3-time2)
    closeAllConnections()
  }
  if(Modeltype=="All"|Modeltype=="RF"){
      time1<-Sys.time()
      Modlist<-Fit_RF(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,test=TRUE)
      time2<-Sys.time()
      print(time2-time1)
      Eval_RF(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=TRUE)
      # Eval_RF(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,Model_list=Modlist,test=TRUE, validation.method="random.split-sample")
      # Eval_RF(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,Model_list=Modlist,test=TRUE, validation.method="random.split-sample",DataSplit="LOO")
      # 
      time3<-Sys.time()
      print(time3-time1)
      print(time3-time2)
      closeAllConnections()
  }
  if(Modeltype=="All"|Modeltype=="GBM"){
    time1<-Sys.time()
      Modlist<-Fit_GBM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,test=TRUE)
      time2<-Sys.time()
      print(time2-time1)
       Eval_GBM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=TRUE)
       time3<-Sys.time()
       print(time3-time1)
       print(time3-time2)
       closeAllConnections()
  }
  # if(Modeltype=="All"|Modeltype=="GAM"){
  #   if(test){
  #     Modlist<-Fit_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,test=TRUE)
  #     Eval_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=TRUE)
  #   }else {
  #     Modlist<-Fit_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,test=FALSE)
  #     Eval_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist)
  #   }
  #     closeAllConnections()
  # }
    
} 
if(!(test)){
  endloop <- ceiling(ncol(OTUdata)/250)
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop} 
  if (endthisloop<=0) { q("no")} #if the array don't contain anything (finished), quit R
  #endthisloop=5
  OTUstoRun <-((arrayID-1)*endloop+1):((arrayID-1)*endloop+endthisloop)
  if(Modeltype=="All"|Modeltype=="GLM"){
      Modlist<-Fit_GLM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,test=FALSE)
      Eval_GLM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=FALSE)
    closeAllConnections()
  }
  if(Modeltype=="All"|Modeltype=="RF"){
      Modlist<-Fit_RF(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,test=FALSE)
      Eval_RF(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=FALSE)
    closeAllConnections()
  }
  if(Modeltype=="All"|Modeltype=="GBM"){
      Modlist<-Fit_GBM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,test=FALSE)
      Eval_GBM(PAAB=PAAB,arrayID=arrayID,OTUdata=OTUdata,OTUstoRun=OTUstoRun,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=FALSE)
    closeAllConnections()
  }
  # if(Modeltype=="All"|Modeltype=="GAM"){
  #   if(test){
  #     Modlist<-Fit_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,test=TRUE)
  #     Eval_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist,test=TRUE)
  #   }else {
  #     Modlist<-Fit_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,test=FALSE)
  #     Eval_GAM(PAAB=PAAB,arrayID=arrayID,GtoM=GtoM,ENVdata=ENVdata,NbRunEval=100,Model_list=Modlist)
  #   }
  #     closeAllConnections()
  # }
  
}

q("no")

######################################################################################################
#Test validation method : Same models fitted, different evaluation methods
load(paste0("PA/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp41.Rda"))
apply(Eval_met_mat[,c("TSS","Kappa","Sensitivity","Specificity","AUC")],2,mean,na.rm=TRUE)

#r.split.sampling PRGLMPA full array 41
# TSS       Kappa Sensitivity Specificity         AUC 
# 0.3266776   0.3022550   0.6711517   0.6836433   0.6804969 

#Bootstrap PRGLMPA full array 41
# TSS       Kappa Sensitivity Specificity         AUC 
# 0.3323616   0.2988006   0.7140921   0.6432774   0.6846073 

#leave one out PRGLMPA full array 41
# TSS       Kappa Sensitivity Specificity         AUC 
# 0.2939630   0.2656627   0.6450164   0.6489465   0.6503283
boxplot(Eval_met_mat[,"TSS"])

     
load(paste0("PA/",GtoM,"/Outputs/RF/Eval_met/Eval_met_temp41.Rda"))
apply(Eval_met_mat[,c("TSS","Kappa","Sensitivity","Specificity","AUC")],2,function(X){mean(X,na.rm=TRUE)})
#Bootstrap PRRFPA full array 41
# TSS       Kappa Sensitivity Specificity         AUC 
# 0.6393342   0.2412015   0.6273809   0.6648474   0.6500535 

#r.split.sampling PRRFPA full array 41
apply(Eval_met_mat[,c("TSS","Kappa","Sensitivity","Specificity","AUC")],2,mean)
# TSS       Kappa Sensitivity Specificity         AUC 
# 0.6316244   0.2340237   0.6180002   0.6583439   0.6462839 

#leave one out PRRFPA full array 41
# TSS       Kappa Sensitivity Specificity         AUC 
# 0.6332075   0.2397675   0.6369084   0.6444118   0.6455337 
boxplot(Eval_met_mat[,"TSS"])

