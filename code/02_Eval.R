###################################################################################################
# Full pipe for evaluation of models. from cluster using arguments                                #                                   #
# e.g. sbatch 01_Fit.txt PR GAM PA   for protists modeled by GAM algo    
# First argument is the batch job (.txt) that tells to run that R script, 
# arg 2 : group to model; tells the script in which folder to read data and find results (mine are called PR FU and BA)
# other arguments are passed to R as options to activate/inactivate parts of the script
# arg 3 : Modelling algorithm to use : either GAM GLM RF LGBM or All; 
# activates the corresponding parts of the code
# arg 4 : either PA AB or both, activates binary (P-A data) or count modelling(Richness, rel Abundance)
# 
# directory should be "Both" folder, 
# with correct arborescence to reach files to read and files to write
# 
# The idea of the script is to run a set of ASV determined by the "arrayID" so that the cluster 
# will run that script in parallel on the cluster (using socketting : each node in the cluster 
# works on one array ID, and each array ID correspond to a subset of all ASV of the modelled 
# microbial group)
###################################################################################################
.libPaths("/work/FAC/FBM/DEE/aguisan/default/rlib/4.0")

library(tidyverse)
library(parallel)
library(doParallel)

#for gam
library(mgcv)

#for evalmetrics
library(modEvA)
library(ecospat)
library(ROCR)

#for GLM
library(data.table)
library(stringi)
library(glmnet)
library(MASS)
library(plyr)

#for RF
library(randomForest)
library(RRF)
# 
#for light GBM
library(lightgbm)

# arrayID=145
# GtoM="FU"
# Modeltype <- "GAM"
# PAAB <- "AB"
#get the modelling options directly from SLURM
args <- commandArgs(trailingOnly = TRUE) 
arrayID= as.numeric(args[[1]])  #will determine which subset of ASVs to model
GtoM <- args[[2]]               #will determine which group to model 
Modeltype <- args[[3]]          # activates the corresponding parts of the code
PAAB <- args[[4]]               # activates the corresponding parts of the code


#functions I wrote that will be used
source("code/Evalmetrics.R")#Function that compute some metrics of model evqluqtion
source("code/create_DataSplitTable.R") #function to split dataset for the CV procedure.

#quiet function by Hadley Wickham. Will be useful to prevent some algos to print to many things in the console
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

#parallelization settings
nthreads <- detectCores()/2
registerDoParallel(cores = nthreads)

#seed for reproducibility
set.seed(seed=08022022) #date of the day I added a seed

load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))

NbRunEval <- 100
binvar <- c()
contvar <- c()
for (i in names(ENVdata)){#i =names(ENVdata) [31]
  if(length(unique(ENVdata[,i]))>2){
    contvar <- c(contvar,i)
  } else {
    binvar <- c(binvar,i)
  }
}

#######################################################################################################################################
#######################################################################################################################################
########################################################   P - A   ####################################################################
#######################################################################################################################################
#######################################################################################################################################

if(PAAB=="PA"|PAAB=="both"){
  load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))
  
  endloop <- ceiling(ncol(OTUdata)/250) #Number of ASV to run
  #endthisloop=10
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop}#correct number of ASV for last loop
  if (endthisloop<=0) { q("no")} #quit r if already finished

  ####################################################################################################################################
  ###########################################################  GAM  ###################################################################
  ####################################################################################################################################
  # Time1=Sys.time()
  if(Modeltype=="All"|Modeltype=="GAM"){
    #load models fitted in 01_Fit.
    load(paste0("PA/",GtoM,"/Outputs/GAM/Models/Models_temp",arrayID,".Rda"))
    
    #prepare lists of results. Predefine length of list?
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      Mgam <- Model_list[[paste0("OTU",OTUtoRun)]] #select the model of the selected ASV
      if(all(is.na(Mgam))){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method : Cross validation  
        #I use bootstrap .632+  Efron & Tibshirani 1997 
        #But split sampling is possible in the create.datasplittable function
        
        # set.seed(143)
        dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
        
        pred_GAM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))
        
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)

          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          
          formCV <- as.formula(paste0("OTUtrain~ ", paste(Mgam$formula)[3]))
          MgamEval <- tryCatch(
            { gam(formCV, data=as.data.frame(ENVtrain), family="binomial", weights=weights, select=T,nthreads=nthreads) #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          if(is.na(MgamEval[1])){
            next
          }else{
            predCV <- predict.gam(MgamEval, newdata = ENVeval, type = "response")
            pred_GAM[rownames(ENVeval), paste0("pred",f)] <- predCV  
            tryCatch(
              { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%names(predCV),OTUtoRun],pred=as.vector(predCV),PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
              }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                message(paste("boot will be ignored, original error :"))
                message(cond)
                return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
              }
            )
          }
        }#end boot    
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GAM)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }#end else
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAM/Eval_met"))){
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GAM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAM/Eval_met"))
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GAM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAM/Pred_data"))){
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GAM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAM/Pred_data"))
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GAM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GAM/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GAM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GAM/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GAM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }#end GAM
  
  
  ####################################################################################################################################
  ###########################################################  GLM  ###################################################################
  ####################################################################################################################################
  #Time1=Sys.time()
  # Sys.time()-Time1
  if(Modeltype=="All"|Modeltype=="GLM"){
    load(paste0("PA/",GtoM,"/Outputs/GLM/Models/Models_temp",arrayID,".Rda"))
    load(paste0("PA/",GtoM,"/Outputs/GLM/VarSel/VarSel_temp",arrayID,".Rda"))
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    #for each species
    for (i in 1:endthisloop){ #i=17 ;endthisloop=5
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-8753
      Mglm <- Model_list[[paste0("OTU",OTUtoRun)]]#plot(Mglm)
      preselected<-eval(parse(text=paste0("VarSel$preselection$OTU",OTUtoRun)))
      if(is.na(Mglm[1])){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method : Cross validation 80% - 20%   
        #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
        
        # set.seed(143)
        dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")

        pred_GLM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "Specificity", "PPV", "NPV", "Jaccar", "TSS", "Kappa", "SEDI")))
        
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)

          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]

          formT<-as.formula(paste0("OTUtrain~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          formE<-as.formula(paste0("OTUeval~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          
          ModMatT <- model.matrix(formT,data = ENVtrain)
          ModMatE <- model.matrix(formE,data = ENVeval)
          MglmEval <- tryCatch(
            { cv.glmnet(ModMatT, OTUtrain, family="binomial", weights=weights, alpha=1, parallel = TRUE, n.cores = nthreads)
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          
          
          if(is.na(MglmEval[1])){
            Eval_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA))
            Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
            next
          } else{
            # Extract results
            lambda<- list(onese=MglmEval$lambda.1se,min=MglmEval$lambda.min)
            if(MglmEval$lambda[1]!=MglmEval$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
              fitted.val <- predict(MglmEval,newx=ModMatE, s=MglmEval$lambda.1se, type="response")
            } else{#take lambda that minimize error
              fitted.val <- predict(MglmEval,newx=ModMatE, s=MglmEval$lambda.min, type="response")
            }
            if(empty(fitted.val)){
              Eval_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA))
              Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
              next
            } else{
              colnames(fitted.val) <- "fit"

              pred_GLM[rownames(ENVeval), paste0("pred",f)] <- fitted.val
              tryCatch(
                { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%rownames(fitted.val),OTUtoRun],pred=fitted.val,PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
                }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                  message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                  message(paste("boot will be ignored"))
                  return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
                }
              )
            }
          }#end what to do if Modeleval exist
        }#end that split
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GLM)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }# end if model exist for that ASV
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/Eval_met"))){
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/Eval_met"))
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/Pred_data"))){
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GLM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/Pred_data"))
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GLM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GLM/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GLM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GLM/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GLM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }
  
  ####################################################################################################################################
  ###########################################################  RF  ###################################################################
  ####################################################################################################################################
  if(Modeltype=="All"|Modeltype=="RF"){
    load(paste0("PA/",GtoM,"/Outputs/RF/Models/Models_temp",arrayID,".Rda"))
    load(paste0("PA/",GtoM,"/Outputs/RF/VarSel/VarSel_temp",arrayID,".Rda"))
    
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      MRF <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(is.na(MRF[1])){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method : Cross validation 80% - 20%   
        #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
        
        # set.seed(143)
        dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
        
        pred_RF <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "Specificity", "PPV", "NPV", "Jaccar", "TSS", "Kappa", "SEDI")))
        preselected <- VarSel$preselection[[paste0("OTU",OTUtoRun)]]
        
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})

          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)
          
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]

          formCV <- as.formula(paste0("factor(OTUtrain)~ ", paste(preselected,collapse= " + ")))
          #################################
          MRFEval <- tryCatch(
            { RRF(formCV, data=as.data.frame(ENVtrain), ntree=100, classwt=c("FALSE"=unique(weights[OTUtrain==FALSE]),"TRUE"=unique(weights[OTUtrain==TRUE])))
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          if(!(is.na(MRFEval[1]))){
            predCVprob <- predict(MRFEval, newdata = ENVeval, type = "prob")[,2]
            predCVresp <- predict(MRFEval, newdata = ENVeval, type = "response")
            pred_RF[rownames(ENVeval), paste0("pred",f)] <- c(proba=predCVprob)
            tryCatch(
              { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%names(predCVprob),OTUtoRun],pred=as.vector(predCVprob),PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
              }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                message(paste("boot will be ignored"))
                return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
              }
            )
          }
        }
        
        
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_RF)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }#end else
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/RF/Eval_met"))){
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/RF/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/RF/Eval_met"))
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/RF/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/RF/Pred_data"))){
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/RF/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/RF/Pred_data"))
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/RF/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/RF/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/RF/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/RF/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/RF/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }#end RF
  
  
  ####################################################################################################################################
  ###########################################################  GBM  ###################################################################
  ####################################################################################################################################
  if(Modeltype=="All"|Modeltype=="GBM"){
    load(paste0("PA/",GtoM,"/Outputs/GBM/Models/Models_temp",arrayID,".Rda"))
    load(paste0("PA/",GtoM,"/Outputs/GBM/VarSel/VarSel_temp",arrayID,".Rda"))
    
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      lgb_model <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(!(exists(lgb_model,mode="environment"))){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method : Cross validation 80% - 20%   
        #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
        
        # set.seed(143)
        dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
        
        pred_GBM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))
        
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
          data_train <-  data.table(ENVtrain,obsval=OTUtrain)
          
          data_valid <-  data.table(ENVeval)
          
          # MODEL #################################
          
          dtrain_data <- as.matrix(data_train[,names(ENVdata),with=FALSE])
          dvalid_data <- as.matrix(data_valid[,names(ENVdata),with=FALSE])
          dtrain_label <- as.matrix(data_train[,.(obsval)])
          
          dtrain <- lgb.Dataset(
            data = dtrain_data,
            label = dtrain_label,
            categorical_feature = which(names(ENVdata)%in%binvar),
            weight=weights)
          
          lgb_model_cv <- lgb.cv( #using biomod default values except ntree = nrounds = 1000
            params = list(
              objective = "binary",
              boosting="gbdt",
              learning_rate = 0.001,
              max.depth = 7,
              bagging_fraction =0.5,
              feature_pre_filter=FALSE),
            verbose = -1L,
            data = dtrain,
            nrounds = 1000,
            early_stopping_rounds = 50,
            eval_freq = 20,
            eval = "auc",
            nfold = 3,
            stratified = TRUE)
          nrounds= lgb_model_cv$best_iter
          
          lgb_model <- lgb.train(
            data=dtrain, 
            verbose=-1, #skip warning messages
            nrounds = nrounds,
            params= list(
              objective = "binary",
              # poisson_max_delta_step = poisson_max_delta_step,
              boosting="gbdt"))
          # if(is.na(lgb_model)){
          #   ranking_list[[paste0("OTU",OTUtoRun)]]<- NA
          #   Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, dev.ratio = NA)
          #   Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
          #   Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
          #   Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
          #   next
          # } else{
          pred_valid <- predict(lgb_model,dvalid_data)
          tryCatch(
            { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%rownames(ENVeval),OTUtoRun],pred=as.vector(pred_valid),PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
              message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
              message(paste("boot will be ignored, original error :"))
              message(cond)
              return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
            }
          )
          
          Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_valid,pred_fit)
          Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
          Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
        }#end else
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met"))){
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met"))
      save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GBM/Pred_data"))){
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GBM/Pred_data"))
      save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }#end GBM
}#end PA












#Eval
NbRunEval <- 100


##############################################################################################################
#### P2 : Model evaluation ##################################################################

#Method : Cross validation 80% - 20%   
#Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 

# set.seed(143)
dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")

pred_GBM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))

for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
  # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
  # ENVtrain <- ENVdata[dataSplits[,f], ]
  # Offset <- log(TotSeqSum[dataSplits[,f]])
  OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
  ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
  ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
  weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
  data_train <-  data.table(ENVtrain,obsval=OTUtrain)
  
  data_valid <-  data.table(ENVeval)
  
  # MODEL #################################
  
  dtrain_data <- as.matrix(data_train[,names(ENVdata),with=FALSE])
  dvalid_data <- as.matrix(data_valid[,names(ENVdata),with=FALSE])
  dtrain_label <- as.matrix(data_curr[,.(obsval)])
  
  dtrain <- lgb.Dataset(
    data = dtrain_data,
    label = dtrain_label,
    categorical_feature = which(names(ENVdata)%in%binvar),
    weight=weights)
  
  # lgb_model <- run_model(dtrain = dtrain, learning_rate = 0.3, num_iterations = 100,
  # min_sum_hessian = 0, poisson_max_delta_step = 0.6 )
  # lgb_model_cv <- lgb.cv(
  #   params = list(
  #     objective = "binary",
  #     boosting="gbdt",
  #     learning_rate = 0.01, 
  #     num_leaves = 25,
  #     num_threads = 2), 
  #   verbose = -1L,
  #   data = dtrain,
  #   nrounds = 1000, 
  #   early_stopping_rounds = 50,
  #   eval_freq = 20, 
  #   eval = "auc", 
  #   nfold = 5, 
  #   stratified = TRUE)
  # best_iter <- lgb_model_cv$best_iter
  # best_iter 
  #did some test number of trees (iterations/rounds) do not need to be that high (100 is already more than enough for most ASV)
  
  lgb_model <- lgb.train(
    data=dtrain, 
    verbose=-1, #skip warning messages
    nrounds = nrounds,
    params= list(
      objective = "binary",
      # poisson_max_delta_step = poisson_max_delta_step,
      boosting="gbdt"))
  # if(is.na(lgb_model)){
  #   ranking_list[[paste0("OTU",OTUtoRun)]]<- NA
  #   Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, dev.ratio = NA)
  #   Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
  #   Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
  #   Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
  #   next
  # } else{
  pred_valid <- predict(lgb_model,dvalid_data)
  tryCatch(
    { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%rownames(ENVeval),OTUtoRun],pred=as.vector(pred_valid),PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
    }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
      message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
      message(paste("boot will be ignored, original error :"))
      message(cond)
      return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
    }
  )
  
  Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_valid,pred_fit)
  Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
  Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
}#end else

}
}



#######################################################################################################################################
#######################################################################################################################################
########################################################  ABUNDANCE ###################################################################
#######################################################################################################################################
#######################################################################################################################################

if(PAAB=="AB"|PAAB=="both"){
  load(paste0("Abundance/",GtoM,"/data/OTUdata.Rda"))
  load(paste0("Abundance/",GtoM,"/data/dataTotSeqSum.Rda"))
  
  

  endloop <- ceiling(ncol(OTUdata)/250)
  #endthisloop=10
  if ((arrayID-1)*endloop+endloop > ncol(OTUdata)){endthisloop <- ncol(OTUdata)-((arrayID-1)*endloop)} else {endthisloop <- endloop}#207 on last loop
  if (endthisloop<=0) { q("no")} #already finished

  ####################################################################################################################################
  ###########################################################  GAM  ###################################################################
  ####################################################################################################################################
  # Time1=Sys.time()
  if(Modeltype=="All"|Modeltype=="GAM"){
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    load(paste0("Abundance/",GtoM,"/Outputs/GAM/Models/Models_temp",arrayID,".Rda"))
    
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      #From 04_
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-29996
      Mgam <- Model_list[[paste0("OTU",OTUtoRun)]]
      ##############################################################################################################
      #### P2 : Model evaluation ##################################################################
      
      #Method  
      #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 

      dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
      
      pred_GAM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
      Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("R2", "D2", "MAE", "MAEs", "RMSE", "RMSEs", "Dspear", "Dpear", "Pdispersion","R2_scaled","MAE_scaled", "RMSE_scaled", "Dspear_scaled", "Dpear_scaled", "Pdisp_scaled")))
      if (is.na(Mgam[1])){
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GAM)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- rep(NA, 15)
        next
      } else {
        for (f in 1:ncol(dataSplits)) { #f=2 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          TotSeqSumtrain <- rep(TotSeqSum,times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
  
          formCV <- as.formula(paste0("OTUtrain~ ", gsub("TotSeqSum","TotSeqSumtrain",paste(summary(Mgam)$formula)[3])))
          MgamEval <- tryCatch(
            { gam(formCV, data=as.data.frame(ENVtrain), family="poisson", select=T,nthreads=nthreads, weights = weights) #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          if(is.na(MgamEval[1])){
            next
          }else{
            predCV <- predict.gam(MgamEval, newdata=data.frame(ENVeval, TotSeqSumtrain = TotSeqSum[dataSplits[ ,f]==FALSE]), type="response")
            pred_GAM[rownames(ENVeval), paste0("pred",f)] <- predCV  
            tryCatch(
              { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%names(predCV),OTUtoRun],pred=predCV,PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
              }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                message(paste("boot will be ignored, original error :"))
                message(cond)
                return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
              }
            )
          }
        }
      }
      
      pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
      
      Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GAM)
      Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
      Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met"))){
      save(Eval_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met"))
      save(Eval_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAM/Pred_data"))){
      save(Pred_data_list, file=paste0("Abundance/",GtoM,"/Outputs/GAM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAM/Pred_data"))
      save(Pred_data_list, file=paste0("Abundance/",GtoM,"/Outputs/GAM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("Abundance/",GtoM,"/Outputs/GAM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }#end GAM
  
  ####################################################################################################################################
  ###########################################################  GLM  ###################################################################
  ####################################################################################################################################
  # Time1=Sys.time()
  # Sys.time()-Time1
  if(Modeltype=="All"|Modeltype=="GLM"){
    load(paste0("Abundance/",GtoM,"/Outputs/GLM/Models/Models_temp",arrayID,".Rda"))
    load(paste0("Abundance/",GtoM,"/Outputs/GLM/VarSel/VarSel_temp",arrayID,".Rda"))
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    for (i in 1:endthisloop){ #i=1 ; endthisloop=10
      # for (i in 1:2){ print(i)
      # print(i)
      #From 04_
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-8753
      Mglm <- Model_list[[paste0("OTU",OTUtoRun)]]
      preselected<-eval(parse(text=paste0("VarSel$preselection$OTU",OTUtoRun)))
      
      dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
      
      pred_GLM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
      Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("R2", "D2", "MAE", "MAEs", "RMSE", "RMSEs", "Dspear", "Dpear", "Pdispersion","R2_scaled","MAE_scaled", "RMSE_scaled", "Dspear_scaled", "Dpear_scaled", "Pdisp_scaled")))
      
      if(!(is.na(Mglm[1]))){#dont bother evaluate if no model
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method  
        #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
        

        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          TotSeqSumtrain <- rep(TotSeqSum,times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
          TotSeqSumeval <- TotSeqSum[dataSplits[,f]==FALSE]
          
          formT<-as.formula(paste0("OTUtrain~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          formE<-as.formula(paste0("OTUeval~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
          
          ModMatT <- model.matrix(formT,data = ENVtrain)
          ModMatE <- model.matrix(formE,data = ENVeval)

          MglmEval <- tryCatch(
            { cv.glmnet(ModMatT, OTUtrain, family="poisson", weights = weights, alpha=1, parallel = TRUE, n.cores = nthreads, offset = log(TotSeqSumtrain))
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )
          # modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          if(is.na(MglmEval[1])){
            Eval_met_list[[paste0("OTU",OTUtoRun)]]  <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
            Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
            next
          } else{
            # Extract results
            lambda<- list(onese=MglmEval$lambda.1se,min=MglmEval$lambda.min)
            if(MglmEval$lambda[1]!=MglmEval$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
              fitted.val <- predict(MglmEval,newx=ModMatE, s=MglmEval$lambda.1se, type="response", newoffset = log(TotSeqSumeval))
            } else{#take lambda that minimize error
              fitted.val <- predict(MglmEval,newx=ModMatE, s=MglmEval$lambda.min, type="response", newoffset = log(TotSeqSumeval))
            }
            if(empty(fitted.val)){
              Eval_met_list[[paste0("OTU",OTUtoRun)]]  <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
              Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
              next
            } else{
              colnames(fitted.val) <- "fit"
              
              pred_GLM[rownames(ENVeval), paste0("pred",f)] <- fitted.val
              tryCatch(
                { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%rownames(fitted.val),OTUtoRun],pred=as.numeric(fitted.val),PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
                }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                  message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                  message(paste("boot will be ignored"))
                  return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
                }
              )
            }
          }#end what to do if Modeleval exist
        }#end that split
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GLM)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met"))){
      save(Eval_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met"))
      save(Eval_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GLM/Pred_data"))){
      save(Pred_data_list, file=paste0("Abundance/",GtoM,"/Outputs/GLM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GLM/Pred_data"))
      save(Pred_data_list, file=paste0("Abundance/",GtoM,"/Outputs/GLM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } 
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("Abundance/",GtoM,"/Outputs/GLM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }#end GLM
  
  ####################################################################################################################################
  ###########################################################  RF  ###################################################################
  ####################################################################################################################################
  if(Modeltype=="All"|Modeltype=="RF"){
    load(paste0("Abundance/",GtoM,"/Outputs/RF/Models/Models_temp",arrayID,".Rda"))
    load(paste0("Abundance/",GtoM,"/Outputs/RF/VarSel/VarSel_temp",arrayID,".Rda"))
    
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      MRF <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(is.na(MRF[1])){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
        
      } else{
        
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method : Cross validation 80% - 20%   
        #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
        
        # set.seed(143)
        dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
        
        pred_RF <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=9, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("R2", "D2", "MAE", "MAEs", "RMSE", "RMSEs", "Dspear", "Dpear", "Pdispersion")))
        preselected <- VarSel$preselection[[paste0("OTU",OTUtoRun)]]
        
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)
          
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          
          formCV <- as.formula(paste0("OTUtrain~ ", paste(preselected,collapse= " + ")))
          #################################
          MRFEval <- tryCatch(
            { RRF(formCV, data=as.data.frame(ENVtrain), ntree=100, weights = weights)
            }, error=function(cond){
              message(paste0("error with boot ",f," of OTU ", OTUtoRun))
              message(cond)
              message(paste("boot will be ignored"))
              return(NA)
            }
          )# modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
          #                data=ENVtrain, family="poisson")
          # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
          #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
          
          predCV <- predict(MRFEval, newdata = ENVeval)
          pred_RF[rownames(ENVeval), paste0("pred",f)] <- c(proba=predCV)
          tryCatch(
            { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%names(predCV),OTUtoRun],pred=predCV,PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
              message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
              message(paste("boot will be ignored"))
              return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
            }
          )
        }
        
        
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_RF)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }#end else
    }#end that ASV
    
    ##############################################################################################################
    ####  Save results      ##################################################################
    # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
    # Model evaluation metrics
    Eval_met_mat <- do.call("rbind", Eval_met_list)
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met"))){
      save(Eval_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met"))
      save(Eval_met_mat, file=paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met/Eval_met_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/RF/Pred_data"))){
      save(Pred_data_list, file=paste0("Abundance/",GtoM,"/Outputs/RF/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/RF/Pred_data"))
      save(Pred_data_list, file=paste0("Abundance/",GtoM,"/Outputs/RF/Pred_data/Pred_data_temp", arrayID, ".Rda"))
    }
    if (file.exists(paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met_allboot/"))){
      save(Eval_values_list, file=paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    } else {
      dir.create(paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met_allboot"))
      save(Eval_values_list, file=paste0("Abundance/",GtoM,"/Outputs/RF/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
    }
  }
  
  
  
  ####################################################################################################################################
  ###########################################################  GBM  ###################################################################
  ####################################################################################################################################
  if(Modeltype=="All"|Modeltype=="GBM"){
    load(paste0("AB/",GtoM,"/Outputs/GBM/Models/Models_temp",arrayID,".Rda"))
    load(paste0("AB/",GtoM,"/Outputs/GBM/VarSel/VarSel_temp",arrayID,".Rda"))
    
    Eval_met_list <- list()
    Pred_data_list <- list()
    Eval_values_list <-list()
    
    #for each species
    for (i in 1:endthisloop){ #i=3
      # for (i in 1:2){ print(i)
      # print(i)
      OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  #OTUtoRun<-14961
      lgb_model <- Model_list[[paste0("OTU",OTUtoRun)]]
      
      if(!(exists(lgb_model,mode="environment"))){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        #Method : Cross validation 80% - 20%   
        #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
        
        # set.seed(143)
        dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap")
        
        pred_GBM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))
        
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),1)})
          data_train <-  data.table(ENVtrain,obsval=OTUtrain)
          
          data_valid <-  data.table(ENVeval)
          
          # MODEL #################################
          
          dtrain_data <- as.matrix(data_train[,names(ENVdata),with=FALSE])
          dvalid_data <- as.matrix(data_valid[,names(ENVdata),with=FALSE])
          dtrain_label <- as.matrix(data_curr[,.(obsval)])
          
          dtrain <- lgb.Dataset(
            data = dtrain_data,
            label = dtrain_label,
            categorical_feature = which(names(ENVdata)%in%binvar),
            weight=weights)
          
          # lgb_model <- run_model(dtrain = dtrain, learning_rate = 0.3, num_iterations = 100,
          # min_sum_hessian = 0, poisson_max_delta_step = 0.6 )
          # lgb_model_cv <- lgb.cv(
          #   params = list(
          #     objective = "binary",
          #     boosting="gbdt",
          #     learning_rate = 0.01,
          #     num_leaves = 25,
          #     num_threads = 2),
          #   verbose = -1L,
          #   data = dtrain,
          #   nrounds = 1000,
          #   early_stopping_rounds = 50,
          #   eval_freq = 20,
          #   eval = "auc",
          #   nfold = 5,
          #   stratified = TRUE)
          # best_iter <- lgb_model_cv$best_iter
          # best_iter
          #did some test number of trees (iterations/rounds) do not need to be that high (100 is already more than enough for most ASV)
          
          lgb_model <- lgb.train(
            data=dtrain, 
            verbose=-1, #skip warning messages
            nrounds = nrounds,
            params= list(
              objective = "binary",
              # poisson_max_delta_step = poisson_max_delta_step,
              boosting="gbdt"))
          # if(is.na(lgb_model)){
          #   ranking_list[[paste0("OTU",OTUtoRun)]]<- NA
          #   Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, dev.ratio = NA)
          #   Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
          #   Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
          #   Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
          #   next
          # } else{
          pred_valid <- predict(lgb_model,dvalid_data)
          tryCatch(
            { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[rownames(OTUdata)%in%rownames(ENVeval),OTUtoRun],pred=as.vector(pred_valid),PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
              message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
              message(paste("boot will be ignored, original error :"))
              message(cond)
              return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
            }
          )
        }#end split
          Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_valid,pred_fit)
          Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
          Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
        }#end else
      }#end that ASV
      
      ##############################################################################################################
      ####  Save results      ##################################################################
      # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
      # Model evaluation metrics
      Eval_met_mat <- do.call("rbind", Eval_met_list)
      if (file.exists(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met"))){
        save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
      } else {
        dir.create(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met"))
        save(Eval_met_mat, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
      }
      if (file.exists(paste0("PA/",GtoM,"/Outputs/GBM/Pred_data"))){
        save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
      } else {
        dir.create(paste0("PA/",GtoM,"/Outputs/GBM/Pred_data"))
        save(Pred_data_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
      }
      if (file.exists(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot/"))){
        save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
      } else {
        dir.create(paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot"))
        save(Eval_values_list, file=paste0("PA/",GtoM,"/Outputs/GBM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
      }
    }#end GBM
}#end if PAAB

q("no")