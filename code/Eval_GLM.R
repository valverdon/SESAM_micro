
Eval_GLM<-function(PAAB,arrayID,GtoM,ENVdata,OTUdata,OTUstoRun, Model_list, NbRunEval=100, DataSplit=0.8, validation.method = "bootstrap",test=FALSE){
  #for GLM
  library(data.table)
  library(stringi)
  library(glmnet)
  library(MASS)
  library(plyr)
  library(parallel)
  library(doParallel)
  # Model_list=Modlist
  # arrayID=41
  # GtoM="PR"
  # PAAB <- "AB"
  # load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
  # load(paste0("data/Outputs/GLM/",PAAB,"/Models/Models_temp",arrayID,".Rda"))
  # NbRunEval=100
  # validation.method = "bootstrap"
  # test=TRUE

  binvar <- c()
  contvar <- c()
  for (covar in names(ENVdata)){#i =names(ENVdata) [31]
    if(length(unique(ENVdata[,covar]))>2){
      contvar <- c(contvar,covar)
    } else {
      binvar <- c(binvar,covar)
      # ENVdata[,covar]<-as.factor(ENVdata[,covar])
    }
  }
  
  #functions I wrote that will be used
  #Function I wrote that compute some metrics of model evaluation
  #part of the P-A metrics code are from Antoine A.'s work
  source("code/Evalmetrics.R")
  source("code/create_DataSplitTable.R") #function to split dataset for the CV procedure.

  nthreads <- 2
  registerDoParallel(cores = nthreads)
  
 #endthisloop=5
  Eval_met_list <- list()
  Pred_data_list <- list()
  Eval_values_list <-list()

  for (OTUtoRun in OTUstoRun){#OTUtoRun=OTUstoRun[1]
    # for (i in 1:2){ print(i)
    # OTUtoRun <- as.numeric(substring(names(Model_list)[i],first=4))
    if(test){
      print(OTUtoRun)
    }
    Mod <- Model_list[[paste0("OTU",OTUtoRun)]]
    
    if(validation.method=="LOO"){
      NbRunEval <- nrow(OTUdata)
    }
    dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = validation.method, DataSplit=DataSplit)

    
    #OTUtoRun<-3
    if(PAAB=="PA"){
      if(all(is.na(Mod))){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{

        pred_GLM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))
        VarSel<-readRDS(paste0("data/Outputs/GLM/",PAAB,"/VarSel/VarSel_",GtoM,"_temp_",arrayID,".rds"))
        preselected<-VarSel$preselection[[paste0("OTU",OTUtoRun)]]
        ##############################################################################################################
        #### Part1 : preselection of covariates #########################################################################
        # Method : Fit  univariate polynomial2 GLMs and take covariates from the 15 best models (avoiding colinear covariates)

        #not sure if using weights is good or not yet
        
        # for each covariate, fit an univariate model (polynomial 2nd degree)
        for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
          # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
          # ENVtrain <- ENVdata[dataSplits[,f], ]
          # Offset <- log(TotSeqSum[dataSplits[,f]])
          # print(f)
          form<-as.formula(paste0("OTUdata~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          ModMat <- model.matrix(form,data = ENVdata)
          
          ModMatT <- ModMat[rep(1:nrow(ModMat),times=dataSplits[,f]),]
          ModMatE <- ModMat[dataSplits[,f]==FALSE,]
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          # ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          # weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
          
          # ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          # ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)

          # ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
# 
#           formT<-as.formula(paste0("OTUtrain~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
#           formE<-as.formula(paste0("OTUeval~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
#           
          # ModMatT <- model.matrix(formT,data = ENVtrain)
          # ModMatE <- model.matrix(formE,data = ENVeval)
          MglmEval <- tryCatch(
            { cv.glmnet(ModMatT, OTUtrain, family="binomial", alpha=1, parallel = TRUE, n.cores = nthreads)
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
              
              pred_GLM[dataSplits[,f]==FALSE, paste0("pred",f)] <- fitted.val
              if(validation.method!="LOO"){
              tryCatch(
                { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUeval,pred=fitted.val,PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
                }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
                  message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
                  message(paste("boot will be ignored"))
                  return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
                }
              )
              }
            }
          }#end what to do if Modeleval exist
        }#end that split
        if(validation.method=="LOO"){
          pred_GLM2<-apply(pred_GLM,1,function(X){mean(X,na.rm=TRUE)})
          tryCatch(
            { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[,OTUtoRun],pred=pred_GLM2,PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
            }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
              message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
              message(paste("boot will be ignored"))
              return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
            }
          )
        }
        pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
        Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GLM)
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
        Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
      }# end if model exist for that ASV
    }#PA
    
#     if(PAAB=="AB"){
#       load(paste0(PAAB,"/",GtoM,"/data/OTUdata.Rda"))
#       load(paste0(PAAB,"/",GtoM,"/data/dataTotSeqSum.Rda")) #total read count data (to go to relative abundance)
# #Mod<-NA
# 
#         pred_GLM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
#         Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("R2", "D2", "MAE", "MAEs", "RMSE", "RMSEs", "Dspear", "Dpear", "Pdispersion","R2_scaled","MAE_scaled", "RMSE_scaled", "Dspear_scaled", "Dpear_scaled", "Pdisp_scaled")))
#         if (is.na(Mod[1])){
#           pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
#           
#           Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,NA)
#           Eval_met_list[[paste0("OTU",OTUtoRun)]] <- rep(NA, 15)
#           next
#         } else {#dont bother evaluate if no model
#           ##############################################################################################################
#           #### P2 : Model evaluation ##################################################################
#           
# 
#           load(paste0("data/Outputs/GLM/",PAAB,"/VarSel/VarSel_temp",arrayID,".Rda"))
#           preselected<-VarSel$preselection[[paste0("OTU",OTUtoRun)]]
#           
#           #Method  
#           #Either split sampling or bootstrap .632+  Efron & Tibshirani 1997 
#           for (f in 1:ncol(dataSplits)) { #f=1 for the different selection of traning/eval plots.
#             # OTUtrain <- OTUdata[dataSplits[,f], OTUtoRun]
#             # ENVtrain <- ENVdata[dataSplits[,f], ]
#             # Offset <- log(TotSeqSum[dataSplits[,f]])
#             # print(f)
#             form<-as.formula(paste0("OTUdata~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + "), "+ offset(log(TotSeqSum))"," -1 "))
#             ModMat <- model.matrix(form,data = ENVdata)
#             
#             ModMatT <- ModMat[rep(1:nrow(ModMat),times=dataSplits[,f]),]
#             ModMatE <- ModMat[dataSplits[,f]==FALSE,]
#             
#             OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
#             # ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
#             weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
#             
#             # ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
#             # ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
#             # sapply(ENVtrain,class)
#             #summary(ENVtrain)
#             
#             # ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
#             OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
#             # 
#             #           formT<-as.formula(paste0("OTUtrain~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
#             #           formE<-as.formula(paste0("OTUeval~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
#             #           
#             # ModMatT <- model.matrix(formT,data = ENVtrain)
#             # ModMatE <- model.matrix(formE,data = ENVeval)
#             MglmEval <- tryCatch(
#               { cv.glmnet(ModMatT, OTUtrain, family="poisson", weights = weights, alpha=1, parallel = TRUE, n.cores = nthreads)
#               }, error=function(cond){
#                 message(paste0("error with boot ",f," of OTU ", OTUtoRun))
#                 message(cond)
#                 message(paste("boot will be ignored"))
#                 return(NA)
#               }
#             )
#             # modGAMp <- gam(OTUtrain ~ s(TcoldQ)+s(PdryM)+s(Trange)+s(sRad)+s(TPI)+s(slope)+s(pH)+s(TOC)+s(clay)+offset(Offset), 
#             #                data=ENVtrain, family="poisson")
#             # modGBM <- gbm(OTUtrain ~ TcoldQ+PdryM+Trange+sRad+TPI+slope+pH+TOC+clay+offset(Offset), 
#             #               data=ENVtrain, distribution="poisson", n.trees = 2000, interaction.depth = 3, shrinkage = 0.01)
#             if(is.na(MglmEval[1])){
#               Eval_met_list[[paste0("OTU",OTUtoRun)]]  <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
#               Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
#               next
#             } else{
#               # Extract results
#               lambda<- list(onese=MglmEval$lambda.1se,min=MglmEval$lambda.min)
#               if(MglmEval$lambda[1]!=MglmEval$lambda.1se){#take lambda1se(lambda with error at 1se from minimum) if its not the first, 
#                 fitted.val <- predict(MglmEval,newx=ModMatE, s=MglmEval$lambda.1se, type="response", newoffset = log(TotSeqSumeval))
#               } else{#take lambda that minimize error
#                 fitted.val <- predict(MglmEval,newx=ModMatE, s=MglmEval$lambda.min, type="response", newoffset = log(TotSeqSumeval))
#               }
#               if(empty(fitted.val)){
#                 Eval_met_list[[paste0("OTU",OTUtoRun)]]  <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
#                 Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
#                 next
#               } else{
#                 colnames(fitted.val) <- "fit"
#                 
#                 pred_GLM[rownames(ModMatE), paste0("pred",f)] <- fitted.val
#                 if(validation.method!="LOO"){
#                   tryCatch(
#                     { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUeval,pred=fitted.val,PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
#                     }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
#                       message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
#                       message(paste("boot will be ignored"))
#                       return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
#                     }
#                   )
#                 }
#               }
#             }#end what to do if Modeleval exist
#           }#end that split
#           if(validation.method=="LOO"){
#             pred_GLM2<-apply(pred_GLM,1,function(X){mean(X,na.rm=TRUE)})
#             tryCatch(
#               { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUdata[,OTUtoRun],pred=pred_GLM2,PAAB="AB") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
#               }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
#                 message(paste0("error evaluating boot ",f," of OTU ", OTUtoRun))
#                 message(paste("boot will be ignored"))
#                 return(Evalvalues[paste0("pred",f),] <- rep(NA, 15))
#               }
#             )
#           }
#           pred_eval <- data.frame(obs=OTUdata[,OTUtoRun])
#           
#           Pred_data_list[[paste0("OTU",OTUtoRun)]]  <- cbind(pred_eval,pred_GLM)
#           Eval_met_list[[paste0("OTU",OTUtoRun)]] <- apply(Evalvalues, 2, function(X){median(X,na.rm=T)})#median of all boots (avoid -inf problems that would appear with mean)
#           Eval_values_list[[paste0("OTU",OTUtoRun)]] <- Evalvalues
#         }
#     }#end AB
  }#end that ASV
  ##############################################################################################################
  ####  Save results      ##################################################################
  # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
  
  Eval_met_mat <- do.call("rbind", Eval_met_list)
  if (file.exists(paste0("data/Outputs/GLM/",PAAB,"/Eval_met"))){
    save(Eval_met_mat, file=paste0("data/Outputs/GLM/",PAAB,"/Eval_met/Eval_met_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0("data/Outputs/GLM/",PAAB,"/Eval_met"))
    save(Eval_met_mat, file=paste0("data/Outputs/GLM/",PAAB,"/Eval_met/Eval_met_temp", arrayID, ".Rda"))
  }
  if (file.exists(paste0("data/Outputs/GLM/",PAAB,"/Pred_data"))){
    save(Pred_data_list, file=paste0("data/Outputs/GLM/",PAAB,"/Pred_data/Pred_data_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0("data/Outputs/GLM/",PAAB,"/Pred_data"))
    save(Pred_data_list, file=paste0("data/Outputs/GLM/",PAAB,"/Pred_data/Pred_data_temp", arrayID, ".Rda"))
  }
  if (file.exists(paste0("data/Outputs/GLM/",PAAB,"/Eval_met_allboot/"))){
    save(Eval_values_list, file=paste0("data/Outputs/GLM/",PAAB,"/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0("data/Outputs/GLM/",PAAB,"/Eval_met_allboot"))
    save(Eval_values_list, file=paste0("data/Outputs/GLM/",PAAB,"/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
  }
}

