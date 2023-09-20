
Eval_GLM_Varsel<-function(PAAB,arrayID,GtoM,ENVdata, Model_list, NbRunEval=100, DataSplit=0.8, validation.method = "bootstrap",test=FALSE){
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
  # PAAB <- "PA"
  # load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
  # load(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Models/Models_temp",arrayID,".Rda"))
  # NbRunEval=100
  # validation.method = "bootstrap"
  # test=TRUE
  load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))
  
  binvar <- c()
  contvar <- c()
  for (covar in names(ENVdata)){#i =names(ENVdata) [31]
    if(length(unique(ENVdata[,covar]))>2){
      contvar <- c(contvar,covar)
    } else {
      binvar <- c(binvar,covar)
      ENVdata[,covar]<-as.factor(ENVdata[,covar])
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
  if(test){
    time1<-Sys.time()
  }
  for (i in 1:length(Model_list)){#i=1
    
    # for (i in 1:2){ print(i)
    # OTUtoRun <- as.numeric(substring(names(Model_list)[i],first=4))
    OTUtoRun_name <- names(Model_list)[i]
    if(test){
      print(OTUtoRun_name)
    }
    OTUtoRun<-as.numeric(gsub("OTU","",OTUtoRun_name))
    Mod <- Model_list[[OTUtoRun_name]]
    
    if(DataSplit=="LOO"){#Leave-One-Out Method
      dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "random.split-sampling",DataSplit=(nrow(OTUdata)-1)/nrow(OTUdata))
    }
    dataSplits <- create.datasplittable(NbSites = nrow(OTUdata), NbRunEval = NbRunEval, validation.method = "bootstrap",DataSplit=DataSplit)
    
    #OTUtoRun<-3
    if(PAAB=="PA"){
      if(all(is.na(Mod))){#dont bother evaluate if no model
        Pred_data_list[[paste0("OTU",OTUtoRun)]] <- NA
        Eval_met_list[[paste0("OTU",OTUtoRun)]] <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        
      } else{
        
        pred_GLM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
        Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("AUC", "AUC.S", "RMSE", "Boyce", "Score", "threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Jaccar", "TSS","accur", "Kappa", "SEDI")))
        load(paste0(PAAB,"/",GtoM,"/Outputs/GLM/VarSel/VarSel_temp",arrayID,".Rda"))
        preselected<-VarSel$preselection[[OTUtoRun_name]]
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
          
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)

          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
          
          explVar_values <- setNames(rep(NA, ncol(ENVtrain)), names(ENVtrain))
          ENVdataforcor<-ENVtrain
          ENVdataforcor[binvar]<-apply(ENVdataforcor[binvar], 2, function(x) as.numeric(as.character(x)))
          var_cor <- cor(ENVdataforcor, use = "pairwise.complete.obs")
          
          
          for (j in names(ENVtrain)) { #j=names(ENVdata)[27]
            if(j %in% contvar){
              tryCatch(
                {
                  # model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j], bs="cs", k=4), family="binomial", weights = weights, na.action="na.omit", control=list(nthreads=nthreads, maxit=500))
                  model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ poly(ENVtrain[,j], 2), family="quasibinomial", weights = weights, na.action="na.omit")
                }, error=function(cond){
                  message(paste0("error with variable ",j," of OTU ", OTUtoRun))
                  return(NA)
                })
            } else { #handle binary variable
              model <- glm(OTUdata[,OTUtoRun] ~ factor(ENVtrain[,j]), family="quasibinomial", weights = weights, na.action="na.omit")
            }
            if (all(is.na(model))){ #if fit failed, drop that variable
              explVar_values[j] <- 0 
            } else { #if everything when right, store variance explianed by the model of that variable
              # model <- bam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j]), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
              explVar_values[j] <- (summary(model)$null.deviance-summary(model)$deviance)/summary(model)$null.deviance*100 
            }
          }
          thresh <- 0.8 #correlation threshold
          
          #get a new object with correlation matrix 
          varcor <- var_cor
          preselected <- c()
          for (k in 1:15){
            explVar_values <- sort(explVar_values,decreasing=T) #sort variable by their explanatory power
            preselected <- c(preselected,names(explVar_values)[1]) #Get the name of the best var
            tokeep <- names(which(abs(varcor[,names(explVar_values)[1]]) < thresh)) #names of variable uncorrelated with selected one;  
            #names(which(abs(varcor[,names(explVar_values)[1]]) > 0.7))
            varcor <- varcor[tokeep,tokeep] # reduce the corr matrix, keeping only the variable not selected and not correlated with selected ones
            explVar_values <- explVar_values[tokeep] #reduce the explVar_values vector the same way
          }
          form<-as.formula(paste0("OTUdata~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar])),collapse=" + ")))
          #glmnet needs to transform the function into a matrix format :
          ModMat <- model.matrix(form,data = ENVdata)
          ModMatT <- ModMat[rep(1:nrow(ModMat),times=dataSplits[,f]),]
          ModMatE <- ModMat[dataSplits[,f]==FALSE,]
          OTUtrain <- rep(OTUdata[,OTUtoRun],times=dataSplits[,f])
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
          # ENVtrain <- data.frame(apply(ENVdata,2,function(X){rep(X,times=dataSplits[,f])}))
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
          
          #fit model
          MglmEval  <- tryCatch(
            {
              cv.glmnet(ModMatT, OTUtrain, family="binomial",alpha=1, weights = weights, parallel = TRUE, n.cores = nthreads)#alpha=1 is lasso =0 is ridge
            }, error=function(cond){
              message(paste0("error fitting OTU ", OTUtoRun))
              message("Original error:")
              message(cond)
              return(NA)
            }
          )
          
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
                { Evalvalues[paste0("pred",f),] <- Evalmetrics(obs=OTUeval,pred=fitted.val,PAAB="PA") #offset log bcz nb family (select = T to flatten smoothers of non important variables)
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
    }#PA
    
    if(PAAB=="AB"){
      load(paste0(PAAB,"/",GtoM,"/data/OTUdata.Rda"))
      load(paste0(PAAB,"/",GtoM,"/data/dataTotSeqSum.Rda")) #total read count data (to go to relative abundance)
      
      
      pred_GLM <- matrix(NA,nrow=nrow(ENVdata), ncol = ncol(dataSplits), dimnames = list(rownames(ENVdata), paste0("pred",1:ncol(dataSplits))))
      Evalvalues <- matrix(NA, nrow=ncol(dataSplits), ncol=15, dimnames = list(paste0("pred",1:ncol(dataSplits)),c("R2", "D2", "MAE", "MAEs", "RMSE", "RMSEs", "Dspear", "Dpear", "Pdispersion","R2_scaled","MAE_scaled", "RMSE_scaled", "Dspear_scaled", "Dpear_scaled", "Pdisp_scaled")))
      if(!(is.na(Mod[1]))){#dont bother evaluate if no model
        ##############################################################################################################
        #### P2 : Model evaluation ##################################################################
        
        
        load(paste0(PAAB,"/",GtoM,"/Outputs/GLM/VarSel/VarSel_temp",arrayID,".Rda"))
        preselected<-VarSel$preselection[[OTUtoRun_name]]
        
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
          weights <- sapply(OTUtrain,function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
          
          ENVtrain[binvar] <- lapply(ENVtrain[binvar],factor)
          ENVtrain[contvar] <- lapply(ENVtrain[contvar],as.numeric)
          # sapply(ENVtrain,class)
          #summary(ENVtrain)
          ENVeval <- ENVdata[dataSplits[,f]==FALSE, ]
          OTUeval <- OTUdata[,OTUtoRun][dataSplits[,f]==FALSE]
          TotSeqSumeval <- TotSeqSum[dataSplits[,f]==FALSE]
          
          formT<-as.formula(paste0("OTUtrain~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          formE<-as.formula(paste0("OTUeval~ ", paste(c(paste0("poly(",preselected[preselected%in%contvar],",2)"),paste0(preselected[preselected%in%binvar],collapse=" + ")),collapse=" + ")," -1"))
          
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
    }#end AB
  }#end that ASV
  if(test){
    time2<-Sys.time()
    print(time2-time1)
  }
  ##############################################################################################################
  ####  Save results      ##################################################################
  # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
  
  Eval_met_mat <- do.call("rbind", Eval_met_list)
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met"))){
    save(Eval_met_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met"))
    save(Eval_met_mat, file=paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met/Eval_met_temp", arrayID, ".Rda"))
  }
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Pred_data"))){
    save(Pred_data_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GLM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Pred_data"))
    save(Pred_data_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GLM/Pred_data/Pred_data_temp", arrayID, ".Rda"))
  }
  if (file.exists(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met_allboot/"))){
    save(Eval_values_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
  } else {
    dir.create(paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met_allboot"))
    save(Eval_values_list, file=paste0(PAAB,"/",GtoM,"/Outputs/GLM/Eval_met_allboot/Eval_met_allboot_temp", arrayID, ".Rda"))
  }
  if(test){
    return(Eval_met_mat)
  }
}