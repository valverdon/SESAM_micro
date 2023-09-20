Fit_GBM<-function(PAAB,arrayID,GtoM,ENVdata,OTUstoRun,OTUdata,test=FALSE){
  #lightGBM bad idea because verfit more than GBM for mall dataset.
  #for light GBM
  library(gbm)
  library(doParallel)
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 

  Hyper.gbm<-function(form,data,distribution="bernoulli",n.cores,n.trees.range=c(10,100,1000,10000),shrinkage.range=c(0.001,0.01,0.1)){
    # data=ENVdata ;n.cores=nthreads
    resultsMatrix<-matrix(NA,4,3)
    rownames(resultsMatrix)<-n.trees.range
    colnames(resultsMatrix)<-shrinkage.range
    for(ntrees in 1:length(n.trees.range)){#ntrees=1
      for(shrink in 1:length(shrinkage.range)){#shrink=1
        resultsMatrix[ntrees,shrink] <- tryCatch(
            {
              cvMod <- gbm(form, data=data, distribution=distribution,n.cores= n.cores,n.trees=n.trees.range[ntrees],shrinkage=shrinkage.range[shrink],bag.fraction=0.5,cv.fold=10)
              cvMod$cv.error[length(cvMod$cv.error)]
              }, error=function(cond){
          return(NA)
            }
        )
      }
    }
    best_shrink<-as.numeric(names(which(apply(resultsMatrix,2,function(X){min(resultsMatrix,na.rm = TRUE)%in%X}))))
    best_ntrees<-as.numeric(names(which(apply(resultsMatrix,1,function(X){min(resultsMatrix,na.rm = TRUE)%in%X}))))
    return(list(best_ntrees=best_ntrees,best_shrink=best_shrink))
  }
  # arrayID=100
  # GtoM="PR"
  # PAAB <- "AB"

  # load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
  # PAAB <- "PA"

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
  nthreads <- 2
  registerDoParallel(cores = nthreads)

  # ENVdataforcor<-ENVdata
  # ENVdataforcor[binvar]<-apply(ENVdataforcor[binvar], 2, function(x) as.numeric(as.character(x)))
  var_cor <- cor(ENVdata, use = "pairwise.complete.obs")
  

  UnivExpl_list <- list()
  preselected_list <- list()
  ranking_list <- list()
  Fit_met_list <- list()
  Fit_data_list <- list()
  Deviances_list <- list()
  Model_list <- list()
  

  for (OTUtoRun in OTUstoRun){ #i=1
    # for (i in 1:2){ print(i)

    if(test){
      print(OTUtoRun)
    }
    # weights <- sapply(OTUdata[,OTUtoRun],function(X){ifelse(X>=1,length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]>=1),length(OTUdata[,OTUtoRun])/sum(OTUdata[,OTUtoRun]<1))})
    #weight of present = nsite/nsite_where_present.
    #weight of absent = nsite/nsite_where_absent.
    explVar_values <- setNames(rep(NA, ncol(ENVdata)), names(ENVdata))
    
    #OTUtoRun<-3
    if(PAAB=="PA"){
      if (sum(OTUdata[,OTUtoRun])/length(OTUdata[,OTUtoRun])>0.95) {#dont bother model it if no absence data (more than 95% presence over all plots)
        UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- setNames(rep(NA, ncol(ENVdata)), names(ENVdata))
        preselected_list[[paste0("OTU",OTUtoRun)]]  <- rep(NA,15)
        ranking_list[[paste0("OTU",OTUtoRun)]] <- list(smooth=NA,param=NA)
        Metrics <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
        Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, Metrics)
        Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- data.frame(obs=OTUdata[,OTUtoRun],fit=rep(1,length(OTUdata[,OTUtoRun])))
        Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, mod.deviance = NA)
        Model_list[[paste0("OTU",OTUtoRun)]] <- NA
      } else{
        
        ##############################################################################################################
        #### Part1 : preselection of covariates #########################################################################
        # Method : Fit  univariate polynomial2 GLMs and take covariates from the 15 best models (avoiding colinear covariates)
        
        #not sure if using weights is good or not yet
        
        # for each covariate, fit an univariate model (polynomial 2nd degree)
        for (j in colnames(ENVdata)) { #j=colnames(ENVdata)[27]
          if(j %in% contvar){
            tryCatch(
              {
                model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ poly(ENVdata[,j], 2), family="binomial", na.action="na.omit")
                # model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ poly(ENVdata[,j], 2), family="quasibinomial", weights = weights, na.action="na.omit")
              }, error=function(cond){
                message(paste0("error with variable ",j," of OTU ", OTUtoRun))
                return(NA)
              })
          } else { #handle binary variable
            model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ factor(ENVdata[,j]), family="binomial", na.action="na.omit")
            # model <- glm(OTUdata[,OTUtoRun] ~ factor(ENVdata[,j]), family="quasibinomial", weights = weights, na.action="na.omit")
          }
          if (all(is.na(model))){ #if fit failed, drop that variable
            explVar_values[j] <- 0 
          } else { #if everything when right, store variance explianed by the model of that variable
            # model <- bam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j]), family="binomial", na.action="na.omit", control=list(nthreads=nthreads))
            explVar_values[j] <- (summary(model)$null.deviance-summary(model)$deviance)/summary(model)$null.deviance*100 
          }
        }
        # store all var expl of univariate models for that ASV
        UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- explVar_values  
        
        #b) preselection of variables, 
        # Method take the most important and remove everything correlated, then take second most important in whats left ...
        thresh <- 0.7 #correlation threshold
        
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
        #store the 15 selected "best" variables for that ASV
        preselected_list[[paste0("OTU",OTUtoRun)]]  <- preselected
        
        
        
        
        
        
        
        ##############################################################################################################
        #### P2 : Model calibration ##################################################################################
        #### a: Fit the full model  #########################################################################
        # Method : Light Generalized boosted regression models
        form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste0(preselected,collapse=" + ")))
        
        hyperparams <- Hyper.gbm(form=form,data=ENVdata,distribution="bernoulli",n.cores=nthreads,n.trees.range=c(10,100,1000,10000),shrinkage.range=c(0.001,0.01,0.1))
        
        
        Mod  <- tryCatch(
          {
            gbm(form, data=ENVdata, distribution="bernoulli",n.cores= nthreads,n.trees=hyperparams$best_ntrees[1],shrinkage=hyperparams$best_shrink[1],bag.fraction=0.5)
          }, error=function(cond){
            message(paste0("error fitting OTU ", OTUtoRun))
            message("Original error:")
            message(cond)
            return(NA)
          }
        )
        
        # Extract results
        if(all(is.na(Mod))){ #if fit failed set all results to NA
          Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(dev_expl=NA, Metrics)
          Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- data.frame(obs=OTUdata[,OTUtoRun],fit=rep(1,length(OTUdata[,OTUtoRun])))
          Model_list[[paste0("OTU",OTUtoRun)]] <- NA
          next
        } else{
        pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=predict.gbm(Mod,ENVdata,n.trees=hyperparams$best_ntrees,type="response"))
        #compute evaluation metrics
        Metrics <- Evalmetrics(obs = pred_expl$obs, pred=pred_expl$fit,PAAB="PA")
        Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- Metrics
        Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
        Model_list[[paste0("OTU",OTUtoRun)]] <- Mod
        }#end ifelse model fit algo worked
      } # end ifelse enough A
    }#end PA
    
    # if(PAAB=="AB"){
    #   load(paste0(PAAB,"/",GtoM,"/data/dataTotSeqSum.Rda")) #total read count data (to go to relative abundance)
    #   
    #   #a) fit univariate models
    #   
    #   for (j in names(ENVdata)) { #j=names(ENVdata)[27]
    #     if(j %in% contvar){
    #       tryCatch(
    #         {
    #           # model <- gam(OTUdata[,OTUtoRun] ~ s(ENVdata[,j], bs="cs", k=4), family="binomial", weights = weights, na.action="na.omit", control=list(nthreads=nthreads, maxit=500))
    #           #I added an offset of the model to convert read count into read relative abundance (log because it has to be the same scale as the response variable that is in log for poisson regression)
    #           model <- glm(as.numeric(OTUdata[,OTUtoRun]) ~ poly(ENVdata[,j], 2) + offset(log(TotSeqSum)), family="poisson", weights = weights, na.action="na.omit")
    #         }, error=function(cond){
    #           message(paste0("error with variable ",j," of OTU ", OTUtoRun))
    #           return(NA)
    #         })
    #     } else {
    #       model <- glm(OTUdata[,OTUtoRun] ~ factor(ENVdata[,j]) + offset(log(TotSeqSum)), family="poisson", weights = weights, na.action="na.omit")
    #     }
    #     if (all(is.na(model))){
    #       explVar_values[j] <- 0 
    #     } else {
    #       explVar_values[j] <- (summary(model)$null.deviance-summary(model)$deviance)/summary(model)$null.deviance*100 
    #     }
    #   }
    #   UnivExpl_list[[paste0("OTU",OTUtoRun)]] <- explVar_values 
    #   
    #   #b ) Fit univariate models 
    #   
    #   thresh <- 0.8
    #   
    #   #get a new object with correlation matrix 
    #   varcor <- var_cor
    #   preselected <- c()
    #   #preselection for all OTUs 1 by one, take the most important and remove everything correlated, then take second most important in whats left ...
    #   for (k in 1:15){
    #     explVar_values <- sort(explVar_values,decreasing=T)
    #     preselected <- c(preselected,names(explVar_values)[1]) #Get the name of the most important var
    #     tokeep <- names(which(abs(varcor[,names(explVar_values)[1]]) < thresh)) #keep non correlated  names(which(abs(varcor[,names(explVar_values)[1]]) > 0.7))
    #     varcor <- varcor[tokeep,tokeep] # reduce the cor matrix
    #     explVar_values <- explVar_values[tokeep] #reduce the explVar_values vector
    #   }
    #   preselected_list[[paste0("OTU",OTUtoRun)]]  <- preselected
    #   
    #   
    #   ##############################################################################################################
    #   #### P2 : Model calibration ##################################################################################
    #   #### a: Fit the full model  #########################################################################
    #   # Method : Light Generalized boosted regression models
    #   form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste0(preselected,collapse=" + "),"+ offset(log(TotSeqSum))"))
    #   
    #   hyperparams <- Hyper.gbm(form=form,data=ENVdata,distribution="poisson",weights=weights,n.cores=nthreads,n.trees.range=c(10,100,1000,10000),shrinkage.range=c(0.001,0.01,0.1))
    #   
    #   
    #   Mod  <- tryCatch(
    #     {
    #       gbm(form, data=ENVdata, distribution="poisson", weights = weights,n.cores= nthreads,n.trees=hyperparams$best_ntrees[1],shrinkage=hyperparams$best_shrink[1],bag.fraction=0.5)
    #     }, error=function(cond){
    #       message(paste0("error fitting OTU ", OTUtoRun))
    #       message("Original error:")
    #       message(cond)
    #       return(NA)
    #     }
    #   )
    # 
    #   # Extract results
    # 
    #   pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=10^(predict.gbm(Mod,ENVdata,n.trees=hyperparams$best_ntrees)+log(TotSeqSum)))
    #   #compute evaluation metrics
    #   Metrics <- Evalmetrics(obs = pred_expl$obs, pred=pred_expl$fit,PAAB="AB")
    #   Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- Metrics
    #   Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
    #   Model_list[[paste0("OTU",OTUtoRun)]] <- Mod
    #   
    #   # } #end what to do with fitted model
    # }#end AB
    }#end that ASV

  ##############################################################################################################
  ####  Save results      ##################################################################
  # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
  
  #Univariate models variance explained
  UnivExpl_mat <- do.call("rbind", UnivExpl_list)
  if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/UnivExpl"))){
    saveRDS(UnivExpl_mat, file=paste0("data/Outputs/GBM/",PAAB,"/UnivExpl/UnivExpl_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/GBM/"))
    dir.create(paste0("data/Outputs/GBM/",PAAB))
    dir.create(paste0("data/Outputs/GBM/",PAAB,"/UnivExpl"))
    saveRDS(UnivExpl_mat, file=paste0("data/Outputs/GBM/",PAAB,"/UnivExpl/UnivExpl_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Selected Variables and shrinkage
  VarSel <- list(preselection=preselected_list,ranking=ranking_list)
  if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/VarSel"))){
    saveRDS(VarSel, file=paste0("data/Outputs/GBM/",PAAB,"/VarSel/VarSel_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/GBM/",PAAB,"/VarSel"))
    saveRDS(VarSel, file=paste0("data/Outputs/GBM/",PAAB,"/VarSel/VarSel_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Model Deviances
  Deviances_mat <- do.call("rbind", Deviances_list)
  if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Deviances"))){
    saveRDS(Deviances_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Deviances/Deviances_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/GBM/",PAAB,"/Deviances"))
    saveRDS(Deviances_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Deviances/Deviances_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Model calibration metrics
  Fit_met_mat <- do.call("rbind", Fit_met_list)
  if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Fit_met"))){
    saveRDS(Fit_met_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Fit_met/Fit_met_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/GBM/",PAAB,"/Fit_met"))
    saveRDS(Fit_met_mat, file=paste0("data/Outputs/GBM/",PAAB,"/Fit_met/Fit_met_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #save fit and cross validation models data (if needed)
  if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Fit_data"))){
    saveRDS(Fit_data_list, file=paste0("data/Outputs/GBM/",PAAB,"/Fit_data/Fit_data_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/GBM/",PAAB,"/Fit_data"))
    saveRDS(Fit_data_list, file=paste0("data/Outputs/GBM/",PAAB,"/Fit_data/Fit_data_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Save model fit
  if (file.exists(paste0("data/Outputs/GBM/",PAAB,"/Models"))){
    saveRDS(Model_list, file=paste0("data/Outputs/GBM/",PAAB,"/Models/Models_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/GBM/",PAAB,"/Models"))
    saveRDS(Model_list, file=paste0("data/Outputs/GBM/",PAAB,"/Models/Models_",GtoM,"_temp_", arrayID, ".rds"))
  }
  return(Model_list)
}#end function
