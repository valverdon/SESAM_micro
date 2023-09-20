
Fit_RF<-function(PAAB,arrayID,OTUdata,GtoM,ENVdata,OTUstoRun,test=FALSE){
  #for RF
  # set.seed(seed=08022022)
  # PAAB=PAAB; arrayID=arrayID; GtoM=GtoM; ENVdata=ENVdata; test=TRUE
  library(randomForest)
  library(RRF)
  Hyper.rf<-function(form,data,n.trees.range=c(10,100,1000,10000),PAAB){
    # data=ENVdata ; weights=c("FALSE"=unique(weights[OTUdata[,OTUtoRun]==FALSE]),"TRUE"=unique(weights[OTUdata[,OTUtoRun]==TRUE]))
    if (PAAB=="PA"){
    modqual<-matrix(NA,100,4)
    for(ntrees in 1:length(n.trees.range)){#ntrees=1
      for (i in 1:100){#i=1
        modqual[i,ntrees]<-tryCatch(
          {
          cvMod<-RRF(form, data=data, replace=TRUE,n.trees=n.trees.range[ntrees])
          cvMod$err.rate[nrow(cvMod$err.rate),1]
          }, error=function(cond){
            return(NA)
          }
        )
      }
    }
    resultsVec<-apply(modqual,2,mean)
    names(resultsVec)<-n.trees.range
    best_ntrees<-as.numeric(names(which(resultsVec==min(resultsVec))))
    return(best_ntrees)
    }
    # if (PAAB=="AB"){
    #   modqual<-matrix(NA,100,4)
    #   for(ntrees in 1:length(n.trees.range)){#ntrees=1
    #     for (i in 1:100){
    #       modqual[i,ntrees]<- tryCatch(
    #         {
    #       cvMod<-RRF(form, data=data, replace=TRUE,n.trees=n.trees.range[ntrees],weights=weights)
    #       sqrt(sum((cvMod$predicted-OTUdata[,OTUtoRun])^2)/length(OTUdata[,OTUtoRun]))
    #         }, error=function(cond){
    #           return(NA)
    #         }
    #       )
    #       }
    #   }
    #   resultsVec<-apply(modqual,2,mean)
    #   names(resultsVec)<-n.trees.range
    #   best_ntrees<-as.numeric(names(which(resultsVec==min(resultsVec))))
    #   return(best_ntrees)
    # }
  }
  
  # arrayID=1
  # GtoM="PR"
  # PAAB <- "AB"
  # load(paste0("PA/",GtoM,"/data/ENVdata.Rda"))
  binvar <- c()
  contvar <- c()
  for (covar in names(ENVdata )){#i =names(ENVdata) [31]
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
  
  for (OTUtoRun in OTUstoRun){ #OTUtoRun=OTUstoRun[3]
    # for (i in 1:2){ print(i);OTUtoRun=8963
    # OTUtoRun <- (arrayID-1)*endloop+i #job 1 will run OTU 1-endloop ; job 2 : endloop+1  -  endloop*2 ...  
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
        # Method : Null-space penalization (shrinkage method) selection.
        
        #tested ntrees, increasing not better, stay at 1000
        form<-as.formula(paste0("factor(OTUdata[,OTUtoRun])~ ", paste0(preselected, collapse = " + ")))
        
        #test param RF
        # tested<-matrix(NA,nrow=100,ncol=3)
        # lwt <- list(c("FALSE"=unique(weights[OTUdata[,OTUtoRun]==FALSE]),"TRUE"=unique(weights[OTUdata[,OTUtoRun]==TRUE])),c("FALSE"=0.001,"TRUE"=1),NULL)
        # for(cwt in 1:3){#cwt=1
        # for (teste in 1:100){
        #   MRF<-RRF(form,data=ENVdata,ntree=1000,classwt=lwt[[cwt]])
        #   MRF$confusion
        #   tested[teste,cwt]<-mean(abs(as.numeric(!(OTUdata[,OTUtoRun]))-MRF$votes[,1]))
        # }
        # }
        # boxplot(tested)
        #keep weigths
        # classwt<-c("FALSE"=unique(weights[OTUdata[,OTUtoRun]==FALSE]),"TRUE"=unique(weights[OTUdata[,OTUtoRun]==TRUE]))
        ntree<-Hyper.rf(PAAB=PAAB, form=form,data=ENVdata,n.trees.range=c(10,100,1000,10000))
        #fit regularized random forest model.
        MRF<-RRF(form,data=ENVdata,replace=TRUE,ntree=ntree[1])

        # if(any(is.na(MRF$votes[,2]))){
        #   browser()
        # }
        if(all(is.na(MRF))){ 
          ranking_list[[paste0("OTU",OTUtoRun)]]<- NA
          Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, dev.ratio = NA)
          Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(auc=NA, aucS=NA, rmse=NA, boyce=NA, score=NA, tre=NA, sensit=NA, specif=NA, pospredval=NA, negpredval=NA, jaccar=NA, TSS=NA, accur=NA, kap=NA, SEDI=NA)
          Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
          Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
          next
        } else{
          ############################################
          
          # Extract results
          #here I created a metric that gives 1 when a datapoint is correctly predicted, and gives 0 when uncorrectly predicted
          #Presence datapoint = 0
          #Absence datapoint = 1
          #Predicted present = ValueOfDatapoint - 1
          #Predicted absent  = ValueOfDatapoint + 0
          #Take Absolute Value of the result
          devs <- c(mean.error  = abs(as.numeric(!(OTUdata[,OTUtoRun]))-MRF$votes[,1]))
          Deviances_list[[paste0("OTU",OTUtoRun)]]  <- devs

          #osbserved vs predicted table.
          fitted.val<-predict(MRF, newdata = ENVdata, type = "prob")[,2]
          pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=MRF$votes[,2], RF.50=MRF$predicted)
          
          #compute evaluation metrics and store everything
          Metrics <- Evalmetrics(obs = na.omit(pred_expl)$obs, pred=na.omit(pred_expl)$fit,PAAB="PA")
          Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- Metrics
          Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
          Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
        }
      } #end what to do with fitted model
    }#end that PA
    
    # if(PAAB=="AB"){
    #   load(paste0(PAAB,"/",GtoM,"/data/OTUdata.Rda"))
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
    #   # Method : Null-space penalization (shrinkage method) selection.
    #   
    #   #tested ntrees, increasing not better, stay at 1000
    #   form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste0(preselected, collapse = " + "), "+offset(log(TotSeqSum))"))
    #   # form<-as.formula(paste0("OTUdata[,OTUtoRun]~ ", paste0(names(ENVdata), collapse = " + ")))
    #   
    #   
    #   #test param RF
    #   # tested<-matrix(NA,nrow=100,ncol=3)
    #   # lwt <- list(c("FALSE"=unique(weights[OTUdata[,OTUtoRun]==FALSE]),"TRUE"=unique(weights[OTUdata[,OTUtoRun]==TRUE])),c("FALSE"=0.001,"TRUE"=1),NULL)
    #   # for(cwt in 1:3){#cwt=1
    #   # for (teste in 1:100){
    #   #   MRF<-RRF(form,data=ENVdata,ntree=1000,classwt=lwt[[cwt]])
    #   #   MRF$confusion
    #   #   tested[teste,cwt]<-mean(abs(as.numeric(!(OTUdata[,OTUtoRun]))-MRF$votes[,1]))
    #   # }
    #   # }
    #   ntree<-Hyper.rf(PAAB=PAAB,form=form,data=ENVdata, weights=weights,n.trees.range=c(10,100,1000,10000))
    #   # boxplot(tested)
    #   #keep weigths
    #   MRF<-RRF(form,data=ENVdata,ntree=ntree[1], weights=weights)
    #   #ntree : takes best one, if more than one, fastest one.
    #   
    #   if(is.na(MRF[1])){
    #     ranking_list[[paste0("OTU",OTUtoRun)]]<- NA
    #     Deviances_list[[paste0("OTU",OTUtoRun)]]  <- c(null.deviance = NA, dev.ratio = NA)
    #     Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- c(R2=NA, D2=NA, MAE=NA, MAEs=NA, RMSE=NA, RMSEs=NA, Dspear=NA, Dpear=NA, Pdispersion=NA, R2_scaled=NA, MAE_scaled=NA, RMSE_scaled=NA, Dspear_scaled=NA, Dpear_scaled=NA, Pdisp_scaled=NA)
    #     Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- NA
    #     Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
    #     next
    #   } else{
    #     ############################################
    #     
    #     # Extract results
    #     devs <- c(rsq = MRF$rsq)
    #     Deviances_list[[paste0("OTU",OTUtoRun)]]  <- devs
    #     pred_expl <- data.frame(obs=OTUdata[,OTUtoRun],fit=MRF$predicted)
    #     
    #     #compute evaluation metrics
    #     Metrics <- Evalmetrics(obs = pred_expl$obs[!(is.na(pred_expl$fit))], pred=pred_expl$fit[!(is.na(pred_expl$fit))],PAAB="AB")
    #     Fit_met_list[[paste0("OTU",OTUtoRun)]]  <- Metrics
    #     Fit_data_list[[paste0("OTU",OTUtoRun)]]  <- pred_expl
    #     Model_list[[paste0("OTU",OTUtoRun)]] <- MRF
    #   }
    # }#end AB
  }#end that ASV
  ##############################################################################################################
  ####  Save results      ##################################################################
  # Method : Save 1 file per tpype of result per array, need to run another script to gather all files into 1.
  
  #Univariate models variance explained
  UnivExpl_mat <- do.call("rbind", UnivExpl_list)
  if (file.exists(paste0("data/Outputs/RF/",PAAB,"/UnivExpl"))){
    saveRDS(UnivExpl_mat, file=paste0("data/Outputs/RF/",PAAB,"/UnivExpl/UnivExpl_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/RF/"))
    dir.create(paste0("data/Outputs/RF/",PAAB))
    dir.create(paste0("data/Outputs/RF/",PAAB,"/UnivExpl"))
    saveRDS(UnivExpl_mat, file=paste0("data/Outputs/RF/",PAAB,"/UnivExpl/UnivExpl_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Selected Variables and shrinkage
  VarSel <- list(preselection=preselected_list,ranking=ranking_list)
  if (file.exists(paste0("data/Outputs/RF/",PAAB,"/VarSel"))){
    saveRDS(VarSel, file=paste0("data/Outputs/RF/",PAAB,"/VarSel/VarSel_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/RF/",PAAB,"/VarSel"))
    saveRDS(VarSel, file=paste0("data/Outputs/RF/",PAAB,"/VarSel/VarSel_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Model Deviances
  Deviances_mat <- do.call("rbind", Deviances_list)
  if (file.exists(paste0("data/Outputs/RF/",PAAB,"/Deviances"))){
    saveRDS(Deviances_mat, file=paste0("data/Outputs/RF/",PAAB,"/Deviances/Deviances_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/RF/",PAAB,"/Deviances"))
    saveRDS(Deviances_mat, file=paste0("data/Outputs/RF/",PAAB,"/Deviances/Deviances_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Model calibration metrics
  Fit_met_mat <- do.call("rbind", Fit_met_list)
  if (file.exists(paste0("data/Outputs/RF/",PAAB,"/Fit_met"))){
    saveRDS(Fit_met_mat, file=paste0("data/Outputs/RF/",PAAB,"/Fit_met/Fit_met_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/RF/",PAAB,"/Fit_met"))
    saveRDS(Fit_met_mat, file=paste0("data/Outputs/RF/",PAAB,"/Fit_met/Fit_met_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #save fit and cross validation models data (if needed)
  if (file.exists(paste0("data/Outputs/RF/",PAAB,"/Fit_data"))){
    saveRDS(Fit_data_list, file=paste0("data/Outputs/RF/",PAAB,"/Fit_data/Fit_data_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/RF/",PAAB,"/Fit_data"))
    saveRDS(Fit_data_list, file=paste0("data/Outputs/RF/",PAAB,"/Fit_data/Fit_data_",GtoM,"_temp_", arrayID, ".rds"))
  }
  #Save model fit
  if (file.exists(paste0("data/Outputs/RF/",PAAB,"/Models"))){
    saveRDS(Model_list, file=paste0("data/Outputs/RF/",PAAB,"/Models/Models_",GtoM,"_temp_", arrayID, ".rds"))
  } else {
    dir.create(paste0("data/Outputs/RF/",PAAB,"/Models"))
    saveRDS(Model_list, file=paste0("data/Outputs/RF/",PAAB,"/Models/Models_",GtoM,"_temp_", arrayID, ".rds"))
  }
  return(Model_list)
}#end function
