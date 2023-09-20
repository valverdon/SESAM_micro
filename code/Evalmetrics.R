#################################################################
# Function to measure distance between data and model         ###
#################################################################
#two part : count model evaluation and P-A model evaluation

#some functions usefull in P-A evaluation : 
# symmetric extremal dependence index. From A.A.
# sensitivity or hit rate
HitRate <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return (a / (a + c))
}

# 1 minus hit rate
OneMHR <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return(1 - (a / (a + c)))
}

# log of hit rate
logHR <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(HitRate(a = a, b = b, c = c, d = d)))
}

# helper
# false alarm
FalseAlarm <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return (b / (b + d))
}

# 1 minus false alarm
OneMFA <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return (1 - (b / (b + d)))
}

# log of false alarm
logFA <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(FalseAlarm(a = a, b = b, c = c, d = d)))
}

# SEDI main function
sedi <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  
  if (isTRUE(x = (max(a,b,c,d) <= 1) && sum(a,b,c,d) <= 1)) { # no value can be computed
    return(list(NaN))  
  }
  
  if (isTRUE(x = (min(a,b,c,d) < 1))) { # catch all matrices which contain at least one zero
    if (isTRUE(x = (sum(a,b) < 1) &&  (sum(c,d) > 1) )) {
      return(list(0))
    }
    if (isTRUE(x = (sum(c,d) < 1) &&  (sum(a,b) > 1) )) {
      return(list(0))
    }
    if (isTRUE(x = (sum(b,c) < 1) &&  (sum(a,d) > 1) )) {
      return(list(1))
    }
    if (isTRUE(x = (sum(a,d) < 1) &&  (sum(b,c) > 1) )) {
      return(list(0))
    }
    
    if (isTRUE(x = (a < 1 && min(b,c,d) >=1))) { # only true presences equal to zero
      return(list(
        (logFA(a = a+1e-09, b = b, c = c, d = d) - logHR(a = a+1e-09, b = b, c = c, d = d) - log(OneMFA(a = a+1e-09, b = b, c = c, d = d)) + log(OneMHR(a = a+1e-09, b = b, c = c, d = d)))
        /
          (logFA(a = a+1e-09, b = b, c = c, d = d) + logHR(a = a+1e-09, b = b, c = c, d = d) + log(OneMFA(a = a+1e-09, b = b, c = c, d = d)) + log(OneMHR(a = a+1e-09, b = b, c = c, d = d)))))  
    }
    
    if (isTRUE(x = (d < 1 && min(a,b,c) >=1))) { # only true absences (also pseudo-absences, incl. background) equal to zero
      return(list(
        (logFA(a = a, b = b, c = c, d = d+1e-09) - logHR(a = a, b = b, c = c, d = d+1e-09) - log(OneMFA(a = a, b = b, c = c, d = d+1e-09)) + log(OneMHR(a = a, b = b, c = c, d = d+1e-09)))
        /
          (logFA(a = a, b = b, c = c, d = d+1e-09) + logHR(a = a, b = b, c = c, d = d+1e-09) + log(OneMFA(a = a, b = b, c = c, d = d+1e-09)) + log(OneMHR(a = a, b = b, c = c, d = d+1e-09)))))  
    }
    
    if (isTRUE(x = (b < 1 && min(a,c,d) >=1))) { # only zero commission errors
      return(list(
        (logFA(a = a, b = b+1e-09, c = c, d = d) - logHR(a = a, b = b+1e-09, c = c, d = d) - log(OneMFA(a = a, b = b+1e-09, c = c, d = d)) + log(OneMHR(a = a, b = b+1e-09, c = c, d = d)))
        /
          (logFA(a = a, b = b+1e-09, c = c, d = d) + logHR(a = a, b = b+1e-09, c = c, d = d) + log(OneMFA(a = a, b = b+1e-09, c = c, d = d)) + log(OneMHR(a = a, b = b+1e-09, c = c, d = d)))))  
    }
    
    if (isTRUE(x = (c < 1 && min(a,b,d) >=1))) { # only zero omission errors
      return(list(
        (logFA(a = a, b = b, c = c+1e-09, d = d) - logHR(a = a, b = b, c = c+1e-09, d = d) - log(OneMFA(a = a, b = b, c = c+1e-09, d = d)) + log(OneMHR(a = a, b = b, c = c+1e-09, d = d)))
        /
          (logFA(a = a, b = b, c = c+1e-09, d = d) + logHR(a = a, b = b, c = c+1e-09, d = d) + log(OneMFA(a = a, b = b, c = c+1e-09, d = d)) + log(OneMHR(a = a, b = b, c = c+1e-09, d = d)))))  
    }
  }
  else { # all regular cases with reasonable values across the confusion matrix
    return (list(
      (logFA(a = a, b = b, c = c, d = d) - logHR(a = a, b = b, c = c, d = d) - log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
      /
        (logFA(a = a, b = b, c = c, d = d) + logHR(a = a, b = b, c = c, d = d) + log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
    ))  
  }
}

#TSS
tss <- function(a,b,c,d){a/(a+c)+d/(b+d)-1}


Evalmetrics <- function(obs, pred,PAAB){ #obs=OTUdata[rownames(OTUdata)%in%names(predCV),OTUtoRun] #pred =fitted.val #coefs=17
  library(modEvA)
  library(ROCR)
  source("code/myboyce.R") #just removed a line that was bugging the function from
  #library(ecospat) #for normal boyce #check with Flavien if boyce works now.
  
  #P-A part : based on confusion matrix a b c d 
  if(PAAB=="PA"){#if obs = PA
    #calculate threshold
    # if(crit=="maxTSS"){ #only option for now, modified from A.A. Valpar work
    

      order_pred <- pred[order(pred,decreasing = T)]
      order_obs <- obs[order(pred,decreasing = T)]
      cpa <- cumsum(order_obs)
      
      # Reduce to sensible range   #why?
      nsns <- which(cpa==max(cpa))
      if (length(nsns)>1) {
        nsns <- nsns[-1]
        order_pred <- order_pred[-nsns]
        order_obs <- order_obs[-nsns]
        cpa <- cpa[-nsns]
      }
      

      #TSS
      tsss <- apply(cbind(1:length(cpa),cpa),1,function(x,y,z){
        out <- tss(x[2],x[1]-x[2],z-x[2],y-x[1]-(z-x[2]))
        return(out)
      },y = length(obs), z = sum(obs))
      #a= cbind(1:length(cpa),cpa)[,2] ; b= cbind(1:length(cpa),cpa)[,1]-cbind(1:length(cpa),cpa)[,2] ; c = sum(obs)-cbind(1:length(cpa),cpa)[,2] ; d = length(obs)-cbind(1:length(cpa),cpa)[,1]-(sum(obs)-cbind(1:length(cpa),cpa)[,2])
      tre <- order_pred[which.max(tsss)]#max tss
      if (length(tre) == 0){#treshold = prevalence
          tre = mean(obs)
        } 

    
    # Calculate threshold-dependent metrics
    pred_bin <- ifelse(pred<=tre,FALSE,TRUE) #binarization
    tb <- table(factor(pred_bin,levels=c(TRUE,FALSE)),factor(as.logical(obs),levels=c(TRUE,FALSE))) #confusion matrix
    a <- tb[1,1] ; b <- tb[1,2] ; c <- tb[2,1] ; d <- tb[2,2] ; n <- a+b+c+d
    
    # Sensitivity
    sensit <- a/(a+c)
    
    # Specificity
    specif <- d/(b+d)
    
    # Positive predictive value
    pospredval <- a/(a+b)
    
    # Negative predictive value
    negpredval <- d/(c+d)
    
    # Jaccard index
    jaccar <- a/(a+b+c)
    
    # True skill statistic
    TSS <- a/(a+c)+d/(b+d)-1
    
    # Accuracy
    accur <- (a+d)/(a+b+c+d)
      
    # Kappa
    kap <- ((a+d)/n-((a+b)*(a+c)+(c+d)*(d+b))/n^2)/(1-((a+b)*(a+c)+(c+d)*(d+b))/n^2)

    #SEDI symmetric extremal dependence index
    SEDI <- sedi(a,b,c,d)
    
    # Boyce
    # obs = predicted suitability values at presence points 
    # fit = predicted suitability values at all points

    # boyce <- try(ecospat.boyce(fit=pred, obs=pred[obs==1], PEplot = FALSE)$cor, TRUE)
    boyce <-try(myboyce(fit=pred, obs=pred[obs==1], PEplot = FALSE)$cor, TRUE)
    if(!is.numeric(boyce)){
      boyce<-NA}

    # AUC   ##############################???????????????????????????????????????????
    z <- prediction(as.numeric(pred),as.numeric(obs))
    auc <- performance(z,measure="auc")@y.values[[1]]
    rmse <- performance(z,measure="rmse")@y.values[[1]]

    # Somer's AUC
    aucS <- (2*auc)-1 # rescale AUC to -1 +1
    
    # # Consensus score
    score <- mean(c(aucS, boyce, TSS), na.rm=TRUE)
  
    
    # Return results
    res <- c(auc, aucS, rmse, boyce, score, tre, sensit, specif, pospredval, negpredval, jaccar, TSS, accur, kap, SEDI)
    names(res)=c("auc", "aucS", "rmse", "boyce", "score", "tre", "sensit", "specif", "pospredval", "negpredval", "jaccar", "TSS", "accur", "kap", "SEDI")
    return(unlist(res))

  
  } #end PA part
  
  #Count part
  if(PAAB=="AB") {#if obs = abundance
    R2 <- 1-((sum((obs-pred)^2))/
               (sum((obs - mean(obs))^2))) #Maybe not the best R2 metric
    D2 <-  tryCatch({
      Dsquared(obs=obs,pred=pred, family = 'poisson') #the one we keep !!!
    }, error=function(cond){#get D2 for the projection on the non-selected sites of the boot.
      return(NA)
    })
    
    #from Waldock etal 2021
    # accuracy : MAE how close the mean prediction is to the data
    MAE <- sum(abs(obs-pred))/length(obs)
    MAEs <- MAE/mean(obs) #Waldock 2021, divide by mean to "scale" the value and have something that can be compared between models 
    #"prediction within x% of observed value"
    RMSE <- sqrt(sum((obs - pred)^2)/length(obs))
    RMSEs <- RMSE/mean(obs) #=coef of variation 
    
    # discrimination
    Dspear <- cor(obs,pred, method="spearman")
    Dpear <-  cor(obs,pred, method="pearson")
    
    # precision not possible to measure ? from Waldock2021
    Pdispersion <- sqrt(var(pred))/sqrt(var(obs))
    
    #same on sccaled abundance, () not very usefull ) ^^
    scaled_obs <- obs/max(obs)
    scaled_pred <-pred/max(obs)
    R2_scaled <- 1-((sum((scaled_obs-scaled_pred)^2))/
               (sum((scaled_obs - mean(scaled_obs))^2)))
    MAE_scaled <- sum(abs(scaled_obs-scaled_pred))/length(scaled_obs)
    RMSE_scaled <- sqrt(sum((scaled_obs - scaled_pred)^2)/length(scaled_obs))
    # MAAPE <- sum(atan(abs(obs-pred)))/length(obs)
    RMSE_normalized <- RMSE/(max(obs)-min(obs))
    #Otto, S.A. (2019, Jan.,7). How to normalize the RMSE [Blog post]. Retrieved from https://www.marinedatascience.co/blog/2019/01/07/normalizing-the-rmse/
    Dspear_scaled <- cor(scaled_obs,scaled_pred, method="spearman")
    Dpear_scaled <-  cor(scaled_obs,scaled_pred, method="pearson")
    Pdisp_scaled <- sqrt(var(scaled_pred))/sqrt(var(scaled_obs))
    
    res <- c(R2=R2, D2=D2, MAE=MAE, MAEs=MAE_scaled, RMSE=RMSE_scaled, RMSEs=RMSEs, Dspear=Dspear, Dpear=Dpear, Pdispersion=Pdispersion, R2_scaled=R2_scaled, MAE_scaled=MAE_scaled, RMSE_scaled=RMSE_scaled, Dspear_scaled=Dspear_scaled, Dpear_scaled=Dpear_scaled, Pdisp_scaled=Pdisp_scaled)
    return(res)
  }
}

#old Evalmetrics
#Evalmetrics <- function(obs, pred){ #obs=OTUdata[rownames(OTUdata)%in%names(predCV),OTUtoRun] #pred =predCV #coefs=17
# R2 <- 1-((sum((obs-pred)^2))/
#            (sum((obs - mean(obs))^2)))
# D2 <-  Dsquared(obs=obs,pred=pred, family = 'poisson') #the one we keep !!!
# #from Waldock etal 2021
# # accuracy : MAE how close the mean prediction is to the data
# MAE <- sum(abs(obs-pred))/length(obs)
# MAEs <- MAE/mean(obs) #Waldock 2021, divide by mean to "scale" the value and have something that can be compared between models 
# #"prediction within x% of observed value"
# RMSE <- sqrt(sum((obs - pred)^2)/length(obs))
# RMSEs <- RMSE/mean(obs^2) #=coef of variation, max = sqrt(length(obs)-1)
# 
# # discrimination
# Dspear <- cor(obs,pred, method="spearman")
# Dpear <-  cor(obs,pred, method="pearson")
# 
# # precision not possible to measure ? from Waldock2021
# Pdispersion <- sqrt(var(pred))/sqrt(var(obs))
# 
# a <- c(R2=R2, D2=D2, MAE=MAE, MAEs=MAEs, RMSE=RMSE, RMSEs=RMSEs, Dspear=Dspear, Dpear=Dpear, Pdispersion=Pdispersion)
# return(a)
# }
