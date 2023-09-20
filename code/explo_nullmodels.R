#Add null models reuslts in Evaluation tabs + check TSS values of Null models across prevalence of data.

GtoM="PR"
Mod="RF"
load(paste0("PA/NULL_",GtoM,"/data/",Mod,"/Eval_Met.Rda"))
load(paste0("PA/NULL_",GtoM,"/data/",Mod,"/Fit_Met.Rda"))
Fitnull<-fitmet
Evalnull<-evalmet


load(paste0("PA/",GtoM,"/data/",Mod,"/Eval_Met.Rda"))
load(paste0("PA/",GtoM,"/data/",Mod,"/Fit_Met.Rda"))
Fit<-fitmet
Eval<-evalmet


# prev1
load(paste0("PA/NULL_",GtoM,"/data/OTUdata.Rda"))
prevs<-unique(apply(OTUdata,2,sum)) #list of prevalences for that group
load(paste0("PA/NULL_",GtoM,"/data/OTUdata.Rda"))
prev_group<-nrow(OTUdata)
# OTUtoRun<-100
# prevOTU<-sum(OTUdata[,OTUtoRun]) #get the prevalence of that ASV
# prevOTUpos<-which(prevs==prevOTU)#position of prevOTU in prevs vector
# 
#which null models to look at ?
# postolook<-(100*(prevOTUpos-1)+1):(100*(prevOTUpos-1)+100)
#
# fit100val<-Fitnull[Fitnull$OTU%in%postolook,]
# eval100val<-Evalnull[Evalnull$OTU%in%postolook,]
# Fit1val<-Fit[Fit$OTU==OTUtoRun,]
# Eval1val<-Eval[Eval$OTU==OTUtoRun,]
# 
# 
# par(mfrow=c(1,2))
# boxplot(fit100val$TSS, ylim=c(0,1),main=paste0("Fit",GtoM,"_",Mod))
# abline(h=fit100val$TSS[order(fit100val$TSS)[95]],col="red")
# abline(h=Fit1val$TSS,lwd=2)
#
# boxplot(eval100val$TSS,  ylim=c(0,1),main=paste0("Eval",GtoM,"_",Mod))
# abline(h=eval100val$TSS[order(eval100val$TSS)[95]],col="red")
# abline(h=Eval1val$TSS,lwd=2)
# par(mfrow=c(1,1))


#What are the null models achieving according to prev values

par(mfrow=c(3,4))
for (GtoM in c("PR","FU","BA","BA_166")){
  # for (GtoM in c("BA_166")){#GtoM="BA
    
  for (Mod in c("GLM","GBM","GAM","RF")){
    # Mod="GLM"
    #load null model results
    if(GtoM=="BA_166"){#BA_166 uses PR plots so uses PR null models
      load(paste0("PA/NULL_PR/data/",Mod,"/Eval_Met.Rda"))
      load(paste0("PA/NULL_PR/data/",Mod,"/Fit_Met.Rda"))
      load(paste0("PA/NULL_PR/data/OTUdata.Rda"))
      prevs<-unique(apply(OTUdata,2,sum)) #list of prevalences for that group
      Fitnull<-fitmet
      Evalnull<-evalmet
    } else {
      load(paste0("PA/NULL_",GtoM,"/data/",Mod,"/Eval_Met.Rda"))
      load(paste0("PA/NULL_",GtoM,"/data/",Mod,"/Fit_Met.Rda"))
      load(paste0("PA/NULL_",GtoM,"/data/OTUdata.Rda"))
      prevs<-unique(apply(OTUdata,2,sum)) #list of prevalences for that group
      Fitnull<-fitmet
      Evalnull<-evalmet
    }
  
  load(paste0("PA/",GtoM,"/data/",Mod,"/Eval_Met.Rda"))
  load(paste0("PA/",GtoM,"/data/",Mod,"/Fit_Met.Rda"))
  Fit<-fitmet
  Eval<-evalmet
  #annoying  Champi = 3019 in protists Evals
  if(nrow(Eval)==3161){
  Eval<-Eval[-3019,]
  Eval$OTU[3020:length(Eval$OTU)]<-3010:(length(Eval$OTU)-1)
  }
  
  # prev1
  
  load(paste0("PA/",GtoM,"/data/OTUdata.Rda"))
  prev_group<-nrow(OTUdata)
  
  Null_means<-c(rep(NA,length(prevs)))
  Null_95s<-c(rep(NA,length(prevs)))
  Null_5s<-c(rep(NA,length(prevs)))
  names(Null_means)<-names(Null_95s)<-as.character(prevs)
  
  
  for (i in 1:length(prevs)){
    #which null models to look at ?
    postolook<-(100*(i-1)+1):(100*(i-1)+100)
    eval100val<-Evalnull[Evalnull$OTU%in%postolook,]
    ordered_eval100val_TSS <- eval100val$TSS[order(eval100val$TSS)]
    Null_means[i]<-mean(ordered_eval100val_TSS,na.rm=TRUE)
    Null_95s[i]<-ordered_eval100val_TSS[round(length(na.omit(ordered_eval100val_TSS))/100*95)]
    Null_5s[i]<-ordered_eval100val_TSS[round(length(na.omit(ordered_eval100val_TSS))/100*5)]
    
  }
  
  plot(prevs/prev_group,Null_means,pch=20,col="blue",ylim=c(-0.1,1),ylab="TSS value",main= paste0(" ",GtoM,Mod), xlab="prevalence")
  points(prevs/prev_group,Null_95s,pch=20,col="red",ylim=c(-0.1,1))
  points(prevs/prev_group,Null_5s,pch=20,col="red",ylim=c(-0.1,1))
  
  
  
  
  #which Eval1Val are above their eval100val 95th best value (is the value better than what a null model would achieve)
  
  
  Eval$TSS_sign<-NA
  Eval$TSS_adj<-NA
  ncol(OTUdata)
  for (i in Eval$OTU){
    OTUtoRun<-i #which ASV to look at
    OTUprev<-sum(OTUdata[,OTUtoRun]) #whats its prevalence
    if (!(OTUprev%in%prevs)){
      next
    }
    prev_to_look<-which(prevs==OTUprev) #position of that value in the prevs vector
    postolook<-(100*(prev_to_look-1)+1):(100*(prev_to_look-1)+100) #in the null model list, where are the models corresponding to that prevalence
    eval100val<-Evalnull[Evalnull$OTU%in%postolook,]
    ordered_eval100val_TSS <- eval100val$TSS[order(eval100val$TSS)]
  
    thresh<-ordered_eval100val_TSS[round(length(na.omit(ordered_eval100val_TSS))/100*95)]
    Eval1val<-Eval[Eval$OTU==OTUtoRun,]
    Eval$TSS_sign[i]<-ifelse(Eval1val$TSS>thresh,TRUE,FALSE)
    Eval$TSS_adj[i]<-(Eval1val$TSS-thresh)/(1-thresh)
  }
  
  #Percentage of ASV that are better predicted than by chance (with 5% percentile ...)
  print(paste0("percentage of ASV better predicted than by chance",GtoM,Mod))
  print(sum(Eval$TSS_sign,na.rm=TRUE)/length(Eval$TSS_sign))
  save(Eval,file=paste0("PA/",GtoM,"/data/",Mod,"/Eval.Rda"))
  }
}
summary(Eval$TSS_adj)
summary(Eval$TSS)
all(is.na(Eval$TSS_adj)==is.na(Eval$TSS))
par(mfrow=c(1,1))