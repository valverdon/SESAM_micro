library(tidyverse)
library(gridExtra)
library(stringi)
#Gathering
PAAB = "AB"
GtoM="PR"
algo="RF"

#TODEBUG TUESDAY
# Fitlist<-list()
# for (i in 1:length(files)) { #for each of the files file=files[2]
#   load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_Met/",files[i])) #load it
#   Fitlist[[i]]<-unlist(Fit_met_mat)  #add what was loaded to the list, the new location takes the same name as the file that was loaded (minus ".Rda")
# }
# which(lapply(Fitlist, ncol)==16)
# Evallist<-list()
# for (i in 1:length(files)) { #for each of the files file=files[2]
#   load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/",files[i])) #load it
#   Evallist[[i]]<-unlist(Eval_met_mat)  #add what was loaded to the list, the new location takes the same name as the file that was loaded (minus ".Rda")
# }

#checking Fit metrics

# files <-list.files(path=paste0("Abundance/",GtoM,"/Outputs/",algo,"/Fit_met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files <- gsub("Fit_met","",files)

for (PAAB in c("PA","AB")){ 
  print(PAAB)
  for (GtoM in c("PR","BA","FU")){ 
    print(GtoM)
    for (algo in c("GAM","GLM","RF","GBM")){
      print(algo)
      # load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/Fit_met",files[1])) #array 1
      load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/Fit_met",files[1])) #array 1
print(ncol(Fit_met_mat))
print(nrow(Fit_met_mat))
# load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/Fit_met",files[3])) #array 1
load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_Met/Fit_met",files[20])) #array 10
print(ncol(Fit_met_mat))
print(nrow(Fit_met_mat))
    }
    }
  }

#checking Eval metrics
# files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files <- gsub("Eval_met","",files)

for (PAAB in c("PA","AB")){
  print(PAAB)
  for (GtoM in c("PR","BA","FU")){
    print(GtoM)
    for (algo in c("GAM","GLM","RF","GBM")){
      print(algo)
      load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/Eval_met",files[1])) #array 1
      print(ncol(Eval_met_mat))
      print(nrow(Eval_met_mat))
      load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/Eval_met",files[3])) #array 10
      print(ncol(Eval_met_mat))
      print(nrow(Eval_met_mat))

    }
    }
  }
files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_Met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
files<- gsub("Fit_met","",files)
for(file in files){
#file=files[1]
  load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_Met/Fit_met",file))
  load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/Eval_met",file))
  print(file)
  print(nrow(Fit_met_mat)==nrow(Eval_met_mat))
}

#actual gathering of data !!!!bug with RF and GBM before 12/10/2022 (date of correction) TSS values are in column "Jaccar" for Rf and GBM (double specif column and no accur column)
for (PAAB in c("PA","AB")){# PAAB="AB"
  for (GtoM in c("PR","BA","FU")){#GtoM="BA"
    print(GtoM)
    for (algo in c("GAM","RF","GLM","GBM")){
     # algo="GLM"
      print(algo)
      files <-list.files(path=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/"), pattern=".Rda") %>% stringr::str_subset(., "temp") #list files ending with ".Rda" and containing "temp"
      files<- gsub("Fit_met","",files)
      Fit_list <- c() #my gathering list Fit metrics
      Eval_list <- c() #my gathering list Eval metrics (median of all boots)
      for (file in files) { #for each of the files file=files[1]
        load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Fit_met/Fit_met",file)) #load it
        Fit_list<-rbind(Fit_list,Fit_met_mat[,(ncol(Fit_met_mat)-14):ncol(Fit_met_mat)])  #add what was loaded to the list, the new location takes the same name as the file that was loaded (minus ".Rda")
        load(paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_met/Eval_met",file)) #load it
        Eval_list<-rbind(Eval_list,unlist(Eval_met_mat))
        # if(nrow(explVar_list)!=nrow(Eval_list)){ #nrow(Fit_list)#nrow(Eval_met_mat) ; nrow(Fit_met_mat)
        #   print(paste0("pb with",file))
        # }
      }
      # nrow(explVar_list) # nrow(Eval_list)  ; rownames(Eval_met_mat) ; rownames(Fit_met_mat)
      fitmet<-as_tibble(Fit_list,.name_repair = "unique")
      fitmet$OTU<-as.numeric(gsub("OTU","",rownames(Fit_list)))
      fitmet <- dplyr::arrange(fitmet,fitmet$OTU)

      evalmet <- as_tibble(Eval_list,.name_repair = "unique")
      # nrow(Eval_list) ;nrow(evalmet); length(rownames(Eval_list))
      evalmet$OTU<-as.numeric(gsub("OTU","",rownames(Eval_list)))
      evalmet <- dplyr::arrange(evalmet,evalmet$OTU)
      print(nrow(fitmet)==nrow(evalmet)) #Model not fitted dont appear in Eval table
      save(fitmet,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
      save(evalmet,file=paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
    }
  }
}
# nrow(Eval_met_mat[str_sub(rownames(Eval_met_mat),1,4)!="pred",])
# Eval_met_mat2<-as.data.frame(rbind(Eval_met_mat[str_sub(rownames(Eval_met_mat),1,4)!="pred",],rep(NA,9)))
# rownames(Eval_met_mat2)[nrow(Eval_met_mat2)] <- "OTU48778"
# Eval_met_mat2$OTU<-as.numeric(gsub("OTU","",rownames(Eval_met_mat2)))
# Eval_met_mat2 <- dplyr::arrange(Eval_met_mat2,Eval_met_mat2$OTU)
# Eval_met_mat <- Eval_met_mat2[1:9]
# Eval_met_mat <- as.matrix(Eval_met_mat)
# save(Eval_met_mat,file=paste0(PAAB,"/",GtoM,"/Outputs/",algo,"/Eval_Met/Eval",file))

# # nrow(Fit_met_mat)
# nrow(Eval_met_mat)
# # Eval_list[10600:10612,]
# # rownames(Eval_list)[10609]
# # rownames(Fit_met_mat) ; rownames(Eval_met_mat)
# # which(str_sub(rownames(Eval_list),1,3)!="OTU")
# # which(rownames(Eval_list)=="OTU10652")
# # rownames(Eval_met_mat)
# # which(apply(Eval_list,1,function(x){all(is.na(x))}))
# # Eval_list[10652,]
# # 29944%%223
# # 29660/223



#checking
for (PAAB in c("PA","AB")){
  for (GtoM in c("PR","BA","FU")){
    print(GtoM)
    for (algo in c("GAM","RF","GLM","GBM")){
      print(algo)
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Fit_Met.Rda"))
      print(nrow(fitmet))
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval_Met.Rda"))
      print(nrow(evalmet))
    }}}
