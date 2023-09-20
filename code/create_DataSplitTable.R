#Create some datasplittable for SESAMALP 
#Modified from Heidi Mod's code, (added bootstrap options)



#ToDo: Block CV ?
##########################################################################
#create a datasplit table NbSite x RbRunEval to prepare for CV procedure #
#two approach for now : random split-sampling and bootstrap .634+

# Create the DataSplitTable # NbSites=250 ; NbRunEval=100 ; DataSplit=80
create.datasplittable <- function(NbSites, #number of sites to choose from
                                  NbRunEval, #number of samplings to do (number of CV loops)
                                  DataSplit, #for split-sampling, percentage of data to put in training set
                                  validation.method #split-sample or .634? 
                                  ){
  #create result matrix
  DataSplitTable <- matrix(data=FALSE, nrow=NbSites, ncol=NbRunEval)

  #Random Split-sample approach
  if(validation.method=="random.split-sample"){
    for(i in 1:NbRunEval){
      #fill matrix with TRUE whereSite have been selected for each CV loop
      DataSplitTable[sample(1:NbSites,round(DataSplit*NbSites), replace=F),i] <- TRUE
    }
  }
  if(validation.method=="LOO"){
    for(i in 1:NbRunEval){
      #fill matrix with TRUE whereSite have been selected for each CV loop
      DataSplitTable[-i,i] <- TRUE
    }
  }
  # #Constrained split-sample : each site will be taken datasplit % for training
  # if(validation.method=="constrained.split-sample"){
  #   grouper <- sample(rep(1:NbRunEval,each=ceiling(NbSites/NbRunEval)),NbSites, replace = FALSE)
  #   iner <- round(DataSplit*NbRunEval/100)
  #   for(i in 1:NbRunEval){
  #     DataSplitTable[which(grouper %in% ((i:(i+iner-1)%%NbRunEval)+1)),i] <- TRUE
  #   }
  # }
  if(validation.method=="bootstrap"){
    for(i in 1:NbRunEval){
      #select Nsite sites with replacement
    Selection <- sample(1:NbSites,NbSites, replace=T)
    for(j in 1:NbSites){
      DataSplitTable[j,i] <- sum(j == Selection)
      #fill the matrix giving the number of time site was selected
      }
    }
  }
  # The matrix is then filled with 0 if site not selected, with the number of time if not
  # max possible value is 1 in split-sample context.
  return(DataSplitTable)
}

