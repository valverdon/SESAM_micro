#Code that do the null model randomisation of data
library(tidyverse)
load("PA/BA/Outputs/OTUdata.Rda")

OTUdata_BA<-OTUdata

load("PA/FU/Outputs/OTUdata.Rda")
OTUdata_FU<-OTUdata

load("PA/PR/Outputs/OTUdata.Rda")
OTUdata_PR<-OTUdata


ASV_npresence_BA<-apply(OTUdata_BA,2,sum)
ASV_npresence_FU<-apply(OTUdata_FU,2,sum)
ASV_npresence_PR<-apply(OTUdata_PR,2,sum)


# sites_BA<-substr(rownames(OTUdata_BA),1,5)
# sites_FU<-substr(rownames(OTUdata_FU),1,5)
# sites_PR<-substr(rownames(OTUdata_PR),1,5)
# sites_PR%in%sites_BA
# sites_FU%in%sites_BA
cas_BA<-unique(ASV_npresence_BA[(nrow(OTUdata_BA)/20) < ASV_npresence_BA & ASV_npresence_BA <    ( nrow(OTUdata_BA)-(nrow(OTUdata_BA)/20)  )  ])[order(unique(ASV_npresence_BA[(nrow(OTUdata_BA)/20) < ASV_npresence_BA & ASV_npresence_BA <    ( nrow(OTUdata_BA)-(nrow(OTUdata_BA)/20)  )  ]))  ]
# all possibility exist between 13 and 237 (5%)
cas_FU<-unique(ASV_npresence_FU[(nrow(OTUdata_FU)/20) < ASV_npresence_FU & ASV_npresence_FU <    ( nrow(OTUdata_FU)-(nrow(OTUdata_FU)/20)  )  ])[order(unique(ASV_npresence_FU[(nrow(OTUdata_FU)/20) < ASV_npresence_FU & ASV_npresence_FU <    ( nrow(OTUdata_FU)-(nrow(OTUdata_FU)/20)  )  ]))  ]
# 
cas_PR<-unique(ASV_npresence_PR[(nrow(OTUdata_PR)/20) < ASV_npresence_PR & ASV_npresence_PR <    ( nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20)  )  ])[order(unique(ASV_npresence_PR[(nrow(OTUdata_PR)/20) < ASV_npresence_PR & ASV_npresence_PR <    ( nrow(OTUdata_PR)-(nrow(OTUdata_PR)/20)  )  ]))  ]
# 166-17


#BA
setlist<-as.list(cas_BA)
names(setlist)<-cas_BA
for (i in cas_BA){#i=cas_BA[1]
  mymatrix<-apply(matrix(NA,nrow=nrow(OTUdata_BA),ncol=100),2,function(X){ifelse(1:nrow(OTUdata_BA)%in%sample(1:nrow(OTUdata_BA),i,replace=FALSE),TRUE,FALSE)})
  rownames(mymatrix)<-rownames(OTUdata_BA)
  colnames(mymatrix)<-paste0("prev",i,"rep",1:100)
  setlist[[which(names(setlist)==i)]]<-mymatrix
}
OTUdata_BA_null<-do.call("cbind",setlist)

OTUdata<-OTUdata_BA_null
save(OTUdata,file="PA/NULL_BA/data/OTUdata.Rda")


#FU
setlist<-as.list(cas_FU)
names(setlist)<-cas_FU
for (i in cas_FU){#i=cas_FU[1]
  mymatrix<-apply(matrix(NA,nrow=nrow(OTUdata_FU),ncol=100),2,function(X){ifelse(1:nrow(OTUdata_FU)%in%sample(1:nrow(OTUdata_FU),i,replace=FALSE),TRUE,FALSE)})
  rownames(mymatrix)<-rownames(OTUdata_FU)
  colnames(mymatrix)<-paste0("prev",i,"rep",1:100)
  setlist[[which(names(setlist)==i)]]<-mymatrix
}
OTUdata_FU_null<-do.call("cbind",setlist)

OTUdata<-OTUdata_FU_null
save(OTUdata,file="PA/NULL_FU/data/OTUdata.Rda")




#PR
setlist<-as.list(cas_PR)
names(setlist)<-cas_PR
for (i in cas_PR){#i=cas_PR[1]
  mymatrix<-apply(matrix(NA,nrow=nrow(OTUdata_PR),ncol=100),2,function(X){ifelse(1:nrow(OTUdata_PR)%in%sample(1:nrow(OTUdata_PR),i,replace=FALSE),TRUE,FALSE)})
  rownames(mymatrix)<-rownames(OTUdata_PR)
  colnames(mymatrix)<-paste0("prev",i,"rep",1:100)
  setlist[[which(names(setlist)==i)]]<-mymatrix
}
OTUdata_PR_null<-do.call("cbind",setlist)

OTUdata<-OTUdata_PR_null
save(OTUdata,file="PA/NULL_PR/data/OTUdata.Rda")




##########################
# BA->PR
OTUdata_BA_166<-OTUdata_BA[substr(rownames(OTUdata_BA),1,5)%in%substr(rownames(OTUdata_PR),1,5),]
OTUdata<-OTUdata_BA_166
save(OTUdata,file="PA/BA_166/data/OTUdata.Rda")

