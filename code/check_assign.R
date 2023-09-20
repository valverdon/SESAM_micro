
.libPaths("/work/FAC/FBM/DEE/aguisan/sometalp/rlib/4.2")
library(dada2)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = '3.16')

library(DECIPHER)
library(RSQLite)
# 
# args <- unlist(commandArgs(trailingOnly = TRUE) )
# # args <- c(41,"PR","All","PA","test")
# # args <- c(200,"PR","GLM","AB")
# 
# arrayID= as.numeric(args[1])  #arrayID=1
# # refpath<-"Y:/FAC/FBM/DEE/aguisan/sometalp/D2c/SometAlp1.0/sh_general_release_dynamic_29.11.2022.fasta"
# # unitedb<-read.table(refpath)
# 
# # load("../../ASV_data/ASV_taxo/FUtaxo.Rda")
# load("../../ASV_data/FU_OTU_232samples.Rdata")
# # 
# # # fastaFung<-paste0(">Seq",rownames(FUtaxo),"\n",FUtaxo$Seq)
# # fastaFung<-paste0(">Seq",1:nrow(OTU_data_fun_232samples),"\n",rownames(OTU_data_fun_232samples))
# # # write(fastaFung,file="../../ASV_data/ASV_taxo/FastaFung.fasta")
# # write(fastaFung,file="FastaFung.fasta")
# 
# # fas<-"../../ASV_data/ASV_taxo/FastaFung.fasta"
fas <- "FastaFung.fasta"
seqs <- readDNAStringSet(fas)
seqs <- RemoveGaps(seqs)

# endloop <- ceiling(length(seqs)/10)
# startthisloop<-(arrayID-1)*endloop+1
# endthisloop<-(arrayID-1)*endloop+endloop
# if (endthisloop > length(seqs)){
#   endthisloop <- length(seqs)-((arrayID-1)*endloop)
# } 
# if (endthisloop<=0) { q("no")} 
# 
# load("../../ASV_data/ASV_taxo/UNITE_v2021_May2021.RData")
# ids <- IdTaxa(seqs[startthisloop:endthisloop,],
#               trainingSet,
#               strand="top", # or "top" if same as trainingSet top2x faster but cant handle reversed sequences (3'-5')
#               threshold=60, # 60 (cautious) or 50 (sensible)
#               processors=NULL) # use all available processors
# 
# save(ids,file=paste0("FU_assign_unite2021_",arrayID,".Rda"))
# 
# q("no")
# seq1<-seqs[1:3,]
# seq2<-seqs[5:7,]
# c(seq1,seq2)
# files <-list.files(path="../../ASV_data/ASV_taxo/", pattern="FU_assign_unite2021_")
files <-list.files(path=getwd(), pattern="FUtaxo2021",full.names = TRUE)

FUtaxo2021_table<-data.frame(Seq=NULL,Kingdom=NULL,k_conf=NULL,Phylum=NULL,p_conf=NULL,Class=NULL,c_conf=NULL,Order=NULL,o_conf=NULL,Family=NULL,f_conf=NULL,Genus=NULL,g_conf=NULL,Species=NULL,s_conf=NULL)

for (i in 1:length(files)){#i=1
  load(files[i])
  
  FUtaxo2021_table<-rbind(FUtaxo2021_table,FUtaxo2021)
}
FUtaxo2021<-FUtaxo2021_table
load("../../ASV_data/FU_OTU_232samples.Rdata")
FUtaxo2021<-FUtaxo2021[match(rownames(OTU_data_fun_232samples),FUtaxo2021$Seq),]
all(test$Seq==rownames(OTU_data_fun_232samples))

save(FUtaxo2021,file=paste0("FUtaxo2021_full.Rda"))

# table(FUtaxo2021$Kingdom)


files <-list.files(path=getwd(), pattern="FU_assign_unite2021",full.names = TRUE)
load(paste0("FU_assign_unite2021_",1,".Rda"))
ids_full<-ids
for (i in 2:length(files)){#i=2
  load(paste0("FU_assign_unite2021_",i,".Rda"))
  ids_full<-c(ids_full,ids)
}
plot(ids_full)

save(ids_full,file="FU_assign_unite2021.Rda")
# print(ids[1:10])
# FUtaxo[FUtaxo$Phylum=="Zygomycota",-1]


##############BACT Silva

files <-list.files(path=getwd(), pattern="BAtaxo2019",full.names = TRUE)

BAtaxo2019_table<-data.frame(Seq=NULL,Kingdom=NULL,k_conf=NULL,Phylum=NULL,p_conf=NULL,Class=NULL,c_conf=NULL,Order=NULL,o_conf=NULL,Family=NULL,f_conf=NULL,Genus=NULL,g_conf=NULL,Species=NULL,s_conf=NULL)

for (i in 1:length(files)){#i=1
  load(files[i])
  
  BAtaxo2019_table<-rbind(BAtaxo2019_table,BAtaxo2019)
}
BAtaxo2019<-BAtaxo2019_table
load("../../ASV_data/BA_OTU_265samples.Rdata")
BAtaxo2019<-BAtaxo2019[match(rownames(OTU_data_bac_265samples),BAtaxo2019$Seq),]
all(BAtaxo2019$Seq==rownames(OTU_data_bac_265samples))

save(BAtaxo2019,file=paste0("BAtaxo2019_full.Rda"))

# table(FUtaxo2021$Kingdom)


files <-list.files(path=getwd(), pattern="BA_assign_silva138_2019",full.names = TRUE)
load(paste0("BA_assign_silva138_2019_",1,".Rda"))
ids_full<-ids
for (i in 2:length(files)){#i=2
  load(paste0("BA_assign_silva138_2019_",i,".Rda"))
  ids_full<-c(ids_full,ids)
}
plot(ids_full)

save(ids_full,file="BA_assign_silva138_2019.Rda")
# print(ids[1:10])
par(mfrow=c(1,1))
load("FU_assign_unite2021.Rda")
plot(ids_full)
