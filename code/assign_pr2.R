#NEED :
# database "../../ASV_data/PR_OTU_179samples.Rdata"
#its fasta version (code to transform below) "../../ASV_data/ASV_taxo/FastaPro.fasta"
# Ref PR2 database "../../ASV_data/ASV_taxo/pr2_version_5.0.0_SSU.decipher.trained.rds"


.libPaths("/work/FAC/FBM/DEE/aguisan/sometalp/rlib/4.2")
library(dada2)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = '3.16')

library(DECIPHER)
library(RSQLite)
# 
args <- unlist(commandArgs(trailingOnly = TRUE) )
# args <- c(41,"PR","All","PA","test")
# args <- c(200,"PR","GLM","AB")

arrayID= as.numeric(args[1])  #arrayID=1

# load("../../ASV_data/PR_OTU_179samples.Rdata")
load("PR_OTU_179samples.Rdata")
database<-OTU_data_pro_179samples
# nrow(database) #11239
# fastaPro<-paste0(">Seq",1:nrow(database),"\n",rownames(database))
# # write(fastaPro,file="../../ASV_data/ASV_taxo/FastaPro.fasta")


# fas<-"../../ASV_data/ASV_taxo/FastaPro.fasta"
fas <- "FastaPro.fasta"
seqs <- readDNAStringSet(fas)
seqs <- RemoveGaps(seqs)

endloop <- ceiling(length(seqs)/10)
startthisloop<-(arrayID-1)*endloop+1
endthisloop<-(arrayID-1)*endloop+endloop
if (endthisloop > length(seqs)){
  endthisloop <- length(seqs)
} 
if (endthisloop<=0) { q("no")} 

# trainingSet<-readRDS("Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/Valentin_Verdon_2023_assign_taxo/pr2_version_5.0.0_SSU.decipher.trained.rds")
trainingSet<-readRDS("pr2_version_5.0.0_SSU.decipher.trained.rds")
plot(trainingSet)
ids <- IdTaxa(seqs,
              trainingSet,
              strand="top", # or "top" if same as trainingSet top2x faster but cant handle reversed sequences (3'-5')
              threshold=60, # 60 (cautious) or 50 (sensible)
              processors=NULL) # use all available processors
# plot(ids)
# ids <- IdTaxa(seqs[1:5,],
#               trainingSet,
#               strand="top", # or "top" if same as trainingSet top2x faster but cant handle reversed sequences (3'-5')
#               threshold=60, # 60 (cautious) or 50 (sensible)
#               processors=NULL) # use all available processors

save(ids,file=paste0("PR_assign_PR2V5_2023_.Rda"))



# load("PR_assign_PR2V5_2023_.Rda")

dbConn <- dbConnect(SQLite(),
":memory:")

# print(ids[1:10])
# FUtaxo[FUtaxo$Phylum=="Zygomycota",-1]
# print(ids[FUtaxo$Phylum=="Zygomycota"])
# plot(ids,trainingSet) #super long
# plot(ids)
# plot(ids[FUtaxo$Phylum=="Zygomycota"],main="Zygomycota")
# plot(ids[FUtaxo$Phylum=="Basidiomycota"],main="Basidiomycota")
# plot(ids[FUtaxo$Phylum=="Ascomycota"],main="Ascomycota")
# plot(ids[FUtaxo$Phylum=="Fungi_unidentified"],main="Fungi_unidentified")
# plot(ids[FUtaxo$Phylum=="Glomeromycota"],main="Glomeromycota")
# plot(ids[FUtaxo$Phylum=="Chytridiomycota"],main="Chytridiomycota")
# plot(ids[FUtaxo$Phylum=="Neocallimastigomycota"],main="Neocallimastigomycota")
# print(ids[FUtaxo$Phylum=="Neocallimastigomycota"])

# ids[,c("rootrank","phylum")]#to look at certain taxo level
# ids[threshold=70] #to change threshold to higher confidence

# kingdom <- sapply(ids,
#                   function(x) {
#                     w <- which(x$rank=="Level 1")
#                     if (length(w) != 1) {
#                       "unknown"
#                     } else {
#                       x$taxon[w]
#                     }
#                   })
domain <- sapply(ids,
                  function(x) {
                    ifelse(is.na(x$taxon[2]),"unknown",x$taxon[2])})
# table(kingdom)
domain_conf  <- sapply(ids,
                        function(x) {
                          ifelse(is.na(x$confidence[2]),"unknown",x$confidence[2])})
kingdom<-sapply(ids,
                    function(x) {
                      ifelse(is.na(x$taxon[3]),"unknown",x$taxon[3])})
kingdom_conf<- sapply(ids,
                          function(x) {
                            x$confidence[3]})

superphylum <- sapply(ids,
                 function(x) {
                   ifelse(is.na(x$taxon[4]),"unknown",x$taxon[4])})
# table(phylum[kingdom =="Fungi"])

superphylum_conf <-sapply(ids,
                     function(x) {
                       x$confidence[4]})

phylum <- sapply(ids,
                function(x) {
                  ifelse(is.na(x$taxon[5]),"unknown",x$taxon[5])})

phylum_conf <-sapply(ids,
                    function(x) {
                      x$confidence[5]})

class <- sapply(ids,
                function(x) {
                  ifelse(is.na(x$taxon[6]),"unknown",x$taxon[6])})
# table(order)
class_conf <-sapply(ids,
                    function(x) {
                      x$confidence[6]})

order <- sapply(ids,
                 function(x) {
                   ifelse(is.na(x$taxon[7]),"unknown",x$taxon[7])})
# table(family)
order_conf <- sapply(ids,
                      function(x) {
                        x$confidence[7]})


family <- sapply(ids,
                function(x) {
                  ifelse(is.na(x$taxon[8]),"unknown",x$taxon[8])})
# table(genus)
family_conf <-  sapply(ids,
                      function(x) {
                        x$confidence[8]})

genus <- sapply(ids,
                  function(x) {
                    ifelse(is.na(x$taxon[9]),"unknown",x$taxon[9])})

genus_conf <- sapply(ids,
                       function(x) {
                         x$confidence[9]})
species <- sapply(ids,
                function(x) {
                  ifelse(is.na(x$taxon[10]),"unknown",x$taxon[10])})

species_conf <- sapply(ids,
                     function(x) {
                       x$confidence[10]})
# table(species)http://127.0.0.1:36749/graphics/plot_zoom_png?width=1958&height=1370
# colnames(FUtaxo)
#check sequences in same order
# all(gsub("ASV","",names(kingdom)) == rownames(FUtaxo))
PRtaxo2023<-data.frame(Seq=rownames(database),Domain=domain,d_conf=domain_conf,Kingdom=kingdom,k_conf=kingdom_conf,Sup_phylum=superphylum,Sup_phylum_conf=superphylum_conf,Phylum2=phylum,p_conf=phylum_conf,Class=class,c_conf=class_conf,Order=order,o_conf=order_conf,Family=family,f_conf=family_conf,Genus=genus,g_conf=genus_conf,Species=species,s_conf=species_conf,Phylum=NA)
# save(FUtaxo2021,file="../../ASV_data/ASV_taxo/FUtaxo2021.Rda")
# nrow(FUtaxo2021)

save(PRtaxo2023,file=paste0("PRtaxo2023_raw.Rda"))


#Before Metazoa, Fungi and Embryophycae removed
#Ill remove them here too, and remove unclassified that are at risk to fall within one of those.
PRtaxo2023[PRtaxo2023$Phylum2=="Metazoa",]
PRtaxo2023[PRtaxo2023$Phylum2=="Fungi",]
table(PRtaxo2023[PRtaxo2023$Sup_phylum=="Opisthokonta",]$Phylum2)#Remove Fungi, Metazoa & unclassified_Opisthokonta at Phylum
PRtaxo2023[PRtaxo2023$Class=="Embryophyceae",]
table(PRtaxo2023[PRtaxo2023$Phylum2=="Streptophyta_X",]$Class)#Remove Embryophyceae & unclassified_Streptophyta_X at Class 

barplot(table(PRtaxo2023$Domain),main="Protist Level 1",cex.names=0.5)
PRtaxo_clean_level1<-PRtaxo2023[PRtaxo2023$Domain!="unclassified_Root",] #remove unclassified

table(PRtaxo_clean_level1$Kingdom)
barplot(table(PRtaxo_clean_level1$Kingdom),main="Protist Level 2",cex.names=0.5)
PRtaxo_clean_level2<-PRtaxo_clean_level1[PRtaxo_clean_level1$Kingdom!="unclassified_Eukaryota",] #remove unclassified

table(PRtaxo_clean_level2$Sup_phylum)
barplot(table(PRtaxo_clean_level2$Sup_phylum),main="Protist Level 3",cex.names=0.5)#nothing to remove here


table(PRtaxo_clean_level2$Phylum2)
barplot(table(PRtaxo_clean_level2$Phylum2),main="Protist Level 4",cex.names=0.5,las=2) ##Remove Fungi, Metazoa & unclassified_Opisthokonta
PRtaxo_clean_level4<-PRtaxo_clean_level2[!(PRtaxo_clean_level2$Phylum2%in%c("unclassified_Opisthokonta","Fungi","Metazoa"))|PRtaxo_clean_level2$Class=="Rozellomycota",] #remove unclassified


table(PRtaxo_clean_level4$Class)
barplot(table(PRtaxo_clean_level4$Class),main="Protist Level 5",cex.names=0.5,las=2) #Remove Embryophyceae & unclassified_Streptophyta_X 
PRtaxo_clean_level5<-PRtaxo_clean_level4[!(PRtaxo_clean_level4$Class%in%c("unclassified_Streptophyta_X","Embryophyceae")),]
barplot(table(PRtaxo_clean_level5$Phylum2),main="Protist Level 4 cleared",cex.names=0.5,las=2) #Remove Embryophyceae & unclassified_Streptophyta_X 
barplot(table(PRtaxo_clean_level5$Kingdom),main="Protist Level 2 cleared",cex.names=0.5)

                                          ########
                                          #Checks#
                                          ########
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Kingdom=="Amoebozoa",]$Sup_phylum) 
# table(PRtaxo_clean_level5[PRtaxo_clean_level5$Kingdom=="Amoebozoa",]$Phylum) 
#Sup_phylum Amoebozoa_X ; Discosea ; Evosea; Tubulinea ;unclassified_Amoebozoa to cluster as Amoebozoa group?
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Apicomplexa",]$Kingdom)
#Phylum Apicomplexa OK (TSAR group, Alveolata)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Apusomonada_X",]$Kingdom)
#Keep phylum Apusomonada ? 3 seq, (sister group to Fungi and Metazoa)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Bigyra",]$Kingdom)
#Phylum Bigyra OK (TSAR group, Stramenopiles)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Breviatea_X",]$Kingdom)
#Keep phylum Breviatea ? 1 seq, (sister group to Fungi and Metazoa)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Centroplasthelida_X",]$Kingdom)
#Keep phylum Centroplasthelida only one group of Haptista (which name to keep?)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Cercozoa",]$Kingdom)
#Phylum Cercozoa OK (TSAR group, Rhizaria)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Chlorophyta_X",]$Kingdom)
#Phylum Chlorophyta OK (plant sister group)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Choanoflagellata",]$Kingdom)
#Phylum Choanoflagellata OK (sister group to Metazoa and Fungi)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Chrompodellids",]$Kingdom)
#Phylum Chrompodellids OK (TSAR group, Alveolata)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Ciliophora",]$Kingdom)
#Phylum Ciliophora OK (TSAR group, Alveolata)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Dinoflagellata",]$Kingdom)
#Kingdom Excavata used as phylum (go into others 1 seq only)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Gyrista",]$Kingdom)
#Phylum Gyrista OK (TSAR group, Stramenopiles)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Ichthyosporea",]$Kingdom)
#Phylum Ichthyosporea OK (sister group to Metazoa and Fungi)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Perkinsea",]$Kingdom)
#Phylum Perkinsea OK (TSAR group,Alveolata)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Streptophyta_X",]$Class)
#Only streptophyta that are not excluded are Klebsormidiophyceae (2seq) & Zygnemophyceae (5seq), I will class them in a Streptophyta Phylum OK (will go in "others")
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Telonemia_X",]$Kingdom)
#Phylum Telonemia OK (TSAR group, T!)

#TSAR group use the 4 main groups, ditch the2 last seq in "others"?

table(PRtaxo_clean_level5$Kingdom)
table(PRtaxo2023[PRtaxo2023$Kingdom=="TSAR",]$Sup_phylum)
table(PRtaxo2023[PRtaxo2023$Phylum2=="Fungi",]$Class)
table(PRtaxo2023[PRtaxo2023$Sup_phylum=="Alveolata",]$Phylum2)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="unclassified_Alveolata",]$Kingdom)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Perkinsea",]$Kingdom)
table(PRtaxo_clean_level5[PRtaxo_clean_level5$Phylum2=="Fungi",]$Class)
PRtaxo2023[PRtaxo2023$Class=="Klebsormidiophyceae",]
PRtaxo2023[PRtaxo2023$Class=="Zygnemophyceae",]

table(PRtaxo_clean_level5$Kingdom)


PRtaxo_clean_level5$Phylum<-factor(apply(PRtaxo_clean_level5,1,function(x){ifelse(x["Kingdom"]=="Amoebozoa",
                                                 "Amoebozoa",
                                                 ifelse(x["Sup_phylum"]%in%c("Chlorophyta","Stramenopiles"),
                                                        x["Sup_phylum"],
                                                        ifelse(x["Phylum2"]%in%c("Apicomplexa","Ciliophora","Cercozoa"),
                                                               x["Phylum2"],
                                                               "Others")))}),
                                   levels=c("Stramenopiles","Apicomplexa","Ciliophora","Cercozoa","Amoebozoa","Chlorophyta","Others"))


plot(table(PRtaxo_clean_level5$Phylum))
PRtaxo2023$Phylum<-apply(PRtaxo2023,1,function(x){ifelse(x["Seq"]%in%PRtaxo_clean_level5$Seq,
                                                         as.character(PRtaxo_clean_level5$Phylum[which(PRtaxo_clean_level5$Seq==x["Seq"])]),
                                                         "Not_Protist")})


save(PRtaxo2023,file="PRtaxo2023.Rda")
