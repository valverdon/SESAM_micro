load("Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/PRtaxo2023_full.Rda")
PRtaxo<-PRtaxo2023
load(paste0("PA/PR/data/OTUdata.Rda"))
table(PRtaxo[match(colnames(OTUdata),PRtaxo$Seq),]$Phylum)
table(PRtaxo$Phylum)
PRtaxo_modeled<-PRtaxo[match(colnames(OTUdata),PRtaxo$Seq),]
save(PRtaxo_modeled,file="Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/PRtaxo2023_modeled.Rda")
# 
# unique(PRtaxo_modeled$Domain)
# unique(PRtaxo_modeled$Phylum)
# plot(table(PRtaxo_modeled$Phylum))
# PRtaxo$Kingdom[PRtaxo$Phylum=="Not_Protist"]

# names_toofew<-names(table(PRtaxo$Phylum)[table(PRtaxo$Phylum)/sum(table(PRtaxo$Phylum))<0.05])
PRtaxo_grouped_phylum<-PRtaxo_modeled
# PRtaxo_grouped_phylum$Phylum[PRtaxo_grouped_phylum$Phylum%in%c(names_toofew)]<-"Others"
plot(table(PRtaxo_grouped_phylum$Phylum))
# if(PRtaxo_grouped_phylum[3019,]$Seq=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT"){
#   # which(colnames(OTUdata)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
#   # any(PRtaxo_grouped_phylum$Seq=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
#   #annoying  Champi = 3019
#   PRtaxo_grouped_phylum<-PRtaxo_grouped_phylum[-3019,]
# }
#nrow(PRtaxo_grouped_phylum)
save(PRtaxo_grouped_phylum,file="../../ASV_data/ASV_taxo/PRtaxo_grouped_phylum.Rda")

# par(mfrow=c(1,1))
load("Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/FUtaxo2023_full.Rda")
FUtaxo<-FUtaxo2021
load(paste0("PA/FU/data/OTUdata.Rda"))
table(FUtaxo$Kingdom)
FUtaxo_modeled<-FUtaxo[match(colnames(OTUdata),FUtaxo$Seq),]
table(FUtaxo_modeled$Kingdom)
FUtaxo<-FUtaxo_modeled
FUtaxo[FUtaxo$Kingdom!="Fungi",2:ncol(FUtaxo)]<-"Not_Fungi"
table(FUtaxo$Phylum)
plot(table(FUtaxo[FUtaxo$Kingdom=="Fungi",]$Phylum))
names_toofew<-names(table(FUtaxo[FUtaxo$Kingdom=="Fungi",]$Phylum)[table(FUtaxo[FUtaxo$Kingdom=="Fungi",]$Phylum)/sum(table(FUtaxo[FUtaxo$Kingdom=="Fungi",]$Phylum))<0.01])
FUtaxo_grouped_phylum<-FUtaxo
FUtaxo_grouped_phylum$Phylum[FUtaxo_grouped_phylum$Phylum%in%names_toofew]<-"Others"
table(FUtaxo_grouped_phylum$Phylum)
# FUtaxo_grouped_phylum[FUtaxo_grouped_phylum$Kingdom!="Fungi",]

# FUtaxo_grouped_phylum[FUtaxo_grouped_phylum$Phylum%in%c("unidentified_60"),]
FUtaxo_grouped_phylum$Phylum[FUtaxo_grouped_phylum$Phylum%in%c("unclassified_Fungi","unidentified_60")]<-"Unidentified_Fungi"
table(FUtaxo_grouped_phylum$Phylum)

barplot(table(FUtaxo_grouped_phylum[FUtaxo_grouped_phylum$Kingdom=="Fungi",]$Phylum))
save(FUtaxo_grouped_phylum,file="../../ASV_data/ASV_taxo/FUtaxo_grouped_phylum.Rda")


######################BA#########################
load("Z:/projects-unil/SOMETALP/30_data/ASVtables_taxo/BAtaxo2019_full.Rda")
BAtaxo<-BAtaxo2019
load(paste0("PA/BA/data/OTUdata.Rda"))
table(BAtaxo$Kingdom)
BAtaxo_modeled<-BAtaxo[match(colnames(OTUdata),BAtaxo$Seq),]
# ARtaxo_grouped_phylum <- BAtaxo
# ARtaxo_grouped_phylum$Phylum[ARtaxo_grouped_phylum$Kingdom=="Bacteria"]<-"Bacteria"
# table(ARtaxo_grouped_phylum$Phylum)
# 
# 
# # #for ASV not assigned to any kingdom, how many per plot (distribution)
# # apply(OTUdata[,colnames(OTUdata)%in%ARtaxo_grouped_phylum[ARtaxo_grouped_phylum$Phylum=="Not_16S",]$Seq],2,summary)[,1:5]
# #  #length of those ASV
# # nchar(ARtaxo_grouped_phylum[ARtaxo_grouped_phylum$Phylum=="Not_16S",]$Seq) #85-150
# # nchar(ARtaxo_grouped_phylum[ARtaxo_grouped_phylum$Phylum=="Bacteria",]$Seq) #85-95
# # #looking at length and the fact it does not align, it should be something else than 16S that was amplified by a primer not specific enough.
# 
# plot(table(BAtaxo$Phylum[ARtaxo_grouped_phylum$Kingdom=="Archaea"]))
# save(ARtaxo_grouped_phylum,file="../../ASV_data/ASV_taxo/ARtaxo_grouped_phylum.Rda")

BAtaxo_grouped_phylum <- BAtaxo_modeled
# BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea"]<-"Archaea"
plot(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria"]))
names_toofewBA<-names(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria"])[table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria"])/sum(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria"]))<0.01])
names_toofewAR<-names(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea"])[table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea"])/sum(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea"]))<0.02])

BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria"&BAtaxo_grouped_phylum$Phylum%in%names_toofewBA]<-"Others"
BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea"&BAtaxo_grouped_phylum$Phylum%in%names_toofewAR]<-"Others"

plot(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria"]))
plot(table(BAtaxo_grouped_phylum$Phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea"]))

save(BAtaxo_grouped_phylum,file="../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda")


