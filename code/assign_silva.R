#NEED :
# database "../../ASV_data/BA_OTU_265samples.Rdata"
#its fasta version (code to transform below) "../../ASV_data/ASV_taxo/FastaBac.fasta"
# Ref silva database "Z:\projects-unil\SOMETALP\30_data\ASVtables_taxo\Valentin_Verdon_2023_assign_taxo\SILVA_SSU_r138_2019.RData"


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

# load("../../ASV_data/BA_OTU_265samples.Rdata")
load("BA_OTU_265samples.Rdata")
database<-OTU_data_bac_265samples
# nrow(database) #60567
# fastaBac<-paste0(">Seq",1:nrow(database),"\n",rownames(database))
# # write(fastaBac,file="../../ASV_data/ASV_taxo/FastaBac.fasta")


# fas<-"../../ASV_data/ASV_taxo/FastaBac.fasta"
fas <- "FastaBac.fasta"
seqs <- readDNAStringSet(fas)
seqs <- RemoveGaps(seqs)

endloop <- ceiling(length(seqs)/20)
startthisloop<-(arrayID-1)*endloop+1
endthisloop<-(arrayID-1)*endloop+endloop
if (endthisloop > length(seqs)){
  endthisloop <- length(seqs)
} 
if (endthisloop<=0) { q("no")} 

# load("../../ASV_data/ASV_taxo/SILVA_SSU_r138_2019.RData")
load("SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(seqs[startthisloop:endthisloop,],
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

save(ids,file=paste0("BA_assign_silva138_2019_",arrayID,".Rda"))



# load("FU_assign_unite2021.Rda")

# dbConn <- dbConnect(SQLite(),
# ":memory:")

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
#load("BA_assign_silva138_2019.Rda")
ids_full<-ids
kingdom <- sapply(ids_full,
                  function(x) {
                    w <- which(x$rank=="domain")
                    if (length(w) != 1) {
                      "unknown"
                    } else {
                      x$taxon[w]
                    }
                  })
# table(kingdom)
kingdom_conf <-sapply(ids_full,
                      function(x) {
                        w <- which(x$rank=="domain")
                        if (length(w) != 1) {
                          NA
                        } else {
                          x$confidence[w]
                        }
                      })

phylum <- sapply(ids_full,
                 function(x) {
                   w <- which(x$rank=="phylum")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
# table(phylum[kingdom =="Fungi"])

phylum_conf <-sapply(ids_full,
                     function(x) {
                       w <- which(x$rank=="phylum")
                       if (length(w) != 1) {
                         NA
                       } else {
                         x$confidence[w]
                       }
                     })

class <- sapply(ids_full,
                function(x) {
                  w <- which(x$rank=="class")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
# table(class)
class <- sapply(ids_full,
                function(x) {
                  w <- which(x$rank=="class")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
class_conf <-sapply(ids_full,
                    function(x) {
                      w <- which(x$rank=="class")
                      if (length(w) != 1) {
                        NA
                      } else {
                        x$confidence[w]
                      }
                    })

order <- sapply(ids_full,
                function(x) {
                  w <- which(x$rank=="order")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
# table(order)
order_conf <-sapply(ids_full,
                    function(x) {
                      w <- which(x$rank=="order")
                      if (length(w) != 1) {
                        NA
                      } else {
                        x$confidence[w]
                      }
                    })

family <- sapply(ids_full,
                 function(x) {
                   w <- which(x$rank=="family")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
# table(family)
family_conf <- sapply(ids_full,
                      function(x) {
                        w <- which(x$rank=="family")
                        if (length(w) != 1) {
                          NA
                        } else {
                          x$taxon[w]
                        }
                      })


genus <- sapply(ids_full,
                function(x) {
                  w <- which(x$rank=="genus")
                  if (length(w) != 1) {
                    "unknown"
                  } else {
                    x$taxon[w]
                  }
                })
# table(genus)
genus_conf <- sapply(ids_full,
                     function(x) {
                       w <- which(x$rank=="genus")
                       if (length(w) != 1) {
                         NA
                       } else {
                         x$taxon[w]
                       }
                     })

species <- sapply(ids_full,
                  function(x) {
                    w <- which(x$rank=="species")
                    if (length(w) != 1) {
                      "unknown"
                    } else {
                      x$taxon[w]
                    }
                  })

species_conf <- sapply(ids_full,
                       function(x) {
                         w <- which(x$rank=="species")
                         if (length(w) != 1) {
                           NA
                         } else {
                           x$taxon[w]
                         }
                       })
# table(species)
# colnames(FUtaxo)
#check sequences in same order
# all(gsub("ASV","",names(kingdom)) == rownames(FUtaxo))
BAtaxo2019<-data.frame(Seq=rownames(database),Kingdom=kingdom,k_conf=kingdom_conf,Phylum=phylum,p_conf=phylum_conf,Class=class,c_conf=class_conf,Order=order,o_conf=order_conf,Family=family,f_conf=family_conf,Genus=genus,g_conf=genus_conf,Species=species,s_conf=species_conf)
# save(FUtaxo2021,file="../../ASV_data/ASV_taxo/FUtaxo2021.Rda")
# nrow(FUtaxo2021)

save(BAtaxo2019,file=paste0("BAtaxo2019.Rda"))

q("no")