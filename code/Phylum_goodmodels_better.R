#For presence- Absence
#Creates figures of proportion of each phylum in overall data and in models better than  null
#Creates figures of proportion of each phylum ASV that have better models than null
#Do the same for TSS > 0.5
#For abundances
#Creates figures of proportion of each phylum in overall data and in models with cor > 0.2
#Creates figures of proportion of each phylum ASV that have model with cor >0.2
#Do the same for cor > 0.5
#whogoodmodels
library(tidyverse)
library(ggplot2)
library(reshape2) 
library(grid)
library(gridExtra)
library(paletteer)

library(palettetown)
# pokedex(51, 11)
# ?pokedex
# paletteer_d("palettetown::vulpix") for BA : zapdos? raticate?
# paletteer_d("palettetown::meowth")
setcolors<-list(BA=c("#F8E080FF","#E0B040FF","#F8F8B0FF","#986800FF","#C87810FF","#E89830FF","#F8D000FF","#D05838FF","#C8A000FF","#F87050FF","#F8F890FF","gray69","#cacaca"),
                FU=c("#5E0B15", "#90323D", "#D9CAB3", "#BC8034", "gray69", "#cacaca"),
                AR=c("#1E441E", "#119822", "#2A7221", "#52A300", "gray69", "#cacaca"),
                PR=c("#03396c","#6497b1","#011f4b", "#b3cde0", "#3B8EA5", "#005b96" ,"#cacaca"),
                dataset=c(BA='#FFD700',AR='#5CB800',FU='#6B4C62',PR='#457EB0'))
#For presence- Absence
PAAB="PA"

for(Mod in c("GBM","GAM","RF","GLM")){#Mod="GLM"
  # group="PR"
  for (group in c("PR","BA","FU")){#group="BA"
    # for (PAAB in c("PA","AB")){

    print(paste0(group,PAAB,Mod))
    load(paste0(PAAB,"/",group,"/data/",Mod,"/Eval.Rda"))
    
    
    if(group=="BA"){
      # load(paste0("../../ASV_data/ASV_taxo/",substr(group,1,2),"taxo.Rda"))
      # Eval<- cbind(BAtaxo,Eval)
      load(paste0("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda"))
      if(nrow(BAtaxo_grouped_phylum)!=nrow(Eval)){#securité
        message("problem nrow")
      }
      evalmet_BA<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria",],Eval[BAtaxo_grouped_phylum$Kingdom=="Bacteria",])
      evalmet_BA$Phylum<-factor(evalmet_BA$Phylum,levels=c(names(table(evalmet_BA$Phylum[!(evalmet_BA$Phylum%in%c("Others","unclassified_Bacteria"))])),"Others","unclassified_Bacteria"))
      evalmet_AR<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea",],Eval[BAtaxo_grouped_phylum$Kingdom=="Archaea",])
      evalmet_AR$Phylum<-factor(evalmet_AR$Phylum,levels=c(names(table(evalmet_AR$Phylum[!(evalmet_AR$Phylum%in%c("Others","unclassified_Archaea"))])),"Others","unclassified_Archaea"))
      #mean TSS
      # mean(evalmet_BA$TSS_adj,na.rm=TRUE) #0.27
      # sd(evalmet_BA$TSS_adj,na.rm=TRUE) # 0.20
      # mean(evalmet_AR$TSS_adj,na.rm=TRUE) # 0.34
      # sd(evalmet_AR$TSS_adj,na.rm=TRUE) #0.16
      # nrow(evalmet_BA)
      # nrow(taxo_grouped_phylum_BA)

      # unique(taxo_grouped_phylum_BA$Species)
      #proportion of each phylum in good modeled ASV pool
      
      goodASV_phylum_OK_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_02_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.2&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.2&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_04_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.4&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.4&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_05_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_06_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.6&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.6&!(is.na(evalmet_BA$TSS))]))
      goodASV_phylum_08_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.8&!(is.na(evalmet_BA$TSS))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.8&!(is.na(evalmet_BA$TSS))]))
      # number instead of proportions
      # summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])
      # Acidobacteriota      Actinobacteriota          Bacteroidota      Bdellovibrionota           Chloroflexi            Firmicutes 
      # 1171                   556                   368                    72                   528                    77 
      # Myxococcota       Patescibacteria       Planctomycetota        Proteobacteria     Verrucomicrobiota                Others 
      # 194                   189                   576                  1404                   285                   221 
      # unclassified_Bacteria 
      # 1334 
      # summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))]) #
      # Crenarchaeota        Euryarchaeota        Halobacterota        Iainarchaeota               Others unclassified_Archaea 
      # 20                    0                    0                    0                    2                    2 
      
      goodASV_phylum_OK_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_02_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.2&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.2&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_04_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.4&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.4&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_05_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_06_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.6&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.6&!(is.na(evalmet_AR$TSS))]))
      goodASV_phylum_08_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.8&!(is.na(evalmet_AR$TSS))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.8&!(is.na(evalmet_AR$TSS))]))
      
      #proportion of each phylum in all modeled ASV pool
      AllASV_phylum_BA<-summary(factor(evalmet_BA$Phylum))/sum(summary(factor(evalmet_BA$Phylum)))
      AllASV_phylum_AR<-summary(factor(evalmet_AR$Phylum))/sum(summary(factor(evalmet_AR$Phylum)))
      
      comparison_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_OK_BA,goodASV_phylum_02_BA,goodASV_phylum_04_BA,goodASV_phylum_06_BA,goodASV_phylum_08_BA))
      comparison_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_OK_AR,goodASV_phylum_02_AR,goodASV_phylum_04_AR,goodASV_phylum_06_AR,goodASV_phylum_08_AR))

      #proportion of good modeled ASV in each phylum
      goodASV_OK_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_sign&!(is.na(evalmet_BA$TSS))])/summary(factor(evalmet_BA$Phylum))
      goodASV_OK_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_sign&!(is.na(evalmet_AR$TSS))])/summary(factor(evalmet_AR$Phylum))
      
      #proportion of verygood modeled ASV in each phylum
      goodASV_05_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$TSS_adj>0.5&!(is.na(evalmet_BA$TSS))])/summary(factor(evalmet_BA$Phylum))
      if(length(goodASV_05_BA)==0){
        goodASV_05_BA<-goodASV_OK_BA*0
      }
      goodASV_05_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$TSS_adj>0.5&!(is.na(evalmet_AR$TSS))])/summary(factor(evalmet_AR$Phylum))
      if(length(goodASV_05_AR)==0){
        goodASV_05_AR<-goodASV_OK_AR*0
      }

      molten_good_BA <- melt(data.frame(goodASV_OK=goodASV_OK_BA,phylum=names(goodASV_OK_BA),goodASV_05=goodASV_05_BA))
      molten_good_BA$phylum<-factor(molten_good_BA$phylum,levels=levels(evalmet_BA$Phylum))
      labels_BA <- paste(unique(molten_good_BA$phylum), " ( n =",c(table(evalmet_BA$Phylum)[-which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))],table(evalmet_BA$Phylum)[which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))]), ")")
      molten_good_AR <- melt(data.frame(goodASV_OK=goodASV_OK_AR,phylum=names(goodASV_OK_AR),goodASV_05=goodASV_05_AR))
      molten_good_AR$phylum<-factor(molten_good_AR$phylum,levels=levels(evalmet_AR$Phylum))
      labels_AR <- paste(unique(molten_good_AR$phylum), " ( n =",c(table(evalmet_AR$Phylum)[-which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))],table(evalmet_AR$Phylum)[which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))]), ")")
      
      p1_BA<-ggplot(comparison_BA, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
        scale_fill_manual("Bacteria phylum",values = setcolors$BA)+ theme_bw() +
        scale_x_discrete(labels=c("AllASV_phylum_BA" = "Whole dataset", "goodASV_phylum_OK_BA" = "TSS >TSSnull", "goodASV_phylum_02_BA" = "TSSadj>0.2", "goodASV_phylum_04_BA" = "TSSadj>0.4", "goodASV_phylum_06_BA" = "TSSadj>0.6", "goodASV_phylum_08_BA" = "TSSadj>0.8"))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.title = element_blank(),
              axis.title.x = element_blank())
      p2_BA<-ggplot(molten_good_BA,aes(x=variable,y=value,fill=phylum,group=phylum))+
        geom_bar(position="dodge",stat="identity")+
        ylim(0,1) + xlab("") + ylab("Proportion of models")+
        scale_fill_manual("Bacteria phylum",values = setcolors$BA,
                          label=labels_BA)+ theme_bw() +
        theme_bw() + scale_x_discrete(labels=c("goodASV_OK" = "TSS >TSSnull", "goodASV_05" = "TSSadj>0.5"))+ 
        theme_classic()+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=6),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0.5,0,0,0),"cm"),
              legend.title = element_text(size=7,colour = setcolors$dataset[group]),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.3,"line"),
              legend.position = c(1,1.1),
              legend.justification = c("right", "top"),
              legend.background = element_blank())
      # grid.arrange(p1_BA,p2_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_02"))
      # p3_BA<-ggplot(comparison_05_BA, aes(x=Var1, y=value, fill=Var2)) +
      #   geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
      #   scale_fill_manual(values = setcolors$BA)
      
      # p4_BA<-ggplot(data.frame(value=goodASV_05_BA,phylum=names(goodASV_05_BA)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
      #   ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
      #   scale_fill_manual(values = setcolors$BA)
      # grid.arrange(p3_BA,p4_BA,nrow=1,top=paste0("BA_",Mod,"_",PAAB,"_05"))
      p1_AR<-ggplot(comparison_AR, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
        scale_fill_manual("Archaea phylum",values = setcolors$AR)+ theme_bw() +
        scale_x_discrete(labels=c("AllASV_phylum_AR" = "Whole dataset", "goodASV_phylum_OK_AR" = "TSS >TSSnull", "goodASV_phylum_02_AR" = "TSSadj>0.2", "goodASV_phylum_04_AR" = "TSSadj>0.4", "goodASV_phylum_06_AR" = "TSSadj>0.6", "goodASV_phylum_08_AR" = "TSSadj>0.8"))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.title = element_blank(),
              axis.title.x = element_blank())
      
      p2_AR<-ggplot(molten_good_AR,aes(x=variable,y=value,fill=phylum,group=phylum))+
        geom_bar(position="dodge",stat="identity")+
        ylim(0,1) + xlab("") + ylab("")+
        scale_fill_manual("Archaea phylum",values = setcolors$AR,
                          label=labels_AR)+ 
        theme_bw() + scale_x_discrete(labels=c("goodASV_OK" = "TSS >TSSnull", "goodASV_05" = "TSSadj>0.5"))+ 
        theme_classic()+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=6),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0.5,0,0,0),"cm"),
              legend.title = element_text(size=7,colour = setcolors$dataset["AR"]),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.3,"line"),
              legend.position = c(1,1.1),
              legend.justification = c("right", "top"))
      # grid.arrange(p1_AR,p2_AR,nrow=1,top=paste0("AR_",Mod,"_",PAAB,"_02"))
      # p3_AR<-ggplot(comparison_05_AR, aes(x=Var1, y=value, fill=Var2)) +
      #   geom_bar(stat="identity",show.legend = FALSE) + xlab("") + ylab("")+
      #   scale_fill_manual(values = setcolors$AR)
      # 
      # p4_AR<-ggplot(data.frame(value=goodASV_05_AR,phylum=names(goodASV_05_AR)),aes(x=phylum,y=value,fill=phylum))+geom_bar(stat="identity")+
      #   ylim(0,1) + theme(legend.key.size = unit(0.2, 'cm'),axis.text.x = element_blank()) + xlab("") + ylab("")+
      #   scale_fill_manual(values = setcolors$AR)
      # grid.arrange(p3_AR,p4_AR,nrow=1,top=paste0("AR_",Mod,"_",PAAB,"_05"))
      
      
    }else{

      load(paste0("../../ASV_data/ASV_taxo/",group,"taxo_grouped_phylum.Rda"))
      taxo_grouped_phylum_temp<-eval(parse(text=paste0(group,"taxo_grouped_phylum")))
      #remove not fungi
      table(taxo_grouped_phylum_temp$Phylum)
      
      if(nrow(taxo_grouped_phylum_temp)!=nrow(Eval)){#securité
        message("problem")
      }
      

        if(group=="PR"){
          evalmet<-cbind(Eval[taxo_grouped_phylum_temp$Phylum!="Not_Protist",],taxo_grouped_phylum_temp[taxo_grouped_phylum_temp$Phylum!="Not_Protist",])
          evalmet$Phylum<-factor(evalmet$Phylum,levels=c(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])[which(!(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])%in%c("Others")))],"Others"))
        }
        if(group=="FU"){
          evalmet<-cbind(Eval[taxo_grouped_phylum_temp$Phylum!="Not_Fungi",],taxo_grouped_phylum_temp[taxo_grouped_phylum_temp$Phylum!="Not_Fungi",])
          evalmet$Phylum<-factor(evalmet$Phylum,levels=c(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])[which(!(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])%in%c("Others","Unidentified_Fungi")))],"Others","Unidentified_Fungi"))
        }
      table(evalmet$Phylum)
      evalmet$Phylum<-droplevels(evalmet$Phylum)
      #mean
      # mean(evalmet$TSS_adj,na.rm=TRUE) #FU 0.19  #PR  0.04
      # sd(evalmet$TSS_adj,na.rm=TRUE) #FU 0.20  #PR   0.14
      # levels(evalmet$Phylum)
      # evalmetGLM<-evalmet
      # fitmetGLM<-fitmet
      goodASV_phylum_OK<-summary(factor(evalmet$Phylum)[evalmet$TSS_sign&!(is.na(evalmet$TSS_sign))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_sign&!(is.na(evalmet$TSS_sign))]))
      goodASV_phylum_02<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.2&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.2&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_04<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.4&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.4&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_05<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_06<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.6&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.6&!(is.na(evalmet$TSS_adj))]))
      goodASV_phylum_08<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.8&!(is.na(evalmet$TSS_adj))])/sum(summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.8&!(is.na(evalmet$TSS_adj))]))
      # number instead of proportions
      # summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS))]) 
      #Fungi GLM
      # Ascomycota     Basidiomycota     Glomeromycota Mortierellomycota            Others      Unidentified 
      # 589                 43                 37                 89                  3                203 
      #Protists GLM
      # SAR      Amoebozoa Archaeplastida         Others           NA's 
      #        1              0              1              0             25 

      
      AllASV_phylum<-summary(factor(evalmet$Phylum))/sum(summary(factor(evalmet$Phylum)))

      comparison<-melt(rbind(AllASV_phylum,goodASV_phylum_OK,goodASV_phylum_02,goodASV_phylum_04,goodASV_phylum_06,goodASV_phylum_08))
      comparison$value[is.na(comparison$value)]<-0
      
      #proportion of good modeled ASV in each phylum
      goodASV_OK<-summary(factor(evalmet$Phylum)[evalmet$TSS_sign&!(is.na(evalmet$TSS_sign))])/summary(factor(evalmet$Phylum))

      #proportion of verygood modeled ASV in each phylum
      goodASV_05<-summary(factor(evalmet$Phylum)[evalmet$TSS_adj>0.5&!(is.na(evalmet$TSS_adj))])/summary(factor(evalmet$Phylum))
      if(length(goodASV_05)==0){
        goodASV_05<-goodASV_OK*0
      }

      
      molten_good <- melt(data.frame(goodASV_OK=goodASV_OK,phylum=names(goodASV_OK),goodASV_05=goodASV_05))
      molten_good$phylum<-factor(molten_good$phylum,levels=levels(evalmet$Phylum))
      # labels <- paste(unique(molten_good$phylum), " (n =",c(table(evalmet$Phylum)[-which(names(table(evalmet$Phylum))%in%c("Others","Fungi_unidentified"))],table(evalmet$Phylum)[which(names(table(evalmet$Phylum))%in%c("Others","Fungi_unidentified"))]), ")",sep="")
      labels <- paste(unique(molten_good$phylum), " (n =",table(evalmet$Phylum), ")",sep="")
      # sum(!(is.na(evalmet$TSS)))
      assign(paste0("p1_", group), ggplot(comparison, aes(x=Var1, y=value, fill=Var2)) +
               geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
               scale_fill_manual(paste0(ifelse(group=="PR","Protist phylum","Fungi phylum")),values = setcolors[[group]])+ theme_bw() +
               scale_x_discrete(labels=c("AllASV_phylum" = "Whole dataset", "goodASV_phylum_OK" = "TSS > TSSnull", "goodASV_phylum_02" = "TSSadj>0.2", "goodASV_phylum_04" = "TSSadj>0.4", "goodASV_phylum_06" = "TSSadj>0.6", "goodASV_phylum_08" = "TSSadj>0.8"))+ 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(size=10),
                     axis.text.y = element_text(size=10),
                     legend.title = element_blank(),
                     axis.title.x = element_blank())
             )#p1_FU
      if(group=="PR"){plotmarg<-unit(c(0,0.5,0,0),"cm")}else{plotmarg<-unit(c(0,0,0,0.5),"cm")}
      assign(paste0("p2_", group), ggplot(molten_good,aes(x=variable,y=value,fill=phylum,group=phylum))+
               geom_bar(position="dodge",stat="identity")+
               ylim(0,1) + xlab("") + ylab(ifelse(group=="PR","","Proportion of models"))+
               scale_fill_manual(paste0(ifelse(group=="PR","Protist phylum","Fungi phylum")),values = setcolors[[group]],
                                 label=labels)+ 
               theme_bw() + scale_x_discrete(labels=c("goodASV_OK" = "TSS > TSSnull", "goodASV_05" = "TSSadj > 0.5"))+ 
               theme_classic()+ 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=10),
                     axis.text.y = element_text(size=6),
                     axis.text.x = element_text(size=10),
                     plot.margin = unit(c(0,0,0,0),"cm"),
                     legend.title = element_text(size=7,colour = setcolors$dataset[group]),
                     legend.text = element_text(size=6),
                     legend.key.size = unit(0.3,"line"),
                     legend.position = c(.95,1),
                     legend.justification = c("right", "top"))
             )
    }#p2_FU
  }

  pdf(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Pylum_proportion",Mod,"_",PAAB,".pdf"),width=1961,height=1500)
  grid.arrange(p2_BA,p2_AR,p2_FU,p2_PR,nrow=2)
  grid.arrange(p1_BA,p1_AR,p1_FU,p1_PR,nrow=4)
  dev.off()
  
  png(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_proportion",Mod,"_",PAAB,".png"),res=300,width=1961,height=1500)
  grid.arrange(p2_BA,p2_AR,p2_FU,p2_PR,nrow=2)
  dev.off()
  
  png(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_relative_richness",Mod,"_",PAAB,".png"),res=300,width=1961,height=1500)
  grid.arrange(p1_BA,p1_AR,p1_FU,p1_PR,nrow=4)
  dev.off()
}

PAAB="AB"
#For Abundances




  # for (PAAB in c("PA","AB")){
for(Mod in c("RF","GLM","GBM","GAM")){#Mod="GLM"
    for (group in c("PR","FU","BA")){#group="PR"
# group="PR"
# 
# Mod="GBM"

    print(paste0(group,PAAB,Mod))
    load(paste0(PAAB,"/",group,"/data/",Mod,"/Eval_Met.Rda"))

    if(group=="BA"){
      load(paste0("../../ASV_data/ASV_taxo/BAtaxo_grouped_phylum.Rda"))
      # load("../../ASV_data/ASV_taxo/BAtaxo.Rda")
      if(nrow(BAtaxo_grouped_phylum)!=nrow(evalmet)){#securité
        message("problem nrow")
        
      }
      evalmet_BA<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Bacteria",],evalmet[BAtaxo_grouped_phylum$Kingdom=="Bacteria",])
      evalmet_BA$Phylum<-factor(evalmet_BA$Phylum,levels=c(names(table(evalmet_BA$Phylum[!(evalmet_BA$Phylum%in%c("Others","unclassified_Bacteria"))])),"Others","unclassified_Bacteria"))
      evalmet_AR<-cbind(BAtaxo_grouped_phylum[BAtaxo_grouped_phylum$Kingdom=="Archaea",],evalmet[BAtaxo_grouped_phylum$Kingdom=="Archaea",])
      evalmet_AR$Phylum<-factor(evalmet_AR$Phylum,levels=c(names(table(evalmet_AR$Phylum[!(evalmet_AR$Phylum%in%c("Others","unclassified_Archaea"))])),"Others","unclassified_Archaea"))
      
      # nrow(evalmet_BA)
      # nrow(taxo_grouped_phylum_BA)
      # mean(evalmet_BA$Dspear,na.rm=TRUE) #0.37
      # sd(evalmet_BA$Dspear,na.rm=TRUE) # 0.16
      # mean(evalmet_AR$Dspear,na.rm=TRUE) # 0.36
      # sd(evalmet_AR$Dspear,na.rm=TRUE) #0.16
      # nrow(evalmet_BA)
      # nrow(taxo_grouped_phylum_BA)
      
      # unique(taxo_grouped_phylum_BA$Species)
      #proportion of each phylum in good modeled ASV pool
    
      
      goodASV_phylum_02_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.2&!(is.na(evalmet_BA$Dspear))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.2&!(is.na(evalmet_BA$Dspear))]))
      goodASV_phylum_04_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.4&!(is.na(evalmet_BA$Dspear))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.4&!(is.na(evalmet_BA$Dspear))]))
      goodASV_phylum_05_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.5&!(is.na(evalmet_BA$Dspear))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.5&!(is.na(evalmet_BA$Dspear))]))
      goodASV_phylum_06_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.6&!(is.na(evalmet_BA$Dspear))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.6&!(is.na(evalmet_BA$Dspear))]))
      goodASV_phylum_08_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.8&!(is.na(evalmet_BA$Dspear))])/sum(summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.8&!(is.na(evalmet_BA$Dspear))]))
      
      
      goodASV_phylum_02_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.2&!(is.na(evalmet_AR$Dspear))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.2&!(is.na(evalmet_AR$Dspear))]))
      goodASV_phylum_04_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.4&!(is.na(evalmet_AR$Dspear))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.4&!(is.na(evalmet_AR$Dspear))]))
      goodASV_phylum_05_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.5&!(is.na(evalmet_AR$Dspear))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.5&!(is.na(evalmet_AR$Dspear))]))
      goodASV_phylum_06_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.6&!(is.na(evalmet_AR$Dspear))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.6&!(is.na(evalmet_AR$Dspear))]))
      goodASV_phylum_08_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.8&!(is.na(evalmet_AR$Dspear))])/sum(summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.8&!(is.na(evalmet_AR$Dspear))]))
      
      
      # number instead of proportions
      # summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.5])
      # Acidobacteriota      Actinobacteriota          Bacteroidota      Bdellovibrionota           Chloroflexi            Firmicutes 
      # 1887                   932                   508                   138                   489                   102 
      # Myxococcota       Patescibacteria       Planctomycetota        Proteobacteria     Verrucomicrobiota                Others 
      # 357                   201                   756                  2487                   715                   296 
      # unclassified_Bacteria                  NA's 
      #            1440
      # summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.5])
       # Crenarchaeota        Euryarchaeota        Halobacterota        Iainarchaeota               Others unclassified_Archaea 
       #            28                    0                    0                    0                    0                    0
      
      #proportion of each phylum in all modeled ASV pool
      AllASV_phylum_BA<-summary(factor(evalmet_BA$Phylum))/sum(summary(factor(evalmet_BA$Phylum)))
      AllASV_phylum_AR<-summary(factor(evalmet_AR$Phylum))/sum(summary(factor(evalmet_AR$Phylum)))
      
      comparison_BA<-melt(rbind(AllASV_phylum_BA,goodASV_phylum_02_BA,goodASV_phylum_04_BA,goodASV_phylum_06_BA,goodASV_phylum_08_BA))
      comparison_BA$value[is.na(comparison_BA$value)]<-0
      comparison_AR<-melt(rbind(AllASV_phylum_AR,goodASV_phylum_02_AR,goodASV_phylum_04_AR,goodASV_phylum_06_AR,goodASV_phylum_08_AR))
      comparison_AR$value[is.na(comparison_AR$value)]<-0
      
      #proportion of good modeled ASV in each phylum
      goodASV_02_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.2&!(is.na(evalmet_BA$Dspear))])/summary(factor(evalmet_BA$Phylum))
      goodASV_02_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.2&!(is.na(evalmet_AR$Dspear))])/summary(factor(evalmet_AR$Phylum))
      goodASV_04_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.4&!(is.na(evalmet_BA$Dspear))])/summary(factor(evalmet_BA$Phylum))
      goodASV_04_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.4&!(is.na(evalmet_AR$Dspear))])/summary(factor(evalmet_AR$Phylum))
      #proportion of verygood modeled ASV in each phylum
      #proportion of verygood modeled ASV in each phylum
      goodASV_05_BA<-summary(factor(evalmet_BA$Phylum)[evalmet_BA$Dspear>0.5&!(is.na(evalmet_BA$Dspear))])/summary(factor(evalmet_BA$Phylum))
      if(length(goodASV_05_BA)==0){
        goodASV_05_BA<-goodASV_02_BA*0
      }
      goodASV_05_AR<-summary(factor(evalmet_AR$Phylum)[evalmet_AR$Dspear>0.5&!(is.na(evalmet_AR$Dspear))])/summary(factor(evalmet_AR$Phylum))
      if(length(goodASV_05_AR)==0){
        goodASV_05_AR<-goodASV_02_AR*0
      }
      
      
      molten_good_BA <- melt(data.frame(goodASV_02=goodASV_02_BA,phylum=names(goodASV_02_BA),goodASV_04=goodASV_04_BA))
      molten_good_BA$phylum<-factor(molten_good_BA$phylum,levels=levels(evalmet_BA$Phylum))
      labels_BA <- paste(unique(molten_good_BA$phylum), " ( n =",c(table(evalmet_BA$Phylum)[-which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))],table(evalmet_BA$Phylum)[which(names(table(evalmet_BA$Phylum))%in%c("Others","unclassified_Bacteria"))]), ")")
      molten_good_AR <- melt(data.frame(goodASV_OK=goodASV_02_AR,phylum=names(goodASV_02_AR),goodASV_04=goodASV_04_AR))
      molten_good_AR$phylum<-factor(molten_good_AR$phylum,levels=levels(evalmet_AR$Phylum))
      labels_AR <- paste(unique(molten_good_AR$phylum), " ( n =",c(table(evalmet_AR$Phylum)[-which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))],table(evalmet_AR$Phylum)[which(names(table(evalmet_AR$Phylum))%in%c("Others","unclassified_Archaea"))]), ")")
      
      
      p1_BA<-ggplot(comparison_BA, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
        scale_fill_manual("Bacteria phylum",values = setcolors$BA)+ theme_bw() +
        scale_x_discrete(labels=c("AllASV_phylum_BA" = "Total", "goodASV_phylum_02_BA" = expression(rho~">0.2"), "goodASV_phylum_04_BA" = expression(rho~">0.4"), "goodASV_phylum_06_BA" = expression(rho~">0.6"), "goodASV_phylum_08_BA" = expression(rho~">0.8")))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.title = element_blank(),
              axis.title.x = element_blank())
      p2_BA<-ggplot(molten_good_BA,aes(x=variable,y=value,fill=phylum,group=phylum))+
        geom_bar(position="dodge",stat="identity")+
        ylim(0,1) + xlab("") + ylab("Proportion of models")+
        scale_fill_manual("Bacteria phylum",values = setcolors$BA,
                          label=labels_BA)+ 
        theme_bw() + scale_x_discrete(labels=c("goodASV_02_BA" = "rho>0.2", "goodASV_04" = "rho>0.4"))+ 
        theme_classic()+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=6),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0.5,0,0,0),"cm"),
              legend.title = element_text(size=7,colour = setcolors$dataset[group]),
              legend.text = element_text(size=5),
              legend.key.size = unit(0.3,"line"),
              legend.position = c(1,1.1),
              legend.justification = c("right", "top"),
              legend.background = element_blank())
      
      p1_AR<-ggplot(comparison_AR, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
        scale_fill_manual("Archaea phylum",values = setcolors$AR)+ theme_bw() +
        scale_x_discrete(labels=c("AllASV_phylum_AR" = "Total", "goodASV_phylum_02_AR" = expression(rho~">0.2"), "goodASV_phylum_04_AR" = expression(rho~">0.4"), "goodASV_phylum_06_AR" = expression(rho~">0.6"), "goodASV_phylum_08_AR" = expression(rho~">0.8")))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=10),
              legend.title = element_blank(),
              axis.title.x = element_blank())
      
      p2_AR<-ggplot(molten_good_AR,aes(x=variable,y=value,fill=phylum,group=phylum))+
        geom_bar(position="dodge",stat="identity")+
        ylim(0,1) + xlab("") + ylab("")+
        scale_fill_manual("Archaea phylum",values = setcolors$AR,
                          label=labels_AR)+ 
        theme_bw() + scale_x_discrete(labels=c("goodASV_02_BA" = "rho>0.2", "goodASV_04" = "rho>0.4"))+ 
        theme_classic()+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.title.y = element_text(size=10),
              axis.text.y = element_text(size=6),
              axis.text.x = element_blank(),
              plot.margin = unit(c(0.5,0,0,0),"cm"),
              legend.title = element_text(size=7,colour = setcolors$dataset["AR"]),
              legend.text = element_text(size=5),
              legend.key.size = unit(0.3,"line"),
              legend.position = c(1,1.1),
              legend.justification = c("right", "top"),
              legend.background = element_blank())

    } else{
      
      load(paste0("../../ASV_data/ASV_taxo/",group,"taxo_grouped_phylum.Rda"))
      taxo_grouped_phylum_temp<-eval(parse(text=paste0(group,"taxo_grouped_phylum")))
      #remove not fungi
      # table(taxo_grouped_phylum_temp$Phylum)
      
      if(nrow(taxo_grouped_phylum_temp)!=nrow(evalmet)){#securité
        message("problem")
        if(nrow(taxo_grouped_phylum_temp)!=nrow(evalmet)&nrow(evalmet)==3161){#securité
          # ncol(OTUdata_PR_AB) #3160
          # which(colnames(OTUdata_PR_AB)=="GCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTCGAAGCTAGAGGCACCGGGGCGCGGGGGTCACTGACCGTCTGCGCCTGCGGGCCTGAACCGTATTCCGGTTTGCGCGTGTGCTTTGTTGTGTGCGTGCGTGCCGGAACGATTACCTTGAGAAAATTAGAGTGTTCAAAGCAGGCAGTTGCTCGAATACATTAGCATGGAATAATAAAAGAGGACTCGGGTTCTTTTTTGTTGGTTTATAGGACCGAGTAATGATTAATAGGAACGGTCGGGGGCATTGGTATTGCTGTGTCAGGAGTGAAATCTGGTGACCATGGTGGGACCAACAAGTGCAAAGGCATTTGCCAAGGACGTTTCCAT")
          #annoying  Champi = 3019
          evalmet<-evalmet[-3019,]
          if(nrow(taxo_grouped_phylum_temp)==nrow(evalmet)){#securité
            message("problem solved")
          }
        }
      }
      
      
      if(group=="PR"){
        
        evalmet<-cbind(evalmet[taxo_grouped_phylum_temp$Phylum!="Not_Protist",],taxo_grouped_phylum_temp[taxo_grouped_phylum_temp$Phylum!="Not_Protist",])
        evalmet$Phylum<-factor(evalmet$Phylum,levels=c(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])[which(!(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])%in%c("Others")))],"Others"))
      }
      if(group=="FU"){
        evalmet<-cbind(evalmet[taxo_grouped_phylum_temp$Phylum!="Not_Fungi",],taxo_grouped_phylum_temp[taxo_grouped_phylum_temp$Phylum!="Not_Fungi",])
        evalmet$Phylum<-factor(evalmet$Phylum,levels=c(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])[which(!(names(table(evalmet$Phylum)[order(table(evalmet$Phylum),decreasing=TRUE)])%in%c("Others","Unidentified_Fungi")))],"Others","Unidentified_Fungi"))
      }
      # table(evalmet$Phylum)
      evalmet$Phylum<-droplevels(evalmet$Phylum)
      
      # mean(evalmet$Dspear,na.rm=TRUE) #FU 0.19  #PR  0.21
      # sd(evalmet$Dspear,na.rm=TRUE) #FU 0.20  #PR   0.13
      
      # #proportion of each phylum in good modeled ASV pool
      goodASV_phylum_02<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.2&!(is.na(evalmet$Dspear))])/sum(summary(factor(evalmet$Phylum)[evalmet$Dspear>0.2&!(is.na(evalmet$Dspear))]))
      # sum(goodASV_phylum_02)
      goodASV_phylum_04<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.4&!(is.na(evalmet$Dspear))])/sum(summary(factor(evalmet$Phylum)[evalmet$Dspear>0.4&!(is.na(evalmet$Dspear))]))
      goodASV_phylum_05<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.5&!(is.na(evalmet$Dspear))])/sum(summary(factor(evalmet$Phylum)[evalmet$Dspear>0.5&!(is.na(evalmet$Dspear))]))
      goodASV_phylum_06<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.6&!(is.na(evalmet$Dspear))])/sum(summary(factor(evalmet$Phylum)[evalmet$Dspear>0.6&!(is.na(evalmet$Dspear))]))
      goodASV_phylum_08<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.8&!(is.na(evalmet$Dspear))])/sum(summary(factor(evalmet$Phylum)[evalmet$Dspear>0.8&!(is.na(evalmet$Dspear))]))
      
      # number instead of proportions
      # summary(factor(evalmet$Phylum)[evalmet$Dspear>0.5])
      # Ascomycota      Basidiomycota      Glomeromycota  Mortierellomycota             Others Unidentified_Fungi     
      #          374                167                 61                194                  4                317 
      # protists
      # Cercozoa    Ciliophora Stramenopiles     Amoebozoa   Chlorophyta   Apicomplexa        Others         
      #      24             0            16             2             3             1             4 
      
      AllASV_phylum<-table(evalmet$Phylum)/sum(table(evalmet$Phylum))

      comparison<-melt(rbind(AllASV_phylum,goodASV_phylum_02,goodASV_phylum_04,goodASV_phylum_06,goodASV_phylum_08))
      comparison$value[is.na(comparison$value)]<-0

      #proportion of good modeled ASV in each phylum
      goodASV_02<-summary(factor(evalmet[evalmet$Dspear>0.2&!(is.na(evalmet$Dspear)),]$Phylum))/summary(factor(evalmet$Phylum))
      goodASV_04<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.4&!(is.na(evalmet$Dspear))])/summary(factor(evalmet$Phylum))
      #proportion of verygood modeled ASV in each phylum
      goodASV_05<-summary(factor(evalmet$Phylum)[evalmet$Dspear>0.5&!(is.na(evalmet$Dspear))])/summary(factor(evalmet$Phylum))
      if(length(goodASV_05)==0){
        goodASV_05<-goodASV_OK*0
      }
      
      
      molten_good <- melt(data.frame(goodASV_02=goodASV_02,phylum=names(goodASV_phylum_02),goodASV_04=goodASV_04))
      molten_good$phylum<-factor(molten_good$phylum,levels=levels(evalmet$Phylum))
      # labels <- paste(unique(molten_good$phylum), " ( n =",c(table(evalmet$Phylum)[-which(names(table(evalmet$Phylum))=="Others_Fungi")],table(evalmet$Phylum)[which(names(table(evalmet$Phylum))=="Others_Fungi")]), ")")
      labels <- paste(unique(molten_good$phylum), " (n =",table(evalmet$Phylum), ")",sep="")
      
      if(group=="PR"){plotmarg<-unit(c(0,0.5,0,0),"cm")}else{plotmarg<-unit(c(0,0,0,0.5),"cm")}
      assign(paste0("p1_", group), ggplot(comparison, aes(x=Var1, y=value, fill=Var2)) +
               geom_bar(stat="identity",show.legend = FALSE) + ylab("")+
               scale_fill_manual(paste0(ifelse(group=="PR","Protist phylum","Fungi phylum")),values = setcolors[[group]])+ theme_bw() +
               scale_x_discrete(labels=c("AllASV_phylum" = "Whole dataset", "goodASV_phylum_02" = expression(rho~">0.2"), "goodASV_phylum_04" = expression(rho~">0.4"), "goodASV_phylum_06" = expression(rho~">0.6"), "goodASV_phylum_08" = expression(rho~">0.8")))+ 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_text(size=10),
                     axis.text.y = element_text(size=10),
                     legend.title = element_blank(),
                     axis.title.x = element_blank())
      )
      assign(paste0("p2_", group), ggplot(molten_good,aes(x=variable,y=value,fill=phylum,group=phylum))+
               geom_bar(position="dodge",stat="identity")+
               ylim(0,1) + xlab("")  + ylab(ifelse(group=="PR","","Proportion of models"))+
               scale_fill_manual(paste0(ifelse(group=="PR","Protist phylum","Fungi phylum")),values = setcolors[[group]],
                                 label=labels)+ 
               theme_bw() + scale_x_discrete(labels=c("goodASV_02" = "rho > 0.2", "goodASV_04" = "rho > 0.4"))+ 
               theme_classic()+ 
               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.title.y = element_text(size=10),
                     axis.text.y = element_text(size=6),
                     axis.text.x = element_text(size=10),
                     plot.margin = unit(c(0,0,0,0),"cm"),
                     legend.title = element_text(size=7,colour = setcolors$dataset[group]),
                     legend.text = element_text(size=6),
                     legend.key.size = unit(0.3,"line"),
                     legend.position = c(.95,1),
                     legend.justification = c("right", "top"),
                     legend.background = element_blank())
      )
    }#p2_AR

  }
  pdf(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_proportion",Mod,"_",PAAB,".pdf"))
  grid.arrange(p2_BA,p2_AR,p2_FU,p2_PR,nrow=2)
  grid.arrange(p1_BA,p1_AR,p1_FU,p1_PR,nrow=4)
  dev.off()
  
  png(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_proportion",Mod,"_",PAAB,".png"),res=300,width=1961,height=1500)
  grid.arrange(p2_BA,p2_AR,p2_FU,p2_PR,nrow=2)
  dev.off()
  
  png(file=paste0("figures/PAAB_selection/test_phylum/goodmodels/Phylum_relative_richness",Mod,"_",PAAB,".png"),res=300,width=1961,height=1500)
  grid.arrange(p1_BA,p1_AR,p1_FU,p1_PR,nrow=4)
  dev.off()
}


