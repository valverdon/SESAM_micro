#figures Varimp

library(ggplot2)
library(tidyverse)
library(reshape2)
library(grid)
library(gridExtra)
library(paletteer)

PAAB="PA"
# PAAB="AB"
GtoM="BA"
algo="GLM"
# load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
#colnames(ENVdata)


climatic<-c("bio1_t", "bio10_twa", "bio11_tco", "bio12_p" ,"bio13_pwe","bio14_pdr","bio15_ps" ,"bio16_pwe","bio17_pdr","bio18_pwa","bio19_pco",
            "bio2_tdr" ,"bio3_tiso","bio4_ts" , "bio5_tmax","bio6_tmin","bio7_tar", "bio8_twet","bio9_tdry","GDD0" , "ETP" ,  "cumday_no","sRadY","noise","MoistVar" )
edaphic<-c("pH","pH.1","bulkSoilW","soilTemp","EC_1_5",    "TotalP",    "Nitrogen",  "Carbon",    "Hydrogen",  "Phyllosil", "Quartz",    "Feldspath","Plagiocla","MassiveLi",
           "Calcite",   "Indoses",   "SiO2",      "TiO2",     "Al2O3",     "Fe2O3" ,"MnO"  , "MgO",   "CaO",   "Na2O","Marlyshal","MarlShale",
           "K2O" ,  "P2O5"  , "OM"  ,  "Cr2O3",  "NiO",   "d15N" ,  "d13C","Silt_clay", "clay"  , "ThinSilt", "ThickSilt", "ThinSand", "ThickSand","Soil_aera","Soil_humu","Soil_mois","Soil_mois.1","Soil_nutr" )
topographic<-c("aspect","slope","Elevation","Altitude")
landcover<-c("forest_ag", "hydro_agg", "lowVeg_ag", "anthropos",   "deciduous")
remote.sensing<-c("ndvi","ndmi")
col_table<-data.frame(categ=c("climatic","edaphic","landcover","remote.sensing","topographic"),color=c("#F5E6A4",'#F5BBAE',"#8CF5CE",'#457EB0','gray69'))

Prop_models_one_cov<-matrix(NA,ncol=5,nrow=2*4*4)
colnames(Prop_models_one_cov)<-c("climatic","edaphic","landcover","remote.sensing","topographic")
rownames(Prop_models_one_cov)<-c("PA_GLM_PR","PA_GLM_FU","PA_GLM_BA","PA_GLM_AR","PA_GAM_PR","PA_GAM_FU","PA_GAM_BA","PA_GAM_AR","PA_GBM_PR","PA_GBM_FU","PA_GBM_BA","PA_GBM_AR","PA_RF_PR","PA_RF_FU","PA_RF_BA","PA_RF_AR",
                                 "AB_GLM_PR","AB_GLM_FU","AB_GLM_BA","AB_GLM_AR","AB_GAM_PR","AB_GAM_FU","AB_GAM_BA","AB_GAM_AR","AB_GBM_PR","AB_GBM_FU","AB_GBM_BA","AB_GBM_AR","AB_RF_PR","AB_RF_FU","AB_RF_BA","AB_RF_AR")

for (PAAB in c("PA","AB")){
  print(PAAB)
  for(algo in c("GLM","GAM","GBM","RF")){#Pb empty Variable_ranks.Rda GAM+GBM PA --> Bug Figure Selected Var
    print(algo)
    for(GtoM in c("PR","FU","BA")){
      print(GtoM)
      load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
      colnames(ENVdata)[colnames(ENVdata)=="Soil_moisture_variability"]<-"MoistVar"
      colnames(ENVdata)<-substr(colnames(ENVdata),1,9)
        ENVdataforcor<-ENVdata
        ENVdataforcor<-apply(ENVdataforcor, 2, function(x) as.numeric(as.character(x)))
        # varcor<-cor(ENVdataforcor, use = "pairwise.complete.obs")
        # varcor[abs(varcor[,"noise"])>0.7,"noise"]
        
        cluster<-hclust(as.dist(1-abs(cor(apply(ENVdata,2,as.numeric), use = "pairwise.complete.obs"))))
        cuttree<-cutree(cluster,h=0.3)
        # groups<-list()
        # for (i in 1:cuttree[length(cuttree)]){
        # groups[[i]]<-names(cuttree)[cuttree==i]
        # }
        # if(sum(unlist(lapply(groups,function(X){any(X=="pH"|X=="pH.1")})))>1){
        #   groups[[which(unlist(lapply(groups,function(X){any(X=="pH.1")})))]]<-NULL
        #   groups[[which(unlist(lapply(groups,function(X){any(X=="pH")})))]]<-c("pH","pH.1")
        # }
        # if(sum(unlist(lapply(groups,function(X){"ndvi"%in%X&"pH"%in%X})))==1){
        #   groups[[which(unlist(lapply(groups,function(X){"ndvi"%in%X&"pH"%in%X})))]]<-c("pH","ndvi","pH.1")
        # }
        #manual reassignation of clusters to accomodate all 3 groups with the same cluster
        groups<-list(c("bio1_t", "bio10_twa", "bio11_tco", "bio12_p",   "bio13_pwe", "bio14_pdr", "bio16_pwe", "bio17_pdr", "bio18_pwa", "bio19_pco" ,"bio4_ts" , 
                       "bio5_tmax" ,"bio6_tmin" ,"bio7_tar",  "bio8_twet" ,"bio9_tdry", "GDD0", "ETP", "cumday_no", "Soil_nutr", "Elevation", "noise","d13C"),#T°
                     c("bio2_tdr","bio3_tiso"),#day_var_temp
                     "bio15_ps",
                     c("sRadY" , "aspect"),#Sol_rad
                     "slope",
                     c("pH" ,  "ndvi" ,"pH.1","Soil_aera" ,"Soil_humu", "MoistVar"),#pH
                     "Soil_mois",
                     "MarlShale",
                     "Marlyshal",
                     "MassiveLi",
                     "Silt_clay",
                     "forest_ag",
                     "hydro_agg",
                     "lowVeg_ag",
                     "anthropos",
                     "deciduous",
                     "ndmi",
                     "soilTemp",
                     c("bulkSoilW" ,"Nitrogen"  ,"Hydrogen","EC_1_5", "Carbon", "OM"  ),#C/N
                     c("TotalP" ,"P2O5"),#Phosphorus
                     "Phyllosil",
                     c("Quartz" ,"SiO2"),#Silicate
                     "Feldspath",
                     "Plagiocla",
                     c("Calcite" ,"CaO"),#Calcite
                     "Indoses",
                     c("TiO2" , "Al2O3", "Fe2O3", "K2O"),#Oxides
                     "MnO",
                     "MgO",
                     "Na2O",
                     "Cr2O3",
                     "NiO",
                     "d15N",
                     "clay",
                     c("ThinSilt", "ThinSand"),#granulo
                     "ThickSand",
                     "ThickSilt")
        # plot(cutree(cluster,h=0.2),cex=0.5)
        # 
        # png(file=paste0("figures/PAAB_selection/VARimp/corr",GtoM,"_",PAAB,".png"),res=300,width=1961,height=1500)
        # 
        # plot(cluster,cex=0.5)
        # abline(h=0.2,col="red")
        # dev.off()

      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_importance.Rda"))
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_ranks.Rda"))
      load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_preselected.Rda"))
      load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo_grouped_phylum.Rda"))
      Taxo<-eval(parse(text=paste0(GtoM,"taxo_grouped_phylum")))
      #checkings

        Variable_importance<- cbind(Variable_importance[Variable_importance$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_importance$seq,])
        all(Variable_importance$Seq==Variable_importance$seq)
        colnames(Variable_importance)[colnames(Variable_importance)=="Soil_moisture_variability"]<-"MoistVar"
        Variable_preselected<- cbind(Variable_preselected[Variable_preselected$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_preselected$seq,])
        all(Variable_preselected$Seq==Variable_preselected$seq)
        colnames(Variable_preselected)[colnames(Variable_preselected)=="Soil_moisture_variability"]<-"MoistVar"
        Variable_ranks<- cbind(Variable_ranks[Variable_ranks$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_ranks$seq,])
        all(Variable_importance$Seq==Variable_importance$seq)
        colnames(Variable_ranks)[colnames(Variable_ranks)=="Soil_moisture_variability"]<-"MoistVar"
        
      names(Variable_importance)<-substr(names(Variable_preselected),1,9)
      names(Variable_ranks)<-substr(names(Variable_ranks),1,9)
      names(Variable_preselected)<-substr(names(Variable_preselected),1,9)


      # summary(Variable_importance)
      # summary(Variable_ranks)
      # apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],1,function(X){sum(!is.na(X))})
      
      # summary(Variable_preselected)
      if(GtoM=="BA"){
        
        Variable_importance_BA <- Variable_importance[Variable_importance$Kingdom=="Bacteria",]
        # Variable_importance_BA <-Variable_importance_BA[apply(Variable_importance_BA[,-which(colnames(Variable_importance_BA)=="seq")],1,function(X){sum(X)>0}),]
        Variable_importance_AR <- Variable_importance[Variable_importance$Kingdom=="Archaea",]
        # Variable_importance_AR <-Variable_importance_AR[apply(Variable_importance_AR[,-which(colnames(Variable_importance_AR)=="seq")],1,function(X){sum(X)>0}),]
        
        Variable_ranks_BA <- Variable_ranks[Variable_ranks$Kingdom=="Bacteria",]
        Variable_ranks_AR <- Variable_ranks[Variable_ranks$Kingdom=="Archaea",]
        Variable_preselected_BA <- Variable_preselected[Taxo$Kingdom=="Bacteria",]
        Variable_preselected_AR <- Variable_preselected[Taxo$Kingdom=="Archaea",]

        
#BA
        
        allpresel<-apply(Variable_preselected_BA[,-(which(colnames(Variable_preselected_BA)=="seq"):ncol(Variable_preselected_BA))],2,function(X){sum(X)/length(X)})
        allpresel2<-allpresel[order(allpresel,decreasing=TRUE)]

        #per variable, proportion of models in which they where preselected
        # barplot(allpresel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Preselection"),las = 2)
        #ggplot
        allpresel3<-data.frame(value=allpresel2,
                               variable=factor(names(allpresel2),levels=names(allpresel2)),
                               categ=ifelse(names(allpresel2)%in%climatic,"climatic",
                                            ifelse(names(allpresel2)%in%edaphic,"edaphic",
                                                   ifelse(names(allpresel2)%in%topographic,"topographic",
                                                          ifelse(names(allpresel2)%in%landcover,"landcover",
                                                                 ifelse(names(allpresel2)%in%remote.sensing,"remote.sensing",names(allpresel2)))))))
        allpresel3$categ<-factor(allpresel3$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        categinplot_allpresel3<-unique(allpresel3$categ)   
        allpresel3_10<-allpresel3[1:10,]
        categinplot_allpresel3_10<-unique(allpresel3[1:10,]$categ)   
        # ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
        #   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"))+ ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))
        assign(paste0("pPres_BA_",algo,"_",PAAB),ggplot(allpresel3, aes(x=variable,y=value,fill=categ))+theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allpresel3])+
                 scale_y_continuous(limits=c(0,1)))


      #pPres_BA_GLM_PA
        #per variable, proportion of models in which they where selected
        allsel_BA<-apply(Variable_ranks_BA[,-(which(colnames(Variable_ranks)=="seq"):ncol(Variable_ranks_BA))],2,function(X){ifelse(is.na(X),NA,1)})
        # summary(allsel_BA)
        allsel_BA2<-matrix(NA,ncol=length(groups),nrow=nrow(allsel_BA))
        rownames(allsel_BA2)<-rownames(allsel_BA)
        colnames(allsel_BA2)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=1
          if(length(groups[[i]])>1){
            allsel_BA2[,groups[[i]][1]]<-apply(allsel_BA[,colnames(allsel_BA)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE)!=0)})
          } else{
            allsel_BA2[,groups[[i]][1]]<-allsel_BA[,colnames(allsel_BA)%in%groups[[i]]]
          }
        }
          
        #proportion of models in which variable selected
        allsel2<-apply(allsel_BA2,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        # barplot(allsel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Selection"),las = 2)
        allsel2<-allsel2[order(allsel2,decreasing=TRUE)]
        #add categories
        allsel3<-data.frame(value=allsel2,
                               variable=factor(names(allsel2),levels=names(allsel2)),
                               categ=ifelse(names(allsel2)%in%climatic,"climatic",
                                            ifelse(names(allsel2)%in%edaphic,"edaphic",
                                                   ifelse(names(allsel2)%in%topographic,"topographic",
                                                          ifelse(names(allsel2)%in%landcover,"landcover",
                                                                 ifelse(names(allsel2)%in%remote.sensing,"remote.sensing",names(allsel2)))))))
        allsel3$categ<-factor(allsel3$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
                # allsel3$variable<-as.character(allsel3$variable)
        # if(allsel3$variable[allsel3$variable=="Altitude"|allsel3$variable=="Elevation"]!=colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]){
        #   allsel3$variable[which(allsel3$variable=="Altitude"|allsel3$variable=="Elevation")]<-colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]
        # }
        # allsel3$variable<-as.factor(allsel3$variable)

        #MERGE categ if same group of correlated vars
        # categ<-vector()
        # for( i in 1:length(allsel3$variable)){
        #   vartest<-allsel3$variable[i]
        #   categ[i]<-paste0(groups[[c(which(unlist(lapply(groups,function(X){vartest%in%X}))))]])
        #   }
        # 
        # nrow(allsel3)
        allsel3$variable<-recode_factor(allsel3$variable,bio1_t="AirTemp",bio2_tdr="day_var_temp",sRadY="solar_rad",bulkSoilW="C_N",TotalP="Phosphorus",Quartz="Silicate",TiO2="Oxides",ThinSilt="granulo")
        allsel3$variable<-with(allsel3,reorder(variable,value,decreasing=TRUE))
        if(sum(as.character(allsel3$variable)=="pH")>1){
          allsel3$variable<-factor(allsel3$variable,levels=as.character(allsel3$variable[-length(allsel3$variable)]))
        } else{allsel3$variable<-factor(allsel3$variable,levels=as.character(allsel3$variable))}

        categinplot_allsel3<-unique(allsel3$categ)
        allsel3_10<-allsel3[1:10,]
        print(allsel3_10)
        levels(allsel3$variable)
        categinplot_allsel3_10<-unique(allsel3[1:10,]$categ)   
        assign(paste0("pSel_BA_",algo,"_",PAAB),ggplot(allsel3, aes(x=variable,y=value,fill=categ))+theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
          scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3])+
            scale_y_continuous(limits=c(0,1)))
        assign(paste0("pSel10_BA_",algo,"_",PAAB),ggplot(allsel3_10, aes(x=variable,y=value,fill=categ))+theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3_10])+
                 scale_y_continuous(limits=c(0,1)))
      #pSel10_BA_GLM_PA
        
        
       #    #list of ranks obtained per variable
       #    allranks<-apply(Variable_ranks_BA[,-which(colnames(Variable_ranks_BA)=="seq")],2,function(X){na.exclude(X)})
       #    # summary(Variable_ranks_BA)
       #    Var_ranks<-lapply(allranks,as.vector)
       #    #per variable, it's mean rank whithin postselected variables
       #    reordering<-with(melt(Var_ranks),                       # Order boxes by median
       #                     reorder(L1,
       #                             value,
       #                             median))
       # 
       #    data_reordered<-melt(Var_ranks)
       #    data_reordered$L1<-factor(data_reordered$L1,
       #                              levels = levels(reordering))
       #    data_reordered$categ<-unlist(lapply(as.vector(data_reordered$L1),function(X){ifelse(X%in%climatic,"climatic",
       #                                                                                        ifelse(X%in%edaphic,"edaphic",
       #                                                                                               ifelse(X%in%topographic,"topographic",
       #                                                                                                      ifelse(X%in%remote.sensing,"remote.sensing",
       #                                                                                                             ifelse(X%in%landcover,"landcover",X)))))}))
       # 
       # categinplot_pRank<-unique(data_reordered$categ)
       # # summary(data_reordered$value)
       # # ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+geom_boxplot()+theme_bw()
       # assign(paste0("pRank_BA_",algo,"_",PAAB),ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+geom_boxplot()+theme_bw()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
       #          theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
       #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       #                legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+
       #          ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
       #          scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_pRank]))
        
      #pRank_BA_GLM_PA
        
        #variable importance corrected by the maximal value of importance across variables (best =1)
        Variable_importance1<-t(apply(Variable_importance_BA[,-(which(colnames(Variable_importance_BA)=="seq"):ncol(Variable_importance_BA))],1,function(X){X/max(X)}))
       summary(Variable_importance_BA);summary(Variable_importance1)
       #regroup variable correlated (sum of all)
        Variable_importance2<-matrix(NA,ncol=length(groups),nrow=nrow(Variable_importance1))
        rownames(Variable_importance2)<-rownames(Variable_importance1)
        colnames(Variable_importance2)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=1
          if(length(groups[[i]])>1){
            Variable_importance2[,groups[[i]][1]]<-apply(Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE))})
          } else{
            Variable_importance2[,groups[[i]][1]]<-ifelse(is.na(Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]]),NA,Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]])
          }
        }
        # summary(Variable_importance2);summary(Variable_importance1)
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        #keep value only if variable is selected in model
        
        allsel4<-apply(Variable_importance2*allsel_BA2,2,function(X){na.exclude(X)})
        # summary(allsel_BA2)
        # allsel4<-Variable_importance2[Variable_importance2==0]
        
        # summary(allsel4)
        # boxplot(allsel4,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_only_when_selected"),las = 2)

        # levels(allsel3$variable)

        data_reordered4<-melt(allsel4)
        data_reordered4$categ<-unlist(lapply(as.vector(data_reordered4$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                            ifelse(X%in%edaphic,"edaphic",
                                                                                                   ifelse(X%in%topographic,"topographic",
                                                                                                          ifelse(X%in%landcover,"landcover",
                                                                                                                 ifelse(X%in%remote.sensing,"remote.sensing",X)))))}))
        data_reordered4$L1<-recode_factor(data_reordered4$L1,bio1_t="AirTemp",bio2_tdr="day_var_temp",sRadY="solar_rad",bulkSoilW="C_N",TotalP="Phosphorus",Quartz="Silicate",TiO2="Oxides",ThinSilt="granulo")
        reordering4<-with(data_reordered4,                       # Order boxes by median
                          reorder(L1,
                                  value,
                                  median,decreasing=TRUE))
        data_reordered4$L1<-factor(data_reordered4$L1,
                                   levels = levels(reordering4))

        data_reordered4$categ<-factor(data_reordered4$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        categinplot_VIMP<-unique(data_reordered4$categ)
        data_reordered4_10<-data_reordered4[data_reordered4$L1%in%levels(data_reordered4$L1)[1:10],]
        categinplot_VIMP_10<-unique(data_reordered4_10$categ)
        # summary(data_reordered4[data_reordered4$L1=="Carbon",]$value)
        assign(paste0("pVIMP_BA_",algo,"_",PAAB),ggplot(data_reordered4,aes(x=L1,y=value,fill=categ))+theme_bw()+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pVIMP10_BA_",algo,"_",PAAB),ggplot(data_reordered4_10,aes(x=L1,y=value,fill=categ))+theme_bw()+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP_10])+
                 scale_y_continuous(limits=c(0,1)))
      #pVIMP_BA_GLM_PA
        
        allsel5_temp<-apply(allsel_BA2,1,function(X){c(any(X[names(X)%in%climatic]>=1),any(X[names(X)%in%edaphic]>=1),any(X[names(X)%in%landcover]>=1), any(X[names(X)=="ndmi"]>=1), any(X[names(X)%in%topographic]>=1))})
        allsel5_temp<-t(allsel5_temp)
        colnames(allsel5_temp)<-c("climatic","edaphic","landcover","remote.sensing","topographic")
        allsel5_temp<-apply(allsel5_temp,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        Prop_models_one_cov[paste0(PAAB,"_",algo,"_BA"),]<-allsel5_temp
        # allsel5_temp<-allsel5_temp[order(allsel5_temp,decreasing=TRUE)]
        allsel5<-data.frame(value=allsel5_temp,grouping=factor(names(allsel5_temp),levels=names(allsel5_temp)))
        print(allsel5)
        # assign(paste0("pgrouped_BA_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping))+ theme_bw() + 
        #          geom_bar(stat='identity', show.legend = FALSE) + 
        #          xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria")) + ylim(c(0,1)) +
        #          theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        #          scale_fill_manual(values=col_table$color[order(match(col_table$categ,levels(allsel5$grouping)))]))
        #     pgrouped_BA_GLM_PA    
        
        # assign(paste0("pgrouped_BA_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping))+ theme_bw() + 
        #          geom_bar(stat='identity', show.legend = FALSE) + 
        #          ylab(ifelse(PAAB=="AB","","Bacteria")) + ylim(c(0,1)) +
        #          theme(axis.text.x=element_text(angle=15,vjust=0.9,hjust=0.6,size=8),plot.margin=unit(c(0,5.5,0,5.5),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #                axis.title.x=element_blank()) +
        # scale_fill_manual(values=rep("grey",length(allsel5$grouping))))
        
        assign(paste0("pgrouped_BA_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping))+ theme_bw() + 
                 geom_bar(stat='identity', show.legend = FALSE) + 
                 ylab(ifelse(PAAB=="AB","","Bacteria")) + ylim(c(0,1)) +
                 theme(axis.text.x=element_blank(),plot.margin=unit(c(2,4,2,4),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank()) +
                 scale_fill_manual(values=rep("grey",length(allsel5$grouping))))
        #AR
        
        allpresel<-apply(Variable_preselected_AR[,-(which(colnames(Variable_preselected_AR)=="seq"):ncol(Variable_preselected_AR))],2,function(X){sum(X)/length(X)})
        allpresel2<-allpresel[order(allpresel,decreasing=TRUE)]
        #list of ranks obtained per variable
        allranks<-apply(Variable_ranks_AR[,-(which(colnames(Variable_ranks_AR)=="seq"):ncol(Variable_ranks_AR))],2,function(X){na.exclude(X)})
        Var_ranks<-lapply(allranks,as.vector)
        #per variable, proportion of models in which they where preselected
        # barplot(allpresel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Preselection"),las = 2)
        #ggplot
        allpresel3<-data.frame(value=allpresel2,variable=factor(names(allpresel2),levels=names(allpresel2)),categ=ifelse(names(allpresel2)%in%climatic,"climatic",
                                                                                                                         ifelse(names(allpresel2)%in%edaphic,"edaphic",
                                                                                                                                ifelse(names(allpresel2)%in%topographic,"topographic",
                                                                                                                                       ifelse(names(allpresel2)%in%landcover,"landcover",
                                                                                                                                              ifelse(names(allpresel2)%in%remote.sensing,"remote.sensing",names(allpresel2)))))))
        allpresel3$categ<-factor(allpresel3$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        
        # allsel3$variable<-as.character(allsel3$variable)
        # if(allsel3$variable[allsel3$variable=="Altitude"|allsel3$variable=="Elevation"]!=colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]){
        #   allsel3$variable[which(allsel3$variable=="Altitude"|allsel3$variable=="Elevation")]<-colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]
        # }
        # allsel3$variable<-as.factor(allsel3$variable)
        #MERGE categ if same group of correlated vars
        categinplot_allpresel3<-unique(allpresel3$categ)
        allpresel3_10<-allpresel3[1:10,]
        categinplot_allpresel3_10<-unique(allpresel3[1:10,]$categ)   
        # ggplot(allpresel3, aes(x=variable,y=value,fill=categ)) + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria"))  +
        #   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"))+ ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))
        assign(paste0("pPres_AR_",algo,"_",PAAB),ggplot(allpresel3, aes(x=variable,y=value,fill=categ))+theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allpresel3])+
                 scale_y_continuous(limits=c(0,1)))
        
      #pPres_AR_GLM_PA
        #per variable, proportion of models in which they where selected
        allsel_AR<-apply(Variable_ranks_AR[,-(which(colnames(Variable_ranks)=="seq"):ncol(Variable_ranks_AR))],2,function(X){ifelse(is.na(X),NA,1)})
        allsel_AR2<-matrix(NA,ncol=length(groups),nrow=nrow(allsel_AR))
        rownames(allsel_AR2)<-rownames(allsel_AR)
        colnames(allsel_AR2)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=1
          if(length(groups[[i]])>1){
            allsel_AR2[,groups[[i]][1]]<-apply(allsel_AR[,colnames(allsel_AR)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE)!=0)})
          }else{
            allsel_AR2[,groups[[i]][1]]<-allsel_AR[,colnames(allsel_AR)%in%groups[[i]]]
          }
        }
        
        allsel2<-apply(allsel_AR2,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        # barplot(allsel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Selection"),las = 2)
        allsel2<-allsel2[order(allsel2,decreasing=TRUE)]
        allsel3<-data.frame(value=allsel2,variable=factor(names(allsel2),levels=names(allsel2)),categ=ifelse(names(allsel2)%in%climatic,"climatic",
                                                                                                             ifelse(names(allsel2)%in%edaphic,"edaphic",
                                                                                                                    ifelse(names(allsel2)%in%topographic,"topographic",
                                                                                                                           ifelse(names(allsel2)%in%landcover,"landcover",
                                                                                                                                  ifelse(names(allsel2)%in%remote.sensing,"remote.sensing",names(allsel2)))))))
        allsel3$categ<-factor(allsel3$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        # allsel3$variable<-as.character(allsel3$variable)
        # if(allsel3$variable[allsel3$variable=="Altitude"|allsel3$variable=="Elevation"]!=colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]){
        #   allsel3$variable[which(allsel3$variable=="Altitude"|allsel3$variable=="Elevation")]<-colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]
        # }
        # allsel3$variable<-as.factor(allsel3$variable)
        #MERGE categ if same group of correlated vars
        # categ<-vector()
        # for( i in 1:length(allsel3$variable)){
        #   vartest<-allsel3$variable[i]
        #   categ[i]<-paste0(groups[[c(which(unlist(lapply(groups,function(X){vartest%in%X}))))]])
        # }
        # 
        # allsel3$categ<-categ
        allsel3$variable<-recode_factor(allsel3$variable,bio1_t="AirTemp",bio2_tdr="day_var_temp",sRadY="solar_rad",bulkSoilW="C_N",TotalP="Phosphorus",Quartz="Silicate",TiO2="Oxides",ThinSilt="granulo")
        allsel3$variable<-with(allsel3,reorder(variable,value,decreasing=TRUE))
        if(sum(as.character(allsel3$variable)=="pH")>1){
        allsel3$variable<-factor(allsel3$variable,levels=as.character(allsel3$variable[-length(allsel3$variable)]))
        } else{allsel3$variable<-factor(allsel3$variable,levels=as.character(allsel3$variable))}
        # levels(allsel3$variable)
        categinplot_allsel3<-unique(allsel3$categ) 
        allsel3_10<-allsel3[1:10,]
        print(allsel3_10)
        categinplot_allsel3_10<-unique(allsel3[1:10,]$categ)   
        assign(paste0("pSel_AR_",algo,"_",PAAB),ggplot(allsel3, aes(x=variable,y=value,fill=categ))+theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pSel10_AR_",algo,"_",PAAB),ggplot(allsel3_10, aes(x=variable,y=value,fill=categ))+theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3_10])+
                 scale_y_continuous(limits=c(0,1)))
        # pSel10_AR_GLM_PA
        #pH 83% of archaea models
        
        # #list of ranks obtained per variable
        # allranks<-apply(Variable_ranks_AR[,-which(colnames(Variable_ranks_AR)=="seq")],2,function(X){na.exclude(X)})
        # # summary(Variable_ranks_BA)
        # Var_ranks<-lapply(allranks,as.vector)
        # #per variable, it's mean rank whithin postselected variables
        # reordering<-with(melt(Var_ranks),                       # Order boxes by median
        #                  reorder(L1,
        #                          value,
        #                          median))
        # 
        # data_reordered<-melt(Var_ranks)
        # data_reordered$L1<-factor(data_reordered$L1,
        #                           levels = levels(reordering))
        # data_reordered$categ<-unlist(lapply(as.vector(data_reordered$L1),function(X){ifelse(X%in%climatic,"climatic",
        #                                                                                     ifelse(X%in%edaphic,"edaphic",
        #                                                                                            ifelse(X%in%topographic,"topographic",
        #                                                                                                   ifelse(X%in%remote_sensing,"remote_sensing",
        #                                                                                                   ifelse(X%in%landcover,"landcover",X)))))}))
        # 
        # categinplot_pRank<-unique(data_reordered$categ)   
        # assign(paste0("pRank_AR_",algo,"_",PAAB),ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+theme_bw() +geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
        #          theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #                legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
        #          ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
        #          scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_pRank])+
        #          scale_y_continuous(limits=c(0,1)))
        # 
        # #pRank_AR_GLM_PA
        
        #variable importance corrected by the maximal value of importance across variables (best =1)
        Variable_importance1<-t(apply(Variable_importance_AR[,-(which(colnames(Variable_importance_AR)=="seq"):ncol(Variable_importance_AR))],1,function(X){X/max(X)}))
        
        Variable_importance2<-matrix(0,ncol=length(groups),nrow=nrow(Variable_importance1))
        rownames(Variable_importance2)<-rownames(Variable_importance1)
        colnames(Variable_importance2)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=1
          if(length(groups[[i]])>1){
            Variable_importance2[,groups[[i]][1]]<-apply(Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE))})
          } else{
            Variable_importance2[,groups[[i]][1]]<-ifelse(is.na(Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]]),NA,Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]])
          }
        }
        Variable_importance2<-matrix(NA,ncol=length(groups),nrow=nrow(Variable_importance1))
        rownames(Variable_importance2)<-rownames(Variable_importance1)
        colnames(Variable_importance2)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=1
          if(length(groups[[i]])>1){
            Variable_importance2[,groups[[i]][1]]<-apply(Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE))})
          } else{
            Variable_importance2[,groups[[i]][1]]<-Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]]
          }
        }
        # summary(Variable_importance2);summary(Variable_importance1)
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        #keep value only if variable is selected in model
        
        
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        allsel4<-apply(Variable_importance2*allsel_AR2,2,function(X){na.exclude(X)})
        # boxplot(allsel4,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_only_when_selected"),las = 2)
        # summary(allsel_BA2)
        # allsel4<-Variable_importance2[Variable_importance2==0]
        
        # summary(allsel4)
        # boxplot(allsel4,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_only_when_selected"),las = 2)
        
        # levels(allsel3$variable)
        
        data_reordered4<-melt(allsel4)
        data_reordered4$categ<-unlist(lapply(as.vector(data_reordered4$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                              ifelse(X%in%edaphic,"edaphic",
                                                                                                     ifelse(X%in%topographic,"topographic",
                                                                                                                   ifelse(X%in%landcover,"landcover",
                                                                                                                          ifelse(X%in%remote.sensing,"remote.sensing",X)))))}))
        data_reordered4$L1<-recode_factor(data_reordered4$L1,bio1_t="AirTemp",bio2_tdr="day_var_temp",sRadY="solar_rad",bulkSoilW="C_N",TotalP="Phosphorus",Quartz="Silicate",TiO2="Oxides",ThinSilt="granulo")
        reordering4<-with(data_reordered4,                       # Order boxes by median
                          reorder(L1,
                                  value,
                                  median,decreasing=TRUE))
        data_reordered4$L1<-factor(data_reordered4$L1,
                                   levels = levels(reordering4))
        data_reordered4$categ<-factor(data_reordered4$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        
        categinplot_VIMP<-unique(data_reordered4$categ)
        data_reordered4_10<-data_reordered4[data_reordered4$L1%in%levels(data_reordered4$L1)[1:10],]
        categinplot_VIMP_10<-unique(data_reordered4_10$categ)

        # data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        assign(paste0("pVIMP_AR_",algo,"_",PAAB),ggplot(data_reordered4,aes(x=L1,y=value,fill=categ))+geom_boxplot()+theme_bw()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pVIMP10_AR_",algo,"_",PAAB),ggplot(data_reordered4_10,aes(x=L1,y=value,fill=categ))+theme_bw() +geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+
                 ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP_10])+
                 scale_y_continuous(limits=c(0,1)))
      #pVIMP_AR_GLM_PA
        allsel5_temp<-apply(allsel_AR2,1,function(X){c(any(X[names(X)%in%climatic]>=1),any(X[names(X)%in%edaphic]>=1),any(X[names(X)%in%landcover]>=1),any(X[names(X)%in%remote.sensing]>=1), any(X[names(X)%in%topographic]>=1))})
        allsel5_temp<-t(allsel5_temp)
        colnames(allsel5_temp)<-c("climatic","edaphic","landcover","remote.sensing","topographic")
        allsel5_temp<-apply(allsel5_temp,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        Prop_models_one_cov[paste0(PAAB,"_",algo,"_AR"),]<-allsel5_temp
        # allsel5_temp<-allsel5_temp[order(allsel5_temp,decreasing=TRUE)]
        allsel5<-data.frame(value=allsel5_temp,grouping=factor(names(allsel5_temp),levels=names(allsel5_temp)))
        print(allsel5)
        # assign(paste0("pgrouped_AR_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping))+theme_bw() + 
        #          geom_bar(stat='identity', show.legend = FALSE) + 
        #          xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  + ylim(c(0,1))+
        #          theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        #          ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
        #          scale_fill_manual(values=col_table$color[order(match(col_table$categ,levels(allsel5$grouping)))]))
        #     pgrouped_AR_GLM_PA      
        # assign(paste0("pgrouped_AR_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping))+theme_bw() + 
        #          geom_bar(stat='identity', show.legend = FALSE) + 
        #          ylab(ifelse(PAAB=="AB","","Archaea"))  + ylim(c(0,1))+
        #          theme(axis.text.x=element_text(angle=15,vjust=0.9,hjust=0.6,size=8),plot.margin=unit(c(0,5.5,0,5.5),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #                axis.title.x=element_blank()) +
        #          ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))+
        # scale_fill_manual(values=rep("grey",length(allsel5$grouping))))
        
        assign(paste0("pgrouped_AR_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping))+theme_bw() + 
                 geom_bar(stat='identity', show.legend = FALSE) + 
                 ylab(ifelse(PAAB=="AB","","Archaea"))  + ylim(c(0,1))+
                 theme(axis.text.x=element_blank(),plot.margin=unit(c(0,4,2,4),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank()) +
                 ggtitle(ifelse(PAAB=="AB","Relative abundance models","Presence-absence models"))+
                 scale_fill_manual(values=rep("grey",length(allsel5$grouping))))
        # pgrouped_AR_GLM_PA
        #variable importance corrected by model performance (attribute the full performance to all variables according to their importance)
          # if(PAAB=="PA"){
          #   #variable importance corrected by model performance (attribute the full performance to all variables according to their importance)
          #     load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval.Rda"))
          #     Eval_BA<-Eval[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
          #     Eval_AR<-Eval[BAtaxo$Kingdom=="Bacteria"&!is.na(BAtaxo$Kingdom),]
          #     TSSadj<-ifelse(Eval_BA$TSS_adj>0,Eval_BA$TSS_adj,NA)
          #     varimp_percentage<-t(apply(Variable_importance_BA[,-which(colnames(Variable_ranks_BA)=="seq")],1,function(X){X/sum(X,na.rm=TRUE)}))
          #     varimp_percentage<-varimp_percentage*allsel_BA
          #     varimp_percentage_weighted<-apply(varimp_percentage,2,function(X){X*TSSadj})
          # 
          #     # boxplot(varimp_percentage_weighted,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_weighted_by_ModPerf"),las = 2)
          #     reordering5<-with(melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),],                       # Order boxes by median
          #                       reorder(Var2,
          #                               value,
          #                               median,decreasing=TRUE))
          #     
          #     data_reordered5<-melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),]
          #     data_reordered5$Var2<-factor(data_reordered5$Var2,
          #                                  levels = levels(reordering5))
          #     assign(paste0("pVIMPw_BA_",algo,"_",PAAB),ggplot(data_reordered5,aes(x=Var2,y=value))+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","","Bacteria")) + 
          #              theme_bw() +
          #              theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
          #                    panel.grid.major = element_blank(), panel.grid.minor = element_blank()))+ ggtitle(ifelse(PAAB=="AB","Rel. Abundance models","Presence-Absence models"))
          #     #pVIMPw_BA_GLM_PA
          #     
          #   TSSadj<-ifelse(Eval_AR$TSS_adj>0,Eval_AR$TSS_adj,NA)
          #   varimp_percentage<-t(apply(Variable_importance_AR[,-which(colnames(Variable_ranks_AR)=="seq")],1,function(X){X/sum(X,na.rm=TRUE)}))
          #   varimp_percentage<-varimp_percentage*allsel_AR
          #   varimp_percentage_weighted<-apply(varimp_percentage,2,function(X){X*TSSadj})
          #   
          #   # boxplot(varimp_percentage_weighted,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_weighted_by_ModPerf"),las = 2)
          #   reordering5<-with(melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),],                       # Order boxes by median
          #                     reorder(Var2,
          #                             value,
          #                             median,decreasing=TRUE))
          #   
          #   data_reordered5<-melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),]
          #   data_reordered5$Var2<-factor(data_reordered5$Var2,
          #                                levels = levels(reordering5))
          #   assign(paste0("pVIMPw_AR_",algo,"_",PAAB),ggplot(data_reordered5,aes(x=Var2,y=value))+geom_boxplot()+ theme_bw()+ xlab("") + ylab(ifelse(PAAB=="AB","","Archaea"))  +
          #            theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
          #                  panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
          # #pVIMPw_AR
          # }
      }
      if(GtoM %in% c("PR","FU")){ #GtoM="FU"
#FU &PR
        if(GtoM=="FU"){
          Variable_preselected<-Variable_preselected[Variable_preselected$Kingdom=="Fungi",]
          Variable_importance<-Variable_importance[Variable_importance$Kingdom=="Fungi",]
          Variable_ranks<-Variable_ranks[Variable_ranks$Kingdom=="Fungi",]
        }
        if(GtoM=="PR"){
          Variable_preselected<-Variable_preselected[Variable_preselected$Phylum!="Not_Protist",]
          Variable_importance<-Variable_importance[Variable_importance$Kingdom!="Not_Protist",]
          Variable_ranks<-Variable_ranks[Variable_ranks$Kingdom!="Not_Protist",]
        }
        allpresel<-apply(Variable_preselected[,-(which(colnames(Variable_preselected)=="seq"):ncol(Variable_preselected))],2,function(X){sum(X,na.rm=TRUE)/length(X)})
        
        allpresel2<-allpresel[order(allpresel,decreasing=TRUE)]
        allranks<-apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],2,function(X){na.exclude(X)})
        Var_ranks<-lapply(allranks,as.vector)
        #per variable, proportion of models in which they where preselected
        # barplot(allpresel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Preselection"),las = 2)
        #ggplot
        allpresel3<-data.frame(value=allpresel2,variable=factor(names(allpresel2),levels=names(allpresel2)),categ=ifelse(names(allpresel2)%in%climatic,"climatic",
                                                                                                                         ifelse(names(allpresel2)%in%edaphic,"edaphic",
                                                                                                                                ifelse(names(allpresel2)%in%topographic,"topographic",
                                                                                                                                       ifelse(names(allpresel2)%in%landcover,"landcover",
                                                                                                                                              ifelse(names(allpresel2)%in%remote.sensing,"remote.sensing",names(allpresel2)))))))
        allpresel3$categ<-factor(allpresel3$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        # allsel3$variable<-as.character(allsel3$variable)
        # if(allsel3$variable[allsel3$variable=="Altitude"|allsel3$variable=="Elevation"]!=colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]){
        #   allsel3$variable[which(allsel3$variable=="Altitude"|allsel3$variable=="Elevation")]<-colnames(ENVdata)[colnames(ENVdata)=="Elevation"|colnames(ENVdata)=="Altitude"]
        # }
        # allsel3$variable<-as.factor(allsel3$variable)
        #MERGE categ if same group of correlated vars
        categinplot_allpresel3<-unique(allpresel3$categ)
        allpresel3_10<-allpresel3[1:10,]
        categinplot_allpresel3_10<-unique(allpresel3[1:10,]$categ)   
        assign(paste0("pPres_",GtoM,"_",algo,"_",PAAB),ggplot(allpresel3, aes(x=variable,y=value,fill=categ))+ theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
            scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allpresel3])+
              scale_y_continuous(limits=c(0,1)))
              #pPres_FU_GLM_PA
        
        #per variable, proportion of models in which they where selected
        allsel<-apply(Variable_ranks[,-(which(colnames(Variable_ranks)=="seq"):ncol(Variable_ranks))],2,function(X){ifelse(is.na(X),NA,1)})
        allsel_group<-matrix(NA,ncol=length(groups),nrow=nrow(allsel))
        rownames(allsel_group)<-rownames(allsel)
        colnames(allsel_group)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=19
          if(length(groups[[i]])>1){
            allsel_group[,groups[[i]][1]]<-apply(allsel[,colnames(allsel)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE)!=0)})
          }else{
            allsel_group[,groups[[i]][1]]<-allsel[,colnames(allsel)%in%groups[[i]]]
          }
        }
        
        allsel2<-apply(allsel_group,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        # allsel2["Soil_mois"]
        # barplot(allsel2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Selection"),las = 2)
        allsel2<-allsel2[order(allsel2,decreasing=TRUE)]
        #proportion of models in which variable selected

        
        allsel3<-data.frame(value=allsel2,variable=factor(names(allsel2),levels=names(allsel2)),categ=ifelse(names(allsel2)%in%climatic,"climatic",
                                                                                                             ifelse(names(allsel2)%in%edaphic,"edaphic",
                                                                                                                    ifelse(names(allsel2)%in%topographic,"topographic",
                                                                                                                           ifelse(names(allsel2)%in%landcover,"landcover",
                                                                                                                                  ifelse(names(allsel2)%in%remote.sensing,"remote.sensing",names(allsel2)))))))
        allsel3$categ<-factor(allsel3$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        # categ<-vector()
        # for( i in 1:length(allsel3$variable)){
        #   vartest<-allsel3$variable[i]
        #   categ[i]<-paste0(groups[[c(which(unlist(lapply(groups,function(X){vartest%in%X}))))]])
        # }
        # 
        # allsel3$categ<-categ
        # allsel3$variable<-recode_factor(allsel3$variable,ndvi="pH",pH.1="pH",bio11_tco="winter_T",Soil_aera="aeration",bulkSoilW="soilwater",EC_1_5="conduct",bio3_tiso="isotherm",bio2_tdr="TempVar",sRadY="solar_rad",bio1_T="AirTemp",Phyllosil="Silicate",bio15_ps="var_prec",lowVeg_ag="lu_lowveg")
        allsel3$variable<-recode_factor(allsel3$variable,bio1_t="AirTemp",bio2_tdr="day_var_temp",sRadY="solar_rad",bulkSoilW="C_N",TotalP="Phosphorus",Quartz="Silicate",TiO2="Oxides",ThinSilt="granulo")
        allsel3$variable<-with(allsel3,reorder(variable,value,decreasing=TRUE))
        if(sum(as.character(allsel3$variable)=="pH")>1){
          allsel3$variable<-factor(allsel3$variable,levels=as.character(allsel3$variable[-length(allsel3$variable)]))
        } else{allsel3$variable<-factor(allsel3$variable,levels=as.character(allsel3$variable))}
        categinplot_allsel3<-unique(allsel3$categ)   
        allsel3_10<-allsel3[1:10,]
        print(allsel3_10)
        categinplot_allsel3_10<-unique(allsel3[1:10,]$categ)   
        assign(paste0("pSel_",GtoM,"_",algo,"_",PAAB),ggplot(allsel3, aes(x=variable,y=value,fill=categ))+ theme_bw() + geom_bar(stat='identity') + xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pSel10_",GtoM,"_",algo,"_",PAAB),ggplot(allsel3_10, aes(x=variable,y=value,fill=categ))+ theme_bw() + geom_bar(stat='identity') + xlab("")+ ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_allsel3_10])+
                 scale_y_continuous(limits=c(0,1)))
              #pSel_FU_GLM_PA
        
        
        # #list of ranks obtained per variable
        # allranks<-apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],2,function(X){na.exclude(X)})
        # # summary(Variable_ranks_BA)
        # Var_ranks<-lapply(allranks,as.vector)
        # #per variable, it's mean rank whithin postselected variables
        # #per variable, it's mean rank whithin postselected variables
        # reordering<-with(melt(Var_ranks),                       # Order boxes by median
        #      reorder(L1,
        #              value,
        #              median))
        # 
        # data_reordered<-melt(Var_ranks)
        # data_reordered$L1<-factor(data_reordered$L1,
        #        levels = levels(reordering))
        # data_reordered$categ<-unlist(lapply(as.vector(data_reordered$L1),function(X){ifelse(X%in%climatic,"climatic",
        #                                                                                     ifelse(X%in%edaphic,"edaphic",
        #                                                                                            ifelse(X%in%topographic,"topographic",
        #                                                                                                   ifelse(X%in%remote_sensing,"remote_sensing",
        #                                                                                                   ifelse(X%in%landcover,"landcover",X)))))}))
        # categinplot_pRank<-unique(data_reordered$categ)   
        # assign(paste0("pRank_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered,aes(x=L1,y=value,fill=categ))+ theme_bw()+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
        #          theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #                legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
        #          scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_pRank])+
        #          scale_y_continuous(limits=c(0,1)))
        # # boxplot(value~testing,melt(Var_ranks),main=paste0(PAAB,"_",GtoM,"_",algo,"_Rank_Distribution_when_selected"),las = 2)
        # # boxplot(value~L1,melt(Var_ranks),main=paste0(PAAB,"_",GtoM,"_",algo,"_Rank_Distribution_when_selected"),las = 2)
        # 
        #       #pRank_PR_GBM_PA
        #careful, possibility to have a low rank in a model with few variables
        
        
        
        
        #per variable, raw variable importance (unit depends on algo), 
        #all model performances weighted the same (true ecological importance depends on model performance)
        # boxplot(Variable_importance[,-which(colnames(Variable_importance)=="seq")],main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_all"),las = 2)
        
        # allsel3<-apply(Variable_importance[,-which(colnames(Variable_ranks)=="seq")]*allsel,2,function(X){na.exclude(X)})
        # boxplot(allsel3,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_only_when_selected"),las = 2)
        
        
        
        #variable importance corrected by the maximal value of importance across variables (best =1)
        Variable_importance1<-t(apply(Variable_importance[,-(which(colnames(Variable_importance)=="seq"):ncol(Variable_importance))],1,function(X){X/max(X)}))

        Variable_importance2<-matrix(NA,ncol=length(groups),nrow=nrow(Variable_importance1))
        rownames(Variable_importance2)<-rownames(Variable_importance1)
        colnames(Variable_importance2)<-unlist(lapply(groups,function(X){X[1]}))
        for (i in 1:length(groups)){#i=1
          if(length(groups[[i]])>1){
            Variable_importance2[,groups[[i]][1]]<-apply(Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE))})
          } else{
            Variable_importance2[,groups[[i]][1]]<-Variable_importance1[,colnames(Variable_importance1)%in%groups[[i]]]
          }
        }
        # summary(Variable_importance2);summary(Variable_importance1)
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        #keep value only if variable is selected in model
        
        
        # boxplot(Variable_importance2,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_standard_all"),las = 2)
        allsel4<-apply(Variable_importance2*allsel_group,2,function(X){na.exclude(X)})
        
        
        data_reordered4<-melt(allsel4)
        data_reordered4$categ<-unlist(lapply(as.vector(data_reordered4$L1),function(X){ifelse(X%in%climatic,"climatic",
                                                                                              ifelse(X%in%edaphic,"edaphic",
                                                                                                     ifelse(X%in%topographic,"topographic",
                                                                                                                   ifelse(X%in%landcover,"landcover",
                                                                                                                          ifelse(X%in%remote.sensing,"remote.sensing",X)))))}))
        data_reordered4$L1<-recode_factor(data_reordered4$L1,bio1_t="AirTemp",bio2_tdr="day_var_temp",sRadY="solar_rad",bulkSoilW="C_N",TotalP="Phosphorus",Quartz="Silicate",TiO2="Oxides",ThinSilt="granulo")
        reordering4<-with(data_reordered4,                       # Order boxes by median
                          reorder(L1,
                                  value,
                                  median,decreasing=TRUE))
        data_reordered4$L1<-factor(data_reordered4$L1,
                                   levels = levels(reordering4))
        
        data_reordered4$categ<-factor(data_reordered4$categ,levels=c("climatic","edaphic","landcover","remote.sensing","topographic"))
        categinplot_VIMP<-unique(data_reordered4$categ)
        data_reordered4_10<-data_reordered4[data_reordered4$L1%in%levels(data_reordered4$L1)[1:10],]
        categinplot_VIMP_10<-unique(data_reordered4_10$categ)
        

        # data_reordered4[data_reordered4$L1%in%levels(reordering4)[1:10],]
        assign(paste0("pVIMP_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered4,aes(x=L1,y=value,fill=categ))+ theme_bw()+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist"))) +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP])+
                 scale_y_continuous(limits=c(0,1)))
        assign(paste0("pVIMP10_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered4_10,aes(x=L1,y=value,fill=categ))+ theme_bw()+geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist"))) +
                 theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       legend.text = element_text(size=5), legend.key.size = unit(0.5,"line"),legend.position = c(0.9,0.8),legend.title = element_blank())+ 
                 scale_fill_manual(values=col_table$color[col_table$categ%in%categinplot_VIMP_10])+
                 scale_y_continuous(limits=c(0,1)))
        #pVIMP10_FU_GLM_PA
        
        #% of models selecting a variable from each group
        allsel5_temp<-apply(allsel_group,1,function(X){c(any(X[names(X)%in%climatic]>=1),any(X[names(X)%in%edaphic]>=1),any(X[names(X)%in%landcover]>=1),any(X[names(X)%in%remote.sensing]>=1), any(X[names(X)%in%topographic]>=1))})
        allsel5_temp<-t(allsel5_temp)
        colnames(allsel5_temp)<-c("climatic","edaphic","landcover","remote.sensing","topographic")
        allsel5_temp<-apply(allsel5_temp,2,function(X){sum(X,na.rm=TRUE)/length(X)})
        Prop_models_one_cov[paste0(PAAB,"_",algo,"_",GtoM),]<-allsel5_temp
        # allsel5_temp<-allsel5_temp[order(allsel5_temp,decreasing=TRUE)]
        allsel5<-data.frame(value=allsel5_temp,grouping=factor(names(allsel5_temp),levels=names(allsel5_temp)))
        print(allsel5)
        # assign(paste0("pgrouped_",GtoM,"_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping)) + theme_bw() + 
        #          geom_bar(stat='identity', show.legend = FALSE) + 
        #          xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist"))) + ylim(c(0,1))+
        #          theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=5),plot.margin=unit(c(0,0,0,0),"pt"),
        #                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        #          scale_fill_manual(values=col_table$color[order(match(col_table$categ,levels(allsel5$grouping)))]))
        # #     pgrouped_FU_GLM_PA    

          assign(paste0("pgrouped_",GtoM,"_",algo,"_",PAAB),ggplot(allsel5, aes(x=grouping,y=value,fill=grouping)) + theme_bw() + 
                   geom_bar(stat='identity', show.legend = FALSE) + 
                   ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist"))) + ylim(c(0,1))+
                   theme(axis.text.x=element_blank(),plot.margin=unit(c(2,4,2,4),"pt"),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         axis.title.x=element_blank()) +
                   scale_fill_manual(values=rep("grey",length(allsel5$grouping))))
        # pgrouped_FU_GLM_PA
        # pgrouped_BA_RF_AB
        #variable importance corrected by model performance (attribute the full performance to all variables according to their importance)
        # if(PAAB=="PA"){
        #   load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Eval.Rda"))
        #   TSSadj<-ifelse(Eval$TSS_adj>0,Eval$TSS_adj,NA)
        #   varimp_percentage<-t(apply(Variable_importance[,-which(colnames(Variable_ranks)=="seq")],1,function(X){X/sum(X,na.rm=TRUE)}))
        #   varimp_percentage<-varimp_percentage*allsel
        #   varimp_percentage_weighted<-apply(varimp_percentage,2,function(X){X*TSSadj})
        #   
        #   # boxplot(varimp_percentage_weighted,main=paste0(PAAB,"_",GtoM,"_",algo,"_Varimp_distribution_weighted_by_ModPerf"),las = 2)
        #   reordering5<-with(melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),],                       # Order boxes by median
        #                     reorder(Var2,
        #                             value,
        #                             median,decreasing=TRUE))
        #   
        #   data_reordered5<-melt(varimp_percentage_weighted)[!is.na(melt(varimp_percentage_weighted)$value),]
        #   data_reordered5$Var2<-factor(data_reordered5$Var2,
        #                              levels = levels(reordering5))
        #   assign(paste0("pVIMPw_",GtoM,"_",algo,"_",PAAB),ggplot(data_reordered5,aes(x=Var2,y=value))+ theme_bw() +geom_boxplot()+ xlab("") + ylab(ifelse(PAAB=="AB","",ifelse(GtoM=="FU","Fungi","Protist")))  +
        #            theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=3),plot.margin=unit(c(0,0,0,0),"pt"),
        #                  panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
        #   #pVIMPw_FU_GLM_PA
        # }
      }#if PRFU
    }
  }
}




# tPA <- textGrob("Presence-Absence models",hjust=-.2)
# tAB <- textGrob("Rel. Abundance models",hjust=-.1)

# plot(Ppresel_GLM)
Ppresel_GLM<-grid.arrange(grobs=list(pPres_AR_GLM_PA,pPres_AR_GLM_AB,pPres_BA_GLM_PA,pPres_BA_GLM_AB,pPres_FU_GLM_PA,pPres_FU_GLM_AB,pPres_PR_GLM_PA,pPres_PR_GLM_AB),
                          nrow=4)
Psel_GLM<-grid.arrange(pSel_AR_GLM_PA,pSel_AR_GLM_AB,pSel_BA_GLM_PA,pSel_BA_GLM_AB,pSel_FU_GLM_PA,pSel_FU_GLM_AB,pSel_PR_GLM_PA,pSel_PR_GLM_AB,nrow=4)
Psel10_GLM<-grid.arrange(pSel10_AR_GLM_PA,pSel10_AR_GLM_AB,pSel10_BA_GLM_PA,pSel10_BA_GLM_AB,pSel10_FU_GLM_PA,pSel10_FU_GLM_AB,pSel10_PR_GLM_PA,pSel10_PR_GLM_AB,nrow=4)
# PRank_GLM<-grid.arrange(pRank_BA_GLM_PA,pRank_BA_GLM_AB,pRank_AR_GLM_PA,pRank_AR_GLM_AB,pRank_FU_GLM_PA,pRank_FU_GLM_AB,pRank_PR_GLM_PA,pRank_PR_GLM_AB,nrow=4)
PVIMP_GLM<-grid.arrange(pVIMP_AR_GLM_PA,pVIMP_AR_GLM_AB,pVIMP_BA_GLM_PA,pVIMP_BA_GLM_AB,pVIMP_FU_GLM_PA,pVIMP_FU_GLM_AB,pVIMP_PR_GLM_PA,pVIMP_PR_GLM_AB,nrow=4)
PVIMP10_GLM<-grid.arrange(pVIMP10_AR_GLM_PA,pVIMP10_AR_GLM_AB,pVIMP10_BA_GLM_PA,pVIMP10_BA_GLM_AB,pVIMP10_FU_GLM_PA,pVIMP10_FU_GLM_AB,pVIMP10_PR_GLM_PA,pVIMP10_PR_GLM_AB,
                          nrow=4,heights=c(2,2,2,2))
Pgrouped_GLM<-grid.arrange(pgrouped_AR_GLM_PA,pgrouped_AR_GLM_AB,pgrouped_BA_GLM_PA,pgrouped_BA_GLM_AB,pgrouped_FU_GLM_PA,pgrouped_FU_GLM_AB,pgrouped_PR_GLM_PA,pgrouped_PR_GLM_AB,
                           nrow=4,heights=c(2.5,2,2,2))

# PVIMPw_GLM<-grid.arrange(pVIMPw_BA_GLM_PA,pVIMPw_AR_GLM_PA,pVIMPw_FU_GLM_PA,pVIMPw_PR_GLM_PA,nrow=4)

Ppresel_GAM<-grid.arrange(pPres_AR_GAM_PA,pPres_AR_GAM_AB,pPres_BA_GAM_PA,pPres_BA_GAM_AB,pPres_FU_GAM_PA,pPres_FU_GAM_AB,pPres_PR_GAM_PA,pPres_PR_GAM_AB,nrow=4)
Psel_GAM<-grid.arrange(pSel_AR_GAM_PA,pSel_AR_GAM_AB,pSel_BA_GAM_PA,pSel_BA_GAM_AB,pSel_FU_GAM_PA,pSel_FU_GAM_AB,pSel_PR_GAM_PA,pSel_PR_GAM_AB,nrow=4)
Psel10_GAM<-grid.arrange(pSel10_AR_GAM_PA,pSel10_AR_GAM_AB,pSel10_BA_GAM_PA,pSel10_BA_GAM_AB,pSel10_FU_GAM_PA,pSel10_FU_GAM_AB,pSel10_PR_GAM_PA,pSel10_PR_GAM_AB,nrow=4)
# PRank_GAM<-grid.arrange(pRank_BA_GAM_PA,pRank_BA_GAM_AB,pRank_AR_GAM_PA,pRank_AR_GAM_AB,pRank_FU_GAM_PA,pRank_FU_GAM_AB,pRank_PR_GAM_PA,pRank_PR_GAM_AB,nrow=4)
PVIMP_GAM<-grid.arrange(pVIMP_AR_GAM_PA,pVIMP_AR_GAM_AB,pVIMP_BA_GAM_PA,pVIMP_BA_GAM_AB,pVIMP_FU_GAM_PA,pVIMP_FU_GAM_AB,pVIMP_PR_GAM_PA,pVIMP_PR_GAM_AB,nrow=4)
PVIMP10_GAM<-grid.arrange(pVIMP10_AR_GAM_PA,pVIMP10_AR_GAM_AB,pVIMP10_BA_GAM_PA,pVIMP10_BA_GAM_AB,pVIMP10_FU_GAM_PA,pVIMP10_FU_GAM_AB,pVIMP10_PR_GAM_PA,pVIMP10_PR_GAM_AB,nrow=4)
Pgrouped_GAM<-grid.arrange(pgrouped_AR_GAM_PA,pgrouped_AR_GAM_AB,pgrouped_BA_GAM_PA,pgrouped_BA_GAM_AB,pgrouped_FU_GAM_PA,pgrouped_FU_GAM_AB,pgrouped_PR_GAM_PA,pgrouped_PR_GAM_AB,nrow=4)

# PVIMPw_GAM<-grid.arrange(pVIMPw_BA_GAM_PA,pVIMPw_AR_GAM_PA,pVIMPw_FU_GAM_PA,pVIMPw_PR_GAM_PA,nrow=4)

Ppresel_GBM<-grid.arrange(pPres_AR_GBM_PA,pPres_AR_GBM_AB,pPres_BA_GBM_PA,pPres_BA_GBM_AB,pPres_FU_GBM_PA,pPres_FU_GBM_AB,pPres_PR_GBM_PA,pPres_PR_GBM_AB,nrow=4)
Psel_GBM<-grid.arrange(pSel_AR_GBM_PA,pSel_AR_GBM_AB,pSel_BA_GBM_PA,pSel_BA_GBM_AB,pSel_FU_GBM_PA,pSel_FU_GBM_AB,pSel_PR_GBM_PA,pSel_PR_GBM_AB,nrow=4)
Psel10_GBM<-grid.arrange(pSel10_AR_GBM_PA,pSel10_AR_GBM_AB,pSel10_BA_GBM_PA,pSel10_BA_GBM_AB,pSel10_FU_GBM_PA,pSel10_FU_GBM_AB,pSel10_PR_GBM_PA,pSel10_PR_GBM_AB,nrow=4)
# PRank_GBM<-grid.arrange(pRank_BA_GBM_PA,pRank_BA_GBM_AB,pRank_AR_GBM_PA,pRank_AR_GBM_AB,pRank_FU_GBM_PA,pRank_FU_GBM_AB,pRank_PR_GBM_PA,pRank_PR_GBM_AB,nrow=4)
PVIMP_GBM<-grid.arrange(pVIMP_AR_GBM_PA,pVIMP_AR_GBM_AB,pVIMP_BA_GBM_PA,pVIMP_BA_GBM_AB,pVIMP_FU_GBM_PA,pVIMP_FU_GBM_AB,pVIMP_PR_GBM_PA,pVIMP_PR_GBM_AB,nrow=4)
PVIMP10_GBM<-grid.arrange(pVIMP10_AR_GBM_PA,pVIMP10_AR_GBM_AB,pVIMP10_BA_GBM_PA,pVIMP10_BA_GBM_AB,pVIMP10_FU_GBM_PA,pVIMP10_FU_GBM_AB,pVIMP10_PR_GBM_PA,pVIMP10_PR_GBM_AB,nrow=4)
Pgrouped_GBM<-grid.arrange(pgrouped_AR_GBM_PA,pgrouped_AR_GBM_AB,pgrouped_BA_GBM_PA,pgrouped_BA_GBM_AB,pgrouped_FU_GBM_PA,pgrouped_FU_GBM_AB,pgrouped_PR_GBM_PA,pgrouped_PR_GBM_AB,nrow=4)
# PVIMPw_GBM<-grid.arrange(pVIMPw_BA_GBM_PA,pVIMPw_AR_GBM_PA,pVIMPw_FU_GBM_PA,pVIMPw_PR_GBM_PA,nrow=4)

Ppresel_RF<-grid.arrange(pPres_AR_RF_PA,pPres_AR_RF_AB,pPres_BA_RF_PA,pPres_BA_RF_AB,pPres_FU_RF_PA,pPres_FU_RF_AB,pPres_PR_RF_PA,pPres_PR_RF_AB,nrow=4)
Psel_RF<-grid.arrange(pSel_AR_RF_PA,pSel_AR_RF_AB,pSel_BA_RF_PA,pSel_BA_RF_AB,pSel_FU_RF_PA,pSel_FU_RF_AB,pSel_PR_RF_PA,pSel_PR_RF_AB,nrow=4)
Psel10_RF<-grid.arrange(pSel10_AR_RF_PA,pSel10_AR_RF_AB,pSel10_BA_RF_PA,pSel10_BA_RF_AB,pSel10_FU_RF_PA,pSel10_FU_RF_AB,pSel10_PR_RF_PA,pSel10_PR_RF_AB,nrow=4)
# PRank_RF<-grid.arrange(pRank_BA_RF_PA,pRank_BA_RF_AB,pRank_AR_RF_PA,pRank_AR_RF_AB,pRank_FU_RF_PA,pRank_FU_RF_AB,pRank_PR_RF_PA,pRank_PR_RF_AB,nrow=4)
PVIMP_RF<-grid.arrange(pVIMP_AR_RF_PA,pVIMP_AR_RF_AB,pVIMP_BA_RF_PA,pVIMP_BA_RF_AB,pVIMP_FU_RF_PA,pVIMP_FU_RF_AB,pVIMP_PR_RF_PA,pVIMP_PR_RF_AB,nrow=4)
PVIMP10_RF<-grid.arrange(pVIMP10_AR_RF_PA,pVIMP10_AR_RF_AB,pVIMP10_BA_RF_PA,pVIMP10_BA_RF_AB,pVIMP10_FU_RF_PA,pVIMP10_FU_RF_AB,pVIMP10_PR_RF_PA,pVIMP10_PR_RF_AB,nrow=4)
Pgrouped_RF<-grid.arrange(pgrouped_AR_RF_PA,pgrouped_AR_RF_AB,pgrouped_BA_RF_PA,pgrouped_BA_RF_AB,pgrouped_FU_RF_PA,pgrouped_FU_RF_AB,pgrouped_PR_RF_PA,pgrouped_PR_RF_AB,nrow=4)
# PVIMPw_RF<-grid.arrange(pVIMPw_BA_RF_PA,pVIMPw_AR_RF_PA,pVIMPw_FU_RF_PA,pVIMPw_PR_RF_PA,nrow=4)




pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Presel_GLMGAMGBMRF.pdf"))
plot(Ppresel_GLM)
plot(Ppresel_GAM)
plot(Ppresel_GBM)
plot(Ppresel_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Presel.png"),res=300,width=1961,height=1500)
plot(Ppresel_RF)
dev.off()

pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Sel_GLMGAMGBMRF.pdf"))
plot(Psel_GLM)
plot(Psel_GAM)
plot(Psel_GBM)
plot(Psel_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Sel.png"),res=300,width=1961,height=1500)
plot(Psel_RF)
dev.off()

pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Sel10_GLMGAMGBMRF.pdf"))
plot(Psel10_GLM)
plot(Psel10_GAM)
plot(Psel10_GBM)
plot(Psel10_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_Sel10.png"),res=300,width=1961,height=1500)
plot(Psel10_RF)
dev.off()


# 
# pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_Rank_GLMGAMGBMRF.pdf"))
# plot(PRank_GLM)
# plot(PRank_GAM)
# plot(PRank_GBM)
# plot(PRank_RF)
# dev.off()
# 
# png(file=paste0("figures/PAAB_selection/VARimp/GLM_Rank.png"),res=300,width=1961,height=1500)
# plot(PRank_GLM)
# dev.off()
# png(file=paste0("figures/PAAB_selection/VARimp/GAM_Rank.png"),res=300,width=1961,height=1500)
# plot(PRank_GAM)
# dev.off()
# png(file=paste0("figures/PAAB_selection/VARimp/GBM_Rank.png"),res=300,width=1961,height=1500)
# plot(PRank_GBM)
# dev.off()
# png(file=paste0("figures/PAAB_selection/VARimp/RF_Rank.png"),res=300,width=1961,height=1500)
# plot(PRank_RF)
# dev.off()


pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_VIMP_GLMGAMGBMRF.pdf"))
plot(PVIMP_GLM)
plot(PVIMP_GAM)
plot(PVIMP_GBM)
plot(PVIMP_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_VIMP.png"),res=300,width=1961,height=1500)
plot(PVIMP_RF)
dev.off()

pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_VIMP10_GLMGAMGBMRF.pdf"))
plot(PVIMP10_GLM)
plot(PVIMP10_GAM)
plot(PVIMP10_GBM)
plot(PVIMP10_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_VIMP10.png"),res=300,width=1961,height=1500)
plot(PVIMP10_RF)
dev.off()

png(file=paste0("figures/PAAB_selection/VARimp/GLM_groups.png"),res=300,width=1961,height=1500)
plot(Pgrouped_GLM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GAM_groups.png"),res=300,width=1961,height=1500)
plot(Pgrouped_GAM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/GBM_groups.png"),res=300,width=1961,height=1500)
plot(Pgrouped_GBM)
dev.off()
png(file=paste0("figures/PAAB_selection/VARimp/RF_groups.png"),res=300,width=1961,height=1500)
plot(Pgrouped_RF)
dev.off()




# pdf(file=paste0("figures/PAAB_selection/VARimp/All_varimp_VIMPw_GLMGAMGBMRF.pdf"))
# plot(PVIMPw_GLM)
# plot(PVIMPw_GAM)
# plot(PVIMPw_GBM)
# plot(PVIMPw_RF)
# dev.off()
# 
# png(file=paste0("figures/PAAB_selection/VARimp/GLM_VIMPw.png"),res=300,width=1961,height=1500)
# plot(PVIMPw_GLM)
# dev.off()
# png(file=paste0("figures/PAAB_selection/VARimp/GAM_VIMPw.png"),res=300,width=1961,height=1500)
# plot(PVIMPw_GAM)
# dev.off()
# png(file=paste0("figures/PAAB_selection/VARimp/GBM_VIMPw.png"),res=300,width=1961,height=1500)
# plot(PVIMPw_GBM)
# dev.off()
# png(file=paste0("figures/PAAB_selection/VARimp/RF_VIMPw.png"),res=300,width=1961,height=1500)
# plot(PVIMPw_RF)
# dev.off()


# ######
# PAAB="PA"
# # PAAB="AB"
# GtoM="BA"
# algo="GLM"
# # load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
# #colnames(ENVdata)
# 
# 
# climatic<-c("bio1_t", "bio10_twa", "bio11_tco", "bio12_p" ,"bio13_pwe","bio14_pdr","bio15_ps" ,"bio16_pwe","bio17_pdr","bio18_pwa","bio19_pco",
#             "bio2_tdr" ,"bio3_tiso","bio4_ts" , "bio5_tmax","bio6_tmin","bio7_tar", "bio8_twet","bio9_tdry","GDD0" , "ETP" ,  "cumday_no","sRadY"  )
# edaphic<-c("pH","pH.1","bulkSoilW","soilTemp","EC_1_5",    "TotalP",    "Nitrogen",  "Carbon",    "Hydrogen",  "Phyllosil", "Quartz",    "Feldspath","Plagiocla","MassiveLi",
#            "Calcite",   "Indoses",   "SiO2",      "TiO2",     "Al2O3",     "Fe2O3" ,"MnO"  , "MgO",   "CaO",   "Na2O","Marlyshal","MarlShale",
#            "K2O" ,  "P2O5"  , "OM"  ,  "Cr2O3",  "NiO",   "d15N" ,  "d13C","Silt_clay", "clay"  , "ThinSilt", "ThickSilt", "ThinSand", "ThickSand","Soil_aera","Soil_humu","Soil_mois","Soil_mois.1","Soil_nutr" )
# topographic<-c("aspect","slope","Elevation","Altitude")
# landcover<-c("forest_ag", "hydro_agg", "lowVeg_ag", "anthropos",   "deciduous")
# remote_sensing<-c("ndmi","ndvi")
# 
# load(paste0(PAAB,"/",GtoM,"/data/ENVdata.Rda"))
# colnames(ENVdata)<-substr(colnames(ENVdata),1,9)
# ENVdataforcor<-ENVdata
# ENVdataforcor<-apply(ENVdataforcor, 2, function(x) as.numeric(as.character(x)))
# cluster<-hclust(as.dist(1-abs(cor(apply(ENVdata,2,as.numeric), use = "pairwise.complete.obs"))))
# cuttree<-cutree(cluster,h=0.2)
# groups<-list()
# for (i in 1:cuttree[length(cuttree)]){
#   groups[[i]]<-names(cuttree)[cuttree==i]
# }
# if(length(unlist(lapply(groups,function(X){any(X=="Soil_mois")})))>1){
#   groups[[which(unlist(lapply(groups,function(X){any(X=="Soil_mois")})))[2]]]<-NULL
# }
# if(sum(unlist(lapply(groups,function(X){any(X=="pH"|X=="pH.1")})))>1){
#   groups[[which(unlist(lapply(groups,function(X){any(X=="pH.1")})))]]<-NULL
#   groups[[which(unlist(lapply(groups,function(X){any(X=="pH")})))]]<-c("pH","pH.1")
# }
# load(paste0(PAAB,"/",GtoM,"/data/",algo,"/Variable_ranks.Rda"))
# load(paste0("../../ASV_data/ASV_taxo/",GtoM,"taxo_grouped_phylum.Rda"))
# Taxo<-eval(parse(text=paste0(GtoM,"taxo_grouped_phylum")))
# Variable_ranks<- cbind(Variable_ranks[Variable_ranks$seq%in%Taxo$Seq,],Taxo[Taxo$Seq%in%Variable_ranks$seq,])
# all(Variable_ranks$Seq==Variable_ranks$seq)
# names(Variable_ranks)<-substr(names(Variable_ranks),1,9)
# allsel<-apply(Variable_ranks[,-which(colnames(Variable_ranks)=="seq")],2,function(X){ifelse(is.na(X),NA,1)})
# allsel_group<-matrix(NA,ncol=length(groups),nrow=nrow(allsel))
# rownames(allsel_group)<-rownames(allsel)
# colnames(allsel_group)<-unlist(lapply(groups,function(X){X[1]}))
# for (i in 1:length(groups)){#i=1
#   if(length(groups[[i]])>1){
#     allsel_group[,groups[[i]][1]]<-apply(allsel[,colnames(allsel)%in%groups[[i]]],1,function(X){ifelse(all(is.na(X)),NA,sum(X,na.rm=TRUE)!=0)})
#   }else{
#     allsel_group[,groups[[i]][1]]<-allsel[,colnames(allsel)%in%groups[[i]]]
#   }
# }

summary(Prop_models_one_cov)
boxplot(Prop_models_one_cov)
