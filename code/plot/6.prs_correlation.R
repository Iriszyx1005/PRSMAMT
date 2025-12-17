require(data.table)
library(dplyr)
library(caret)
library(pROC)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggcorrplot)
library(corrplot)
###32prs_european_corr_heatmap####
setwd("/home/sshen/Disk_m2/PRS/plot")
total=NULL

#######onco_white
phe=fread("/home/sshen/Disk_m2/PRS_yxzhang/train/white/pheno_29095.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/Onco_white_prs.txt")
LC=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/LC_Onco_EUR_PRS_29095.txt"))
LC=LC[,-"FID"]
prs=merge(prs,LC,by="IID")
onco=prs[,2:34]
onco[,"PGS004691"] = -onco[,"PGS004691"]
onco=as.matrix(onco)
onco = onco[,order(colnames(onco))]

onco=onco[,1:32]
d = apply(onco,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'onco' = 'value')


total=rbind(total,data_cor)

#######PLCO
phe=fread("/home/sshen/Disk_m2/PRS_yxzhang/validation/PLCO/pheno_98651.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/PLCO_white_prs.txt")
LC=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/LC_PLCO_PRS_98651.txt"))
LC=LC[,-"FID"]
prs=merge(prs,LC,by="IID")
PLCO=prs[,2:34]
PLCO[,"PGS004691"] = -PLCO[,"PGS004691"]
PLCO=as.matrix(PLCO)
PLCO = PLCO[,order(colnames(PLCO))]
PLCO=PLCO[,1:32]
d = apply(PLCO,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'PLCO' = 'value')

total=merge(total,data_cor,by=c("x","y"))

#######TRICL_white
load("/home/sshen/Disk_m2/PRS/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==1,]
LC=fread("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/TRICL_white.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/Onco_TRICL_white_prs.txt")

prs=merge(prs,LC,by="IID")
TRICL=prs[,2:34]
TRICL[,"PGS004691"] = -TRICL[,"PGS004691"]
TRICL=as.matrix(TRICL)
TRICL = TRICL[,order(colnames(TRICL))]
TRICL=TRICL[,1:32]
d = apply(TRICL,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")

data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'TRICL' = 'value')

total=merge(total,data_cor,by=c("x","y"))


#######AOU_white
load("/home/sshen/Disk_m2/PRS/AOU/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eur",]
phe$Lungcancer = phe$lung_ca
load("/home/sshen/Disk_m2/PRS/AOU/PRS_pub_score_AOU.RData")
LC=prs
LC=LC[,-"#FID"]
load("/home/sshen/Disk_m2/PRS/AOU/LungCa_V3_eur.RData")
prs=prs[,c("IID","score")]
setnames(prs,"score","PRS")

prs=merge(prs,LC,by="IID")
AOU=prs[,2:34]
AOU[,"PGS004691"] = -AOU[,"PGS004691"]
AOU=as.matrix(AOU)
AOU = AOU[,order(colnames(AOU))]
AOU=AOU[,1:32]
d = apply(AOU,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'AOU' = 'value')

total=merge(total,data_cor,by=c("x","y"))

total$Corr=(total$onco+total$PLCO+total$TRICL+total$AOU)/4
total=total[,c("x","y","Corr")]

list <- rownames(cor)
list <- factor(list,levels = list)


data_cor_tri <- subset(total, match(x, list) >= match(y, list))


tiff(paste0("/home/sshen/Disk_m2/PRS/plot/32prs_european_corr_heatmap.tiff"),width = 12,height = 10,units = "in",res = 300)
p5=ggplot(data_cor_tri,aes(factor(x,levels = list),
                           factor(y,levels = list), 
                           fill=Corr))+  
  geom_tile(colour = "black")+  
  # geom_text(aes(label=Corr),size=2,colour="white")+coord_equal()+
  scale_fill_gradientn(colours = c("white",
                                   brewer.pal(7,"Set1")[1]),na.value = NA,
                       limits=c(0,1),breaks=c(0,0.5,1))+
  labs(x=NULL,y=NULL,title = "32 PRSs correlation in European")+
  scale_y_discrete(position = "right") +
  theme(panel.background=element_blank(),
        axis.text.x = element_text(angle=90,hjust=1,size=12,color="black",face="italic"),
        axis.text.y = element_text(hjust=1,size=12,color="black",face="italic"),
        legend.text = element_text(size = 15),     
        legend.title = element_text(size = 20),   
        axis.ticks.y.right = element_line(), 
        axis.title = element_text(size=20,color="black",face="italic"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        plot.margin = margin(5, 5, 5, 32))
p5
dev.off()

###32prs_asian_corr_heatmap####
setwd("/home/sshen/Disk_m2/PRS/plot")
total=NULL
#######onco_asian
phe=fread("/home/sshen/Disk_m2/PRS_yxzhang/train/asian/pheno_1974.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/ALL_onco_asian_prs.txt")
LC=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/LC_Onco_Asian_PRS_1974.txt"))
LC=LC[,-"FID"]
prs=merge(prs,LC,by="IID")
onco=prs[,2:34]
onco[,"PGS004691"] = -onco[,"PGS004691"]
onco[,"PGS003392"] = -onco[,"PGS003392"]
onco=as.matrix(onco)
onco = onco[,order(colnames(onco))]

onco=onco[,1:32]
d = apply(onco,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'onco' = 'value')


total=rbind(total,data_cor)

#######FLCCA
phe=fread("/home/sshen/Disk_m2/PRS_yxzhang/validation/FLCCA/pheno_8807.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/FLCCA_asian_prs.txt")
LC=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/LC_FLCCA_PRS_8807.txt"))
LC=LC[,-"FID"]
prs=merge(prs,LC,by="IID")
FLCCA=prs[,2:34]
FLCCA[,"PGS004691"] = -FLCCA[,"PGS004691"]
FLCCA[,"PGS003392"] = -FLCCA[,"PGS003392"]
FLCCA=as.matrix(FLCCA)
FLCCA = FLCCA[,order(colnames(FLCCA))]

FLCCA=FLCCA[,1:32]
d = apply(FLCCA,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'FLCCA' = 'value')

total=merge(total,data_cor,by=c("x","y"))

#######TRICL_asian
load("/home/sshen/Disk_m2/PRS/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==6,]
LC=fread("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/TRICL_asian.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/ALL_TRICL_asian_prs.txt")

prs=merge(prs,LC,by="IID")
TRICL=prs[,2:34]
TRICL[,"PGS004691"] = -TRICL[,"PGS004691"]
TRICL=as.matrix(TRICL)
TRICL = TRICL[,order(colnames(TRICL))]
TRICL=TRICL[,1:32]
d = apply(TRICL,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")

data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'TRICL' = 'value')

total=merge(total,data_cor,by=c("x","y"))


#######AOU_asian lung_ca
load("/home/sshen/Disk_m2/PRS/AOU/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eas",]
phe$Lungcancer = phe$lung_ca
load("/home/sshen/Disk_m2/PRS/AOU/PRS_pub_score_AOU.RData")
LC=prs
LC=LC[,-"#FID"]
load("/home/sshen/Disk_m2/PRS/AOU/LungCa_V3_eas.RData")
prs=prs[,c("IID","score")]
setnames(prs,"score","PRS")

prs=merge(prs,LC,by="IID")
AOU=prs[,2:34]
AOU[,"PGS004691"] = -AOU[,"PGS004691"]
AOU[,"PGS000156"] = -AOU[,"PGS000156"]
AOU[,"PGS000395"] = -AOU[,"PGS000395"]
AOU[,"PGS000396"] = -AOU[,"PGS000396"]
AOU[,"PGS000397"] = -AOU[,"PGS000397"]
AOU[,"PGS004442"] = -AOU[,"PGS004442"]
AOU[,"PGS004512"] = -AOU[,"PGS004512"]
AOU=as.matrix(AOU)
AOU = AOU[,order(colnames(AOU))]

AOU=AOU[,1:32]
d = apply(AOU,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'AOU' = 'value')

total=merge(total,data_cor,by=c("x","y"))

total$Corr=(total$onco+total$FLCCA+total$TRICL+total$AOU)/4
total=total[,c("x","y","Corr")]

#write.table(total,paste0("/home/sshen/Disk_m2/PRS/plot/32prs_asian.txt"),row.names = F,quote=F,sep="\t")

#p.mat <- cor_pmat(d)

list <- rownames(cor)
list <- factor(list,levels = list)


data_cor_tri <- subset(total, match(x, list) >= match(y, list))


tiff(paste0("/home/sshen/Disk_m2/PRS/plot/32prs_asian_corr_heatmap.tiff"),width = 12,height = 10,units = "in",res = 300)
p5=ggplot(data_cor_tri,aes(factor(x,levels = list),
                           factor(y,levels = list), 
                           fill=Corr))+  
  geom_tile(colour = "black")+  
  # geom_text(aes(label=Corr),size=2,colour="white")+coord_equal()+
  scale_fill_gradientn(colours = c("white",
                                   brewer.pal(7,"Set1")[1]),na.value = NA,
                       limits=c(0,1),breaks=c(0,0.5,1))+
  labs(x=NULL,y=NULL,title = "32 PRSs correlation in Asian")+
  scale_y_discrete(position = "right") +
  theme(panel.background=element_blank(),
        axis.text.x = element_text(angle=90,hjust=1,size=12,color="black",face="italic"),
        axis.text.y = element_text(hjust=1,size=12,color="black",face="italic"),
        legend.text = element_text(size = 15),     
        legend.title = element_text(size = 20),   
        axis.ticks.y.right = element_line(), 
        axis.title = element_text(size=20,color="black",face="italic"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        plot.margin = margin(5, 5, 5, 32))
p5
dev.off()

###32prs_African_corr_heatmap####
setwd("/home/sshen/Disk_m2/PRS/plot")
total=NULL

#######onco_african
phe=fread("/home/sshen/Disk_m2/PRS_yxzhang/train/black/pheno_532.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/ALL_onco_black_prs.txt")
LC=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/LC_Onco_AFR_PRS_532.txt"))
LC=LC[,-"FID"]
prs=merge(prs,LC,by="IID")
prs1=prs[,2:34]
prs1[,"PGS004691"] = -prs1[,"PGS004691"]
prs1[,"PGS004512"] = -prs1[,"PGS004512"]
prs1[,"PGS004442"] = -prs1[,"PGS004442"]
prs1[,"PGS003393"] = -prs1[,"PGS003393"]
prs1[,"PGS003391"] = -prs1[,"PGS003391"]
prs1[,"PGS000397"] = -prs1[,"PGS000397"]
onco=as.matrix(prs1)
onco = onco[,order(colnames(onco))]
onco=onco[,1:32]
d = apply(onco,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'onco' = 'value')


total=rbind(total,data_cor)
#######NCI
phe=fread("/home/sshen/Disk_m2/PRS_yxzhang/validation/NCI/pheno_5578.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/NCI_black_prs.txt")
LC=fread(paste0("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/LC_NCI_PRS_5578.txt"))
LC=LC[,-"FID"]
prs=merge(prs,LC,by="IID")
prs1=prs[,2:34]
prs1[,"PGS004691"] = -prs1[,"PGS004691"]
NCI=as.matrix(prs1)
NCI = NCI[,order(colnames(NCI))]
NCI=NCI[,1:32]
d = apply(NCI,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'NCI' = 'value')
total=merge(total,data_cor,by=c("x","y"))

#######TRICL_african
load("/home/sshen/Disk_m2/PRS/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==5,]
LC=fread("/home/sshen/Disk_m2/PRS_yxzhang/train/LC/TRICL_black.txt")
prs=fread("/home/sshen/Disk_m2/PRS/LungCa_V3/ALL_TRICL_black_prs.txt")

prs=merge(prs,LC,by="IID")
prs1=prs[,2:34]
prs1[,"PGS004691"] = -prs1[,"PGS004691"]
prs1[,"PGS000078"] = -prs1[,"PGS000078"]
prs1[,"PGS000156"] = -prs1[,"PGS000156"]
prs1[,"PGS000388"] = -prs1[,"PGS000388"]
prs1[,"PGS000393"] = -prs1[,"PGS000393"]
prs1[,"PGS000394"] = -prs1[,"PGS000394"]
prs1[,"PGS000396"] = -prs1[,"PGS000396"]
prs1[,"PGS000740"] = -prs1[,"PGS000740"]
prs1[,"PGS000789"] = -prs1[,"PGS000789"]
prs1[,"PGS002808"] = -prs1[,"PGS002808"]
prs1[,"PGS003391"] = -prs1[,"PGS003391"]
prs1[,"PGS003392"] = -prs1[,"PGS003392"]
prs1[,"PGS004442"] = -prs1[,"PGS004442"]
prs1[,"PGS004512"] = -prs1[,"PGS004512"]
prs1[,"PGS004884"] = -prs1[,"PGS004884"]
prs1[,"PGS004955"] = -prs1[,"PGS004955"]
TRICL=as.matrix(prs1)
TRICL = TRICL[,order(colnames(TRICL))]
TRICL=TRICL[,1:32]
d = apply(TRICL,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")

data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'TRICL' = 'value')

total=merge(total,data_cor,by=c("x","y"))


#######AOU_african
load("/home/sshen/Disk_m2/PRS/AOU/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="afr",]
phe$Lungcancer = phe$lung_ca
load("/home/sshen/Disk_m2/PRS/AOU/PRS_pub_score_AOU.RData")
LC=prs
LC=LC[,-"#FID"]
load("/home/sshen/Disk_m2/PRS/AOU/LungCa_V3_afr.RData")
prs=prs[,c("IID","score")]
setnames(prs,"score","PRS")

prs=merge(prs,LC,by="IID")
prs1=prs[,2:34]
prs1[,"PGS004691"] = -prs1[,"PGS004691"]
prs1[,"PGS000396"] = -prs1[,"PGS000396"]
prs1[,"PGS004512"] = -prs1[,"PGS004512"]
AOU=as.matrix(prs1)
AOU = AOU[,order(colnames(AOU))]
AOU=AOU[,1:32]
d = apply(AOU,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'AOU' = 'value')

total=merge(total,data_cor,by=c("x","y"))

total$Corr=(total$onco+total$NCI+total$TRICL+total$AOU)/4
total=total[,c("x","y","Corr")]

#write.table(total,paste0("/home/sshen/Disk_m2/PRS/plot/32prs_african.txt"),row.names = F,quote=F,sep="\t")


#p.mat <- cor_pmat(d)

list <- rownames(cor)
list <- factor(list,levels = list)


data_cor_tri <- subset(total, match(x, list) >= match(y, list))


tiff(paste0("/home/sshen/Disk_m2/PRS/plot/32prs_african_corr_heatmap.tiff"),width = 12,height = 10,units = "in",res = 300)
p5=ggplot(data_cor_tri,aes(factor(x,levels = list),
                           factor(y,levels = list), 
                           fill=Corr))+  
  geom_tile(colour = "black")+  
  # geom_text(aes(label=Corr),size=2,colour="white")+coord_equal()+
  scale_fill_gradientn(colours = c("white",
                                   brewer.pal(7,"Set1")[1]),na.value = NA,
                       limits=c(0,1),breaks=c(0,0.5,1))+
  labs(x=NULL,y=NULL,title = "32 PRSs correlation in African")+
  scale_y_discrete(position = "right") +
  theme(panel.background=element_blank(),
        axis.text.x = element_text(angle=90,hjust=1,size=12,color="black",face="italic"),
        axis.text.y = element_text(hjust=1,size=12,color="black",face="italic"),
        legend.text = element_text(size = 15),     
        legend.title = element_text(size = 20),   
        axis.ticks.y.right = element_line(), 
        axis.title = element_text(size=20,color="black",face="italic"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        plot.margin = margin(5, 5, 5, 32))
p5
dev.off()


###pearson correlation European#####
setwd("/Users/zhangyixin/Desktop/PRS/trait")
#####Onco_white#####
total=NULL
load("v4_all_trait+LC+LC_meta_prscs/onco_white.RData")
phe=fread("v3_all_trait+LC/pheno_29095.txt")
prs = prs[match(phe$IID,rownames(prs)),]

prs <- data.table(prs)
prs=rena(prs)
onco=prs
d = apply(onco,2,scale)
cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'onco' = 'value')

total=rbind(total,data_cor)

####PLCO_white#####
load("v4_all_trait+LC+LC_meta_prscs/PLCO_white.RData")
phe=fread("v3_all_trait+LC/pheno_98651.txt")
prs = prs[match(phe$IID,rownames(prs)),]

prs <- data.table(prs)
prs=rena(prs)
PLCO=prs
d = apply(PLCO,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'PLCO' = 'value')

total=merge(total,data_cor,by=c("x","y"))
####TRICL_white#####
load("v4_all_trait+LC+LC_meta_prscs/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==1,]
load("v4_all_trait+LC+LC_meta_prscs/TRICL_white.RData")
prs = prs[match(phe$individual_ID,rownames(prs)),]

prs <- data.table(prs)
prs=rena(prs)
TRICL=prs
d = apply(TRICL,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")

data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'TRICL' = 'value')

total=merge(total,data_cor,by=c("x","y"))
####AoU_white#####
load("v4_all_trait+LC+LC_meta_prscs/AOU_white.RData")
load("v4_all_trait+LC+LC_meta_prscs/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eur",]
phe$Lungcancer = phe$lung_ca

prs = prs[match(phe$person_id,rownames(prs)),]

prs <- data.table(prs)
prs=rena(prs)
AOU=prs
d = apply(AOU,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'AOU' = 'value')

total=merge(total,data_cor,by=c("x","y"))

########cor#####
total$Corr=(total$onco+total$PLCO+total$TRICL+total$AOU)/4
total=total[,c("x","y","Corr")]

list <- rownames(cor)
list <- factor(list,levels = list)


data_cor_tri <- subset(total, match(x, list) >= match(y, list))
write.table(data_cor_tri,paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/mtprs_cor_eur.txt"),row.names = F,quote=F,sep="\t")


tiff(paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/trait_eur.tiff"),width = 12,height = 10,units = "in",res = 300)
p5=ggplot(data_cor_tri,aes(factor(x,levels = list),
                           factor(y,levels = list), 
                           fill=Corr))+  
  geom_tile(colour = "black")+  
  # geom_text(aes(label=Corr),size=2,colour="white")+coord_equal()+
  scale_fill_gradientn(colours = c(brewer.pal(7,"Set1")[2],"white",
                                   brewer.pal(7,"Set1")[1]),na.value = NA,
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL,title = "Cross-trait PRSs correlation in European")+
  scale_y_discrete(position = "right") +
  theme(panel.background=element_blank(),
        axis.text.x = element_text(angle=90,hjust=1,size=15,color="black",face="italic"),
        axis.text.y = element_text(hjust=1,size=15,color="black",face="italic"),
        legend.text = element_text(size = 15),     
        legend.title = element_text(size = 15),   
        axis.ticks.y.right = element_line(), 
        axis.title = element_text(size=15,color="black",face="italic"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        plot.margin = margin(5, 5, 5, 32))
p5
dev.off()

###pearson correlation Asian#####
total=NULL
####Onco_asian#####
load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/onco_asian.RData")
phe=fread("/Users/zhangyixin/Desktop/PRS/trait/v3_all_trait+LC/pheno_1974.txt")
prs = prs[match(phe$IID,rownames(prs)),]
tmp=c("F17_ref_afr","F17_ref_eur","fev1_ref_eas","fev1_ref_eur","fev1fvc_ref_eas",
      "fvc_ref_eas","Cigarettes_per_day_ref_eas","Cigarettes_per_day_ref_eur",
      "smoking_cessation_ref_eas","smoking_cessation_ref_eur","smoking_initiation_ref_afr",
      "smoking_initiation_ref_eas","smoking_initiation_ref_eur","Asthma_ref_multi")
for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
onco=prs
d = apply(onco,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'onco' = 'value')


total=rbind(total,data_cor)

####FLCCA_asian#####
load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/FLCCA_asian.RData")
phe=fread("/Users/zhangyixin/Desktop/PRS/trait/v3_all_trait+LC/pheno_8807.txt")
prs = prs[match(phe$IID,rownames(prs)),]
tmp=c("CRP_ref_afr","CRP_ref_eas","emphysema_combine_ref_eur","F17_ref_afr",
      "F17_ref_eas","F17_ref_eur","fev1_ref_afr","fev1_ref_eas","fev1_ref_eur",
      "fev1fvc_ref_afr","fvc_ref_eur","age_of_smoking_ref_afr","age_of_smoking_ref_eas",
      "smoking_cessation_ref_eas","smoking_cessation_ref_eur","smoking_initiation_ref_eas",
      "smoking_initiation_ref_eur","IPF_ref_multi","fvc_ref_afr","fvc_ref_eas")

for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
FLCCA=prs
d = apply(FLCCA,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'FLCCA' = 'value')

total=merge(total,data_cor,by=c("x","y"))

####TRICL_asian#####
load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==6,]

load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/TRICL_white.RData")
prs = prs[match(phe$individual_ID,rownames(prs)),]

tmp=c("CRP_ref_afr","CRP_ref_eas" ,"CRP_ref_eur","emphysema_combine_ref_eur",
      "F17_2_ref_eur","F17_ref_afr","F17_ref_eur","fev1_ref_eur","fev1fvc_ref_eas",
      "fvc_ref_afr","age_of_smoking_ref_afr","age_of_smoking_ref_eas",
      "Cigarettes_per_day_ref_afr","Cigarettes_per_day_ref_eas","Cigarettes_per_day_ref_eur",
      "smoking_cessation_ref_eas","smoking_initiation_ref_afr","smoking_initiation_ref_eas",
      "smoking_initiation_ref_eur","Asthma_ref_multi")

prs=as.data.frame(prs)
for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
TRICL=prs
d = apply(TRICL,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")

data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'TRICL' = 'value')

total=merge(total,data_cor,by=c("x","y"))

####AoU_asian#####
phe=dta
phe=phe[phe$ancestry_pred=="eas",]
phe$Lungcancer = phe$lung_ca

prs = prs[match(phe$person_id,rownames(prs)),]
tmp=c("CRP_ref_afr","CRP_ref_eas","CRP_ref_eur","F17_ref_afr","F17_ref_eur",
      "fev1fvc_ref_eas","fev1fvc_ref_eur","fvc_ref_afr","age_of_smoking_ref_eas",
      "age_of_smoking_ref_eur","Cigarettes_per_day_ref_afr","smoking_cessation_ref_afr",
      "smoking_cessation_ref_eas","Asthma_ref_multi","COPD_ref_multi")

for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
AOU=prs
d = apply(AOU,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'AOU' = 'value')

total=merge(total,data_cor,by=c("x","y"))

########cor#####
total$Corr=(total$onco+total$FLCCA+total$TRICL+total$AOU)/4
total=total[,c("x","y","Corr")]

list <- rownames(cor)
list <- factor(list,levels = list)


data_cor_tri <- subset(total, match(x, list) >= match(y, list))
write.table(data_cor_tri,paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/mtprs_cor_eas.txt"),row.names = F,quote=F,sep="\t")

tiff(paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/trait_eas.tiff"),width = 12,height = 10,units = "in",res = 300)
p5=ggplot(data_cor_tri,aes(factor(x,levels = list),
                           factor(y,levels = list), 
                           fill=Corr))+  
  geom_tile(colour = "black")+  
  # geom_text(aes(label=Corr),size=2,colour="white")+coord_equal()+
  scale_fill_gradientn(colours = c(brewer.pal(7,"Set1")[2],"white",
                                   brewer.pal(7,"Set1")[1]),na.value = NA,
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL,title = "Cross-trait PRSs correlation in Asian")+
  scale_y_discrete(position = "right") +
  theme(panel.background=element_blank(),
        axis.text.x = element_text(angle=90,hjust=1,size=15,color="black",face="italic"),
        axis.text.y = element_text(hjust=1,size=15,color="black",face="italic"),
        legend.text = element_text(size = 15),     
        legend.title = element_text(size = 15),   
        axis.ticks.y.right = element_line(), 
        axis.title = element_text(size=15,color="black",face="italic"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        plot.margin = margin(5, 5, 5, 32))
p5
dev.off()

###pearson correlation African#####
total=NULL
####Onco_african####
load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/onco_black.RData")
phe=fread("/Users/zhangyixin/Desktop/PRS/trait/v3_all_trait+LC/pheno_532.txt")
prs = prs[match(phe$IID,rownames(prs)),]
tmp=c("CRP_ref_eur","F17_2_ref_eur","F17_ref_eur","fev1fvc_ref_afr","fvc_ref_afr",
      "fvc_ref_eur","Cigarettes_per_day_ref_afr","Cigarettes_per_day_ref_eas",
      "Cigarettes_per_day_ref_eur","smoking_cessation_ref_afr","smoking_cessation_ref_eas",
      "smoking_cessation_ref_eur","smoking_initiation_ref_afr","smoking_initiation_ref_eas",
      "smoking_initiation_ref_eur","Asthma_ref_multi")
for(j in tmp){
  prs[,j] = -prs[,j]
}

prs <- data.table(prs)
prs=rena(prs)
onco=prs
d = apply(onco,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'onco' = 'value')


total=rbind(total,data_cor)

####NCI_african#####
load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/NCI_black.RData")
phe=fread("/Users/zhangyixin/Desktop/PRS/trait/v3_all_trait+LC/pheno_5578.txt")
prs = prs[match(phe$IID,rownames(prs)),]
tmp=c("CRP_ref_afr","CRP_ref_eas","CRP_ref_eur","emphysema_combine_ref_eur",
      "F17_2_ref_eur","F17_ref_afr","F17_ref_eur","fev1fvc_ref_eas","fev1fvc_ref_eur",
      "fvc_ref_afr","age_of_smoking_ref_afr","Cigarettes_per_day_ref_afr",
      "Cigarettes_per_day_ref_eas","Cigarettes_per_day_ref_eur","smoking_cessation_ref_afr",
      "smoking_cessation_ref_eas","smoking_cessation_ref_eur","smoking_initiation_ref_afr",
      "smoking_initiation_ref_eas","smoking_initiation_ref_eur","IPF_ref_multi")

for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
NCI=prs
d = apply(NCI,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'NCI' = 'value')

total=merge(total,data_cor,by=c("x","y"))

####TRICL_african#####
load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==5,]

load("/Users/zhangyixin/Desktop/PRS/trait/v4_all_trait+LC+LC_meta_prscs/TRICL_white.RData")
prs = prs[match(phe$individual_ID,rownames(prs)),]
tmp=c("CRP_ref_afr","emphysema_combine_ref_eur","F17_2_ref_eur","F17_ref_afr",
      "F17_ref_eas","F17_ref_eur","fev1_ref_eas","fvc_ref_eas","age_of_smoking_ref_eur",
      "smoking_cessation_ref_afr","smoking_cessation_ref_eas","smoking_cessation_ref_eur",
      "smoking_initiation_ref_afr","Asthma_ref_multi","IPF_ref_multi")

prs=as.data.frame(prs)
for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
TRICL=prs
d = apply(TRICL,2,scale)

cor=cor(d,method = "pearson",use="complete.obs")

data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'TRICL' = 'value')

total=merge(total,data_cor,by=c("x","y"))

####AoU_african#####
phe=dta
phe=phe[phe$ancestry_pred=="afr",]
phe$Lungcancer = phe$lung_ca

prs = prs[match(phe$person_id,rownames(prs)),]
tmp=c("CRP_ref_afr","emphysema_combine_ref_eur","F17_2_ref_eur","F17_ref_eas",
      "fev1_ref_afr","fev1_ref_eur","fev1fvc_ref_eur","fvc_ref_eur","age_of_smoking_ref_eas",
      "age_of_smoking_ref_eur","Cigarettes_per_day_ref_afr","smoking_cessation_ref_eas",
      "smoking_initiation_ref_afr","smoking_initiation_ref_eas")

for(j in tmp){
  prs[,j] = -prs[,j]
}


prs <- data.table(prs)
prs=rena(prs)
AOU=prs
d = apply(AOU,2,scale)


cor=cor(d,method = "pearson",use="complete.obs")
data_cor <- as.data.frame(cor) %>%
  mutate(x = rownames(cor)) %>%
  reshape2::melt(id = 'x') %>%
  rename('y' = 'variable', 'AOU' = 'value')

total=merge(total,data_cor,by=c("x","y"))

########cor#####


total$Corr=(total$onco+total$NCI+total$TRICL+total$AOU)/4
total=total[,c("x","y","Corr")]

list <- rownames(cor)
list <- factor(list,levels = list)


data_cor_tri <- subset(total, match(x, list) >= match(y, list))
write.table(data_cor_tri,paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/mtprs_cor_afr.txt"),row.names = F,quote=F,sep="\t")


tiff(paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/trait_afr.tiff"),width = 12,height = 10,units = "in",res = 300)
p5=ggplot(data_cor_tri,aes(factor(x,levels = list),
                           factor(y,levels = list), 
                           fill=Corr))+  
  geom_tile(colour = "black")+  
  # geom_text(aes(label=Corr),size=2,colour="white")+coord_equal()+
  scale_fill_gradientn(colours = c(brewer.pal(7,"Set1")[2],"white",
                                   brewer.pal(7,"Set1")[1]),na.value = NA,
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL,title = "Cross-trait PRSs correlation in African")+
  scale_y_discrete(position = "right") +
  theme(panel.background=element_blank(),
        axis.text.x = element_text(angle=90,hjust=1,size=15,color="black",face="italic"),
        axis.text.y = element_text(hjust=1,size=15,color="black",face="italic"),
        legend.text = element_text(size = 15),     
        legend.title = element_text(size = 15),   
        axis.ticks.y.right = element_line(), 
        axis.title = element_text(size=15,color="black",face="italic"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        plot.margin = margin(5, 5, 5, 32))
p5
dev.off()

###gene_cor_traits######
prs=fread("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/gene_cor.txt")
prs$label=gsub("emphysema_combine", "Emphysema",prs$label)
prs$label=gsub("F17_2", "Nicotinedependence",prs$label)
prs$label=gsub("crp", "CRP",prs$label)
prs$label=gsub("fev1", "Fev1",prs$label)
prs$label=gsub("fev1fvc", "Fev1fvc",prs$label)
prs$label=gsub("fvc", "Fvc",prs$label)
prs$t1=sapply(prs$label,function(x) strsplit(x,"_")[[1]][1])
prs$t1_race=sapply(prs$label,function(x) strsplit(x,"_")[[1]][2])
prs$t2=sapply(prs$label,function(x) strsplit(x,"_")[[1]][3])
prs$t2_race=sapply(prs$label,function(x) strsplit(x,"_")[[1]][4])
prs$x=paste0(prs$t1,"_",prs$t1_race)
prs$y=paste0(prs$t2,"_",prs$t2_race)

prs$x=gsub("AgeSmk", "Age_of_smoking",prs$x)
prs$y=gsub("AgeSmk", "Age_of_smoking",prs$y)
prs$x=gsub("CigDay", "Cigarettes_per_day",prs$x)
prs$y=gsub("CigDay", "Cigarettes_per_day",prs$y)
prs$x=gsub("Nicotinedependence_eur", "Nicotine_dependence_eur",prs$x)
prs$y=gsub("Nicotinedependence_eur", "Nicotine_dependence_eur",prs$y)
prs$x=gsub("F17", "Tobacco_disorder",prs$x)
prs$y=gsub("F17", "Tobacco_disorder",prs$y)
prs$x=gsub("SmkCes", "Smoking_cessation",prs$x)
prs$y=gsub("SmkCes", "Smoking_cessation",prs$y)
prs$x=gsub("SmkInit", "Smoking_initiation",prs$x)
prs$y=gsub("SmkInit", "Smoking_initiation",prs$y)
a=unique(prs$x)
b=a
all=NULL
for(i in a){
  for(j in b){
    tmp=data.frame(i,j)
    all=rbind(all,tmp)
  }
}


prs1=prs
prs1$y=prs$x
prs1$x=prs$y
total=rbind(prs,prs1)
total=total[,c("x","y","Cor","P")]
total=unique(total)
total=merge(total,all,by.x=c("x","y"),by.y=c("i","j"),all=T)

total$Cor=ifelse(total$Cor>1 |total$Cor< -1,3,total$Cor)
total$Cor=ifelse(is.na(total$Cor)==T,3,total$Cor)
total$P=ifelse(total$Cor==3,0,total$P)
total$Cor[total$x==total$y]=1
total$P[total$x==total$y]=0

list <- unique(c(total$x, total$y))
data_cor_tri <- subset(total, match(x, list) >= match(y, list))

# 绘图
p5 <- ggplot(data_cor_tri, aes(factor(x, levels = list), 
                               factor(y, levels = list), 
                               fill = Cor)) +
  geom_tile(data = subset(data_cor_tri, Cor=3),
            fill = "grey70", colour = "black") +
  geom_tile(colour = "black") +
  geom_text(aes(label = ifelse(P > 0.05, "*", "")), color = "black", size = 5) +
  coord_equal() +
  scale_fill_gradientn(colours = c(brewer.pal(7,"Set1")[2], "white", brewer.pal(7,"Set1")[1]),
                       na.value = NA, limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1)) +
  labs(x = NULL, y = NULL, title = "Traits genetic correlation") +
  scale_y_discrete(position = "right") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black", face = "italic"),
        axis.text.y = element_text(hjust = 1, size = 15, color = "black", face = "italic"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.ticks.y.right = element_line(),
        axis.title = element_text(size = 10, color = "black", face = "italic"),
        plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 32))

p5
ggsave(paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/Figure/genetic_correlation.tiff"),width = 12,height = 10,units = "in",dpi = 300)

dev.off()

