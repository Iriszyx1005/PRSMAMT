library(data.table)
library(dplyr)
library(pROC)
library(openxlsx)
####EAS####
#####FLCCA####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/FLCCA_eas_prscsx.RData")
prs$IID=rownames(prs)
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_8807.txt")
prs = prs[match(phe$IID,rownames(prs)),]
y_FLCCA=phe$Lungcancer
prs_FLCCA=prs[,c("IID","LC_prscsx")]

#####TRICL####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/TRICL_prscsx.RData")
prs$IID=rownames(prs)
load("/Users/zhangyixin/Desktop/PRS/pheno/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==6,]
prs=prs[match(phe$individual_ID,rownames(prs)),]
y_TRICL=phe$Lungcancer
prs_TRICL=prs[,c("IID","LC_prscsx")]

#####OncoArray#####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/oncoarray_eas_prscsx.RData")
prs$IID=rownames(prs)
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_1974.txt")
prs = prs[match(phe$IID,rownames(prs)),]
y_OncoArray=phe$Lungcancer
prs_OncoArray=prs[,c("IID","LC_prscsx")]

#####AoU####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/AOU_prscsx.RData")
prs$IID=rownames(prs)
load("/Users/zhangyixin/Desktop/PRS/pheno/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eas",]
prs = prs[match(phe$person_id,prs$IID),]
rm(dta)
y_AoU=phe$lung_ca
prs_AoU=prs[,c("IID","LC_prscsx")]

###############
load("/Users/zhangyixin/Desktop/PRS/Bayesian_Network/final_PRS_eas.RData")
re=data.frame(label,y1)
unique(re$label)
#"onco"  "FLCCA" "TRICL" "aou" 
y_total=c(y_OncoArray,y_FLCCA,y_TRICL,y_AoU)
prs_total=rbind(prs_OncoArray,prs_FLCCA,prs_TRICL,prs_AoU)
re=cbind(re,y_total,prs_total)
colnames(re)
colnames(re)=c("Cohort","PRSMT","Lungcancer","IID","PRSMA")
re=re[,c("Lungcancer","PRSMA","PRSMT","IID","Cohort")]
re$Cohort=ifelse(re$Cohort=="onco","OncoArray",ifelse(re$Cohort=="aou","AoU",re$Cohort))
unique(re$Cohort)
re=re[is.na(re$PRSMA)==F,]
result=re
save(result,file = "/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_eas.RData")

#"OncoArray" "FLCCA"     "TRICL"     "AoU"
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_eas.RData")
total = NULL
id_list = c("TRICL","OncoArray","AoU")
for(i in id_list){
  re = NULL
  prs1=result[result$Cohort==i,]
  prs1=prs1[is.na(prs1$PRSMA)==F,]
  prs1$PRSMA <- (prs1$PRSMA - mean(prs1$PRSMA[prs1$Lungcancer==0]))/sd(prs1$PRSMA[prs1$Lungcancer==0])
  prs1$PRSMT <- (prs1$PRSMT - mean(prs1$PRSMT[prs1$Lungcancer==0]))/sd(prs1$PRSMT[prs1$Lungcancer==0])
  
  model = glm(Lungcancer~PRSMA+PRSMT,family = binomial,data = prs1);summary(model)
  prs1$MAMT= coef(model)[2]*prs1$PRSMA+coef(model)[3]*prs1$PRSMT
  prs1$MAMT <- (prs1$MAMT - mean(prs1$MAMT[prs1$Lungcancer==0]))/sd(prs1$MAMT[prs1$Lungcancer==0])
  
  roc(prs1$Lungcancer,prs1$PRSMA)
  
  roc(prs1$Lungcancer,prs1$PRSMT)
  
  roc(prs1$Lungcancer,prs1$MAMT)
  ##PRSMA
  model = glm(Lungcancer~prs1$PRSMA,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$PRSMA,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  ##PRSMT
  model = glm(Lungcancer~prs1$PRSMT,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$PRSMT,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  
  ####PRSMAMT
  model = glm(Lungcancer~prs1$MAMT,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$MAMT,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  
  rownames(re)=NULL
  re=re[order(re$Metric),]  
  total=rbind(total,re)
}

total$Cl=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Cu,total$Cl)
total$Cu=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Cl,total$Cu)
total$Effect=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Effect,total$Effect)

total$Cl=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Cu,total$Cl)
total$Cu=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Cl,total$Cu)
total$Effect=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Effect,total$Effect)

write.xlsx(total,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_perform_eas.xlsx"),rowNames = F,quote=F,sep="\t")

####AFR####
#####NCI####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/NCI_afr_prscsx.RData")
prs$IID=rownames(prs)
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_5578.txt")
prs = prs[match(phe$IID,rownames(prs)),]
y_NCI=phe$Lungcancer
prs_NCI=prs[,c("IID","LC_prscsx")]

#####TRICL####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/TRICL_prscsx.RData")
prs$IID=rownames(prs)
load("/Users/zhangyixin/Desktop/PRS/pheno/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==5,]
prs=prs[match(phe$individual_ID,rownames(prs)),]
y_TRICL=phe$Lungcancer
prs_TRICL=prs[,c("IID","LC_prscsx")]

#####OncoArray#####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/oncoarray_afr_prscsx.RData")
prs$IID=rownames(prs)
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_532.txt")
prs = prs[match(phe$IID,rownames(prs)),]
y_OncoArray=phe$Lungcancer
prs_OncoArray=prs[,c("IID","LC_prscsx")]

#####AoU####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/AOU_prscsx.RData")
prs$IID=rownames(prs)
load("/Users/zhangyixin/Desktop/PRS/pheno/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="afr",]
prs = prs[match(phe$person_id,rownames(prs)),]
rm(dta)
y_AoU=phe$lung_ca
prs_AoU=prs[,c("IID","LC_prscsx")]

###############
load("/Users/zhangyixin/Desktop/PRS/Bayesian_Network/final_PRS_afr.RData")
re=data.frame(label,y1)
unique(re$label)
#"onco"  "NCI"   "TRICL" "aou"  
y_total=c(y_OncoArray,y_NCI,y_TRICL,y_AoU)
prs_total=rbind(prs_OncoArray,prs_NCI,prs_TRICL,prs_AoU)
re=cbind(re,y_total,prs_total)
colnames(re)
colnames(re)=c("Cohort","PRSMT","Lungcancer","IID","PRSMA")
re=re[,c("Lungcancer","PRSMA","PRSMT","IID","Cohort")]
re$Cohort=ifelse(re$Cohort=="onco","OncoArray",ifelse(re$Cohort=="aou","AoU",re$Cohort))
unique(re$Cohort)
re=re[is.na(re$PRSMA)==F,]
result=re
save(result,file = "/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_afr.RData")

#OncoArray" "NCI"       "TRICL"     "AoU" 
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_afr.RData")
total = NULL
id_list = c("TRICL","OncoArray","AoU")
for(i in id_list){
  re = NULL
  prs1=result[result$Cohort==i,]
  prs1=prs1[is.na(prs1$PRSMA)==F,]
  prs1$PRSMA <- (prs1$PRSMA - mean(prs1$PRSMA[prs1$Lungcancer==0]))/sd(prs1$PRSMA[prs1$Lungcancer==0])
  prs1$PRSMT <- (prs1$PRSMT - mean(prs1$PRSMT[prs1$Lungcancer==0]))/sd(prs1$PRSMT[prs1$Lungcancer==0])
  
  model = glm(Lungcancer~PRSMA+PRSMT,family = binomial,data = prs1);summary(model)
  prs1$MAMT= coef(model)[2]*prs1$PRSMA+coef(model)[3]*prs1$PRSMT
  prs1$MAMT <- (prs1$MAMT - mean(prs1$MAMT[prs1$Lungcancer==0]))/sd(prs1$MAMT[prs1$Lungcancer==0])
  
  roc(prs1$Lungcancer,prs1$PRSMA)
  
  roc(prs1$Lungcancer,prs1$PRSMT)
  
  roc(prs1$Lungcancer,prs1$MAMT)
  ##PRSMA
  model = glm(Lungcancer~prs1$PRSMA,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$PRSMA,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  ##PRSMT
  model = glm(Lungcancer~prs1$PRSMT,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$PRSMT,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  
  ####PRSMAMT
  model = glm(Lungcancer~prs1$MAMT,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$MAMT,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  
  rownames(re)=NULL
  re=re[order(re$Metric),]  
  total=rbind(total,re)
}

total$Cl=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Cu,total$Cl)
total$Cu=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Cl,total$Cu)
total$Effect=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Effect,total$Effect)

total$Cl=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Cu,total$Cl)
total$Cu=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Cl,total$Cu)
total$Effect=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Effect,total$Effect)

write.xlsx(total,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_perform_afr.xlsx"),rowNames = F,quote=F,sep="\t")

####EUR####
#####PLCO####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/PLCO_eur_prscsx.RData")
prs$IID=rownames(prs)
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_98651.txt")
prs = prs[match(phe$IID,rownames(prs)),]
y_PLCO=phe$Lungcancer
prs_PLCO=prs[,c("IID","LC_prscsx")]

#####TRICL####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/TRICL_prscsx.RData")
prs$IID=rownames(prs)
load("/Users/zhangyixin/Desktop/PRS/pheno/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==1,]
prs=prs[match(phe$individual_ID,rownames(prs)),]
y_TRICL=phe$Lungcancer
prs_TRICL=prs[,c("IID","LC_prscsx")]

#####OncoArray#####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/oncoarray_eur_prscsx.RData")
prs$IID=rownames(prs)
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_29095.txt")
prs = prs[match(phe$IID,rownames(prs)),]
y_OncoArray=phe$Lungcancer
prs_OncoArray=prs[,c("IID","LC_prscsx")]

#####AoU####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/AOU_prscsx.RData")
prs$IID=rownames(prs)
load("/Users/zhangyixin/Desktop/PRS/pheno/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eur",]
prs = prs[match(phe$person_id,rownames(prs)),]
rm(dta)
y_AoU=phe$lung_ca
prs_AoU=prs[,c("IID","LC_prscsx")]

###############
load("/Users/zhangyixin/Desktop/PRS/Bayesian_Network/final_PRS_eur.RData")
re=data.frame(label,y1)
unique(re$label)
#"onco"  "PLCO"  "TRICL" "aou" 
y_total=c(y_OncoArray,y_PLCO,y_TRICL,y_AoU)
prs_total=rbind(prs_OncoArray,prs_PLCO,prs_TRICL,prs_AoU)
re=cbind(re,y_total,prs_total)
colnames(re)
colnames(re)=c("Cohort","PRSMT","Lungcancer","IID","PRSMA")
re=re[,c("Lungcancer","PRSMA","PRSMT","IID","Cohort")]
re$Cohort=ifelse(re$Cohort=="onco","OncoArray",ifelse(re$Cohort=="aou","AoU",re$Cohort))
unique(re$Cohort)
re=re[is.na(re$PRSMA)==F,]
result=re
save(result,file = "/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_eur.RData")

#OncoArray" "PLCO"       "TRICL"     "AoU" 
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_eur.RData")
total = NULL
id_list = c("TRICL","OncoArray","AoU")
for(i in id_list){
  re = NULL
  prs1=result[result$Cohort==i,]
  prs1=prs1[is.na(prs1$PRSMA)==F,]
  prs1$PRSMA <- (prs1$PRSMA - mean(prs1$PRSMA[prs1$Lungcancer==0]))/sd(prs1$PRSMA[prs1$Lungcancer==0])
  prs1$PRSMT <- (prs1$PRSMT - mean(prs1$PRSMT[prs1$Lungcancer==0]))/sd(prs1$PRSMT[prs1$Lungcancer==0])
  
  model = glm(Lungcancer~PRSMA+PRSMT,family = binomial,data = prs1);summary(model)
  prs1$MAMT= coef(model)[2]*prs1$PRSMA+coef(model)[3]*prs1$PRSMT
  prs1$MAMT <- (prs1$MAMT - mean(prs1$MAMT[prs1$Lungcancer==0]))/sd(prs1$MAMT[prs1$Lungcancer==0])
  
  roc(prs1$Lungcancer,prs1$PRSMA)
  
  roc(prs1$Lungcancer,prs1$PRSMT)
  
  roc(prs1$Lungcancer,prs1$MAMT)
  ##PRSMA
  model = glm(Lungcancer~prs1$PRSMA,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$PRSMA,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  ##PRSMT
  model = glm(Lungcancer~prs1$PRSMT,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$PRSMT,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  
  ####PRSMAMT
  model = glm(Lungcancer~prs1$MAMT,family = binomial,data = prs1);summary(model)
  proc = roc(prs1$Lungcancer,prs1$MAMT,ci = TRUE)
  re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="AUC",
                         Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
  re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="OR/SD",
                         Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
  
  rownames(re)=NULL
  re=re[order(re$Metric),]  
  total=rbind(total,re)
}

total$Cl=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Cu,total$Cl)
total$Cu=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Cl,total$Cu)
total$Effect=ifelse(total$Metric=="AUC" & total$Effect<0.5,1-total$Effect,total$Effect)

total$Cl=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Cu,total$Cl)
total$Cu=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Cl,total$Cu)
total$Effect=ifelse(total$Metric=="OR/SD" & total$Effect<1,1/total$Effect,total$Effect)

write.xlsx(total,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_perform_eur.xlsx"),rowNames = F,quote=F,sep="\t")

###R2#####
re=NULL
###EUR####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_eur.RData")
race="EUR"
######TRICL####
study="TRICL"
load("/Users/zhangyixin/Desktop/PRS/pheno/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==1,]
phe=phe[,c("individual_ID","age_age","gender","typesmok")]
colnames(phe)=c("IID","age","sex","typesmok")

######OncoArray####
study="OncoArray"
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_29095.txt")
phe=phe[,c("IID","age","sex","typesmok")]

#####AoU####
study="AoU"
load("/Users/zhangyixin/Desktop/PRS/pheno/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eur",]
phe=phe[,c("person_id","age","sex","smoke")]
colnames(phe)=c("IID","age","sex","typesmok")

###EAS####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_eas.RData")
race="EAS"
######TRICL####
study="TRICL"
load("/Users/zhangyixin/Desktop/PRS/pheno/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==6,]
phe=phe[,c("individual_ID","age_age","gender","typesmok")]
colnames(phe)=c("IID","age","sex","typesmok")

######OncoArray####
study="OncoArray"
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_1974.txt")
phe=phe[,c("IID","age","sex","typesmok")]

#####AoU####
study="AoU"
load("/Users/zhangyixin/Desktop/PRS/pheno/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="eas",]
phe=phe[,c("person_id","age","sex","smoke")]
colnames(phe)=c("IID","age","sex","typesmok")

###AFR####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_PRS_afr.RData")
race="AFR"
######TRICL####
study="TRICL"
load("/Users/zhangyixin/Desktop/PRS/pheno/TRICL_multiethnic.RData")
phe=pheno
phe$Lungcancer = phe$disease
phe=phe[phe$Lungcancer!= -9] 
phe = phe[phe$ethnicity1==5,]
phe=phe[,c("individual_ID","age_age","gender","typesmok")]
colnames(phe)=c("IID","age","sex","typesmok")

######OncoArray####
study="OncoArray"
phe=fread("/Users/zhangyixin/Desktop/PRS/pheno/pheno_532.txt")
phe=phe[,c("IID","age","sex","typesmok")]

#####AoU####
study="AoU"
load("/Users/zhangyixin/Desktop/PRS/pheno/dta_final.RData")
phe=dta
phe=phe[phe$ancestry_pred=="afr",]
phe=phe[,c("person_id","age","sex","smoke")]
colnames(phe)=c("IID","age","sex","typesmok")



prs1=result[result$Cohort==study,]
prs1=merge(prs1,phe,by="IID")
colnames(prs1)
prs1$smoking_status=ifelse(prs1$typesmok==-9 | prs1$typesmok==4,NA,prs1$typesmok)
prs1$PRSMA <- (prs1$PRSMA - mean(prs1$PRSMA[prs1$Lungcancer==0]))/sd(prs1$PRSMA[prs1$Lungcancer==0])
prs1$PRSMT <- (prs1$PRSMT - mean(prs1$PRSMT[prs1$Lungcancer==0]))/sd(prs1$PRSMT[prs1$Lungcancer==0])

model = glm(Lungcancer~PRSMA+PRSMT,family = binomial,data = prs1);summary(model)
prs1$MAMT= coef(model)[2]*prs1$PRSMA+coef(model)[3]*prs1$PRSMT
prs1$MAMT <- (prs1$MAMT - mean(prs1$MAMT[prs1$Lungcancer==0]))/sd(prs1$MAMT[prs1$Lungcancer==0])

m0=glm(Lungcancer~age+factor(sex)+factor(smoking_status),family = binomial,data = prs1);summary(m0)
m1 = glm(Lungcancer~age+factor(sex)+factor(smoking_status)+MAMT,family = binomial,data = prs1)
delta=fmsb::NagelkerkeR2(m1)$R2- fmsb::NagelkerkeR2(m0)$R2
re=rbind(re,data.frame(Cohort=study,Ancestry=race,
                       Model="Lungcancer ~ age+sex+smoking_status",
                       NagelkerkeR2=fmsb::NagelkerkeR2(m0)$R2,
                       Incremental_R2=""))
re=rbind(re,data.frame(Cohort=study,Ancestry=race,
                       Model="Lungcancer ~ age+sex+smoking_status+PRSMAMT",
                       NagelkerkeR2=fmsb::NagelkerkeR2(m1)$R2,
                       Incremental_R2=delta))

write.xlsx(re,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/Incremental_R2.xlsx"),rowNames = F,quote=F,sep="\t")


###UKB####
load("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/UKB/UKB_prscsx.RData")
new=prs[,c("IID","PRSMA")]
load("/Users/zhangyixin/Desktop/PRS/ukb/ukb_prsmamt_0522.RData")
prs=prs[,-c("PRSMA","MAMT")]
prs=merge(prs,new,by="IID")
i="UKB"
prs1=prs
re=NULL
prs1=prs1[is.na(prs1$PRSMA)==F,]

model = glm(Lungcancer~PRSMA+PRSMT,family = binomial,data = prs1);summary(model)
prs1$MAMT=  predict(model)
# prs1$MAMT <- (prs1$MAMT - mean(prs1$MAMT[prs1$Lungcancer==0]))/sd(prs1$MAMT[prs1$Lungcancer==0])

roc(prs1$Lungcancer,prs1$PRSMA)

roc(prs1$Lungcancer,prs1$PRSMT)

roc(prs1$Lungcancer,prs1$MAMT)
save(prs1,file="/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/UKB/ukb_prsmamt_1116.RData")

##PRSMA
model = glm(Lungcancer~prs1$PRSMA,family = binomial,data = prs1);summary(model)
proc = roc(prs1$Lungcancer,prs1$PRSMA,ci = TRUE)
re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="AUC",
                       Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
re=rbind(re,data.frame(Study=i,PRS="PRSMA",Metric="OR/SD",
                       Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))
##PRSMT
model = glm(Lungcancer~prs1$PRSMT,family = binomial,data = prs1);summary(model)
proc = roc(prs1$Lungcancer,prs1$PRSMT,ci = TRUE)
re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="AUC",
                       Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
re=rbind(re,data.frame(Study=i,PRS="PRSMT",Metric="OR/SD",
                       Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))

####PRSMAMT
model = glm(Lungcancer~prs1$MAMT,family = binomial,data = prs1);summary(model)
proc = roc(prs1$Lungcancer,prs1$MAMT,ci = TRUE)
re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="AUC",
                       Cl=proc$ci[1],Effect=proc$auc,Cu=proc$ci[3]))
re=rbind(re,data.frame(Study=i,PRS="PRSMAMT",Metric="OR/SD",
                       Cl=exp(confint(model)[2,1]),Effect=exp(coef(model)[2]),Cu=exp(confint(model)[2,2])))

rownames(re)=NULL
re=re[order(re$Metric),]  

write.xlsx(re,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/UKB/final_perform_ukb.xlsx"),rowNames = F,quote=F,sep="\t")

#####R2#####
re=read.xlsx("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/Incremental_R2.xlsx")
study="UKB"
race="Multi"
prs1$smoking_status=ifelse(prs1$smoke==2,1,0)

m0=glm(Lungcancer~age+factor(gender)+factor(smoking_status),family = binomial,data = prs1);summary(m0)
m1 = glm(Lungcancer~age+factor(gender)+factor(smoking_status)+MAMT,family = binomial,data = prs1)
delta=fmsb::NagelkerkeR2(m1)$R2- fmsb::NagelkerkeR2(m0)$R2
re=rbind(re,data.frame(Cohort=study,Ancestry=race,
                       Model="Lungcancer ~ age+sex+smoking_status",
                       NagelkerkeR2=fmsb::NagelkerkeR2(m0)$R2,
                       Incremental_R2=""))
re=rbind(re,data.frame(Cohort=study,Ancestry=race,
                       Model="Lungcancer ~ age+sex+smoking_status+PRSMAMT",
                       NagelkerkeR2=fmsb::NagelkerkeR2(m1)$R2,
                       Incremental_R2=delta))

write.xlsx(re,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/Incremental_R2.xlsx"),rowNames = F,quote=F,sep="\t")

####plco#####
re=read.xlsx("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/Incremental_R2.xlsx")
study="UKB"
race="Multi"

m0=glm(Lungcancer~plco,family = binomial,data = prs);summary(m0)
m1 = glm(Lungcancer~plco+MAMT,family = binomial,data = prs);summary(m1)
delta=fmsb::NagelkerkeR2(m1)$R2- fmsb::NagelkerkeR2(m0)$R2
re=rbind(re,data.frame(Cohort=study,Ancestry=race,
                       Model="Lungcancer ~ PLCO",
                       NagelkerkeR2=fmsb::NagelkerkeR2(m0)$R2,
                       Incremental_R2=""))
re=rbind(re,data.frame(Cohort=study,Ancestry=race,
                       Model="Lungcancer ~ PLCO+PRSMAMT",
                       NagelkerkeR2=fmsb::NagelkerkeR2(m1)$R2,
                       Incremental_R2=delta))

write.xlsx(re,paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/Incremental_R2.xlsx"),rowNames = F,quote=F,sep="\t")



