rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(pROC))
setwd("/home/sshen/Disk_m2/PRS/UKB")

transform_intensity <- function(x) {
  ((x / 10)^(-1)) - 0.4021541613
}

# model
plcom2012_risk <- function(age, race, education, bmi, copd, personal_cancer,
                           family_history, smoking_status, intensity,
                           duration, quit_time) {
  #Centralization
  age <- age - 62
  education <- education - 4
  bmi <- bmi - 27
  duration <- duration - 27
  quit_time <- quit_time - 10
  intensity <- transform_intensity(intensity)
  
  # Race Corresponding Regression Coefficient
  race_beta <- c(
    "White" = 0,
    "Black" = 0.3944778,
    "Hispanic" = -0.7434744,
    "Asian" = -0.466585
    # "Pacific" = 0,
    # "American" = 1.027152
  )
  race_coef <- race_beta[race]
  
  if (is.na(race_coef)) race_coef <- 0
  
  logit <- age * 0.0778868 +
    race_coef +
    education * -0.0812744 +
    bmi * -0.0274194 +
    copd * 0.3553063 +
    personal_cancer * 0.4589971 +
    family_history * 0.587185 +
    smoking_status * 0.2597431 +
    intensity * (-1.822606) +
    duration * 0.0317321 +
    quit_time * (-0.0308572) +
    -4.532506
  
  risk <- exp(logit) / (1 + exp(logit))
  return(risk)
}

load("~/UKB/ukb20240117.RData")
cov = fread("UKB/ukb_smoking.csv")
cov = cov[match(ukb$id_shen,cov$`Participant ID`),-1]
colnames(cov)
cov[cov == ""] <- NA

matrix_merge = function(x,numeric = F){
  x1 = x[,1]
  for(i in 2:(ncol(x))){x1[is.na(x1)] =  x[,i][is.na(x1)]}
  if (numeric) x1 = as.numeric(x1)
  return(x1)
}
cpd1 = matrix_merge(as.matrix(cov[,1:4]),numeric = T)
cpd2 = matrix_merge(as.matrix(cov[,5:8]),numeric = T)
age_smoke1 = matrix_merge(as.matrix(cov[,9:12]),numeric = T)
age_smoke2 = matrix_merge(as.matrix(cov[,13:16]),numeric = T)
age_stop = matrix_merge(as.matrix(cov[,17:20]),numeric = T)
edu = matrix_merge(as.matrix(cov[,21:24]),numeric = F)

cpd = cpd1;cpd[is.na(cpd)] =  cpd2[is.na(cpd)]
age_smoke = age_smoke1;age_smoke[is.na(age_smoke)] =  age_smoke2[is.na(age_smoke)]

edu1= NA
edu1[grep("nursing",edu)] = 1
edu1[grep("NVQ",edu)] = 2
edu1[grep("CSEs",edu)] = 3
edu1[grep("GCSEs",edu)] = 4
edu1[grep("A levels",edu)] = 5
edu1[grep("College",edu)] = 6
table(edu1,useNA = "ifany")

ukb$copd_emphysema_cbronchitis=ifelse(ukb$emphysema_combine==1 |ukb$copd_combine | ukb$cbronchitis==1,1,0)
ukb$family_history=ifelse(ukb$history_siblings==3 |
                            ukb$history_father==3 |
                            ukb$history_mother==3,1,0)
ukb$smoking_status=ifelse(ukb$smoke==2,1,0)
ukb$personal_cancer=0
ukb$smkyears=ifelse(ukb$smoke==2,(ukb$age-age_smoke),ifelse(ukb$smoke==1,(age_stop-age_smoke),0))
ukb$smkyears[ukb$smkyears<0] = 0

ukb$race = NA
ukb$race[ukb$ethnic %in% c("1001","1002","1003","1")] = "White"
ukb$race[ukb$ethnic %in% c("3001","3002","3003","3004","3")] = "Asian"
ukb$race[ukb$ethnic %in% c("4001","4002","4003","4")] = "Black"


ukb_plco = data.frame(ukb[,c("age","race","bmi","copd_emphysema_cbronchitis","personal_cancer","family_history","smoking_status","smkyears","quityr","smoke","lung_ca","personyear","id_shen")],edu1,cpd)
ukb_plco = subset(ukb_plco,smoke!=0)
IID=ukb_plco$id_shen
ukb_plco=select(ukb_plco,-c("id_shen"))

imp_mice = function(x,select=1,max=3){
  if(!require(mice)){install.packages("mice");require(mice)}
  tmp1 <- mice(x,m=max,maxit=10,meth='pmm',seed=500)
  tmp2 <- complete(tmp1,select)
  return(tmp2)
}
ukb_plco = imp_mice(ukb_plco)

ukb_plco$plco <- mapply(plcom2012_risk,
                        age = ukb_plco$age,
                        race = ukb_plco$race,
                        education = ukb_plco$edu1,
                        bmi = ukb_plco$bmi,
                        copd = ukb_plco$copd_emphysema_cbronchitis,
                        personal_cancer = ukb_plco$personal_cancer,
                        family_history = ukb_plco$family_history,
                        smoking_status = ukb_plco$smoking_status,
                        intensity = ukb_plco$cpd,
                        duration = ukb_plco$smkyears,
                        quit_time = ukb_plco$quityr)

df=cbind(IID,ukb_plco)
roc(df$lung_ca,df$plco)

save(df,file="ukb_plcom2012_0522.RData")
