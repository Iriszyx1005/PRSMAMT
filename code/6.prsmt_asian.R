require(data.table)
library(pROC)

setwd("/home/sshen/Disk_m2/PRS")
#### Asian
covar_onco_eas = fread("../PRS_yxzhang/train/asian/pheno_1974.txt")
load("trait_PRS_data/onco_asian.RData")
prs = prs[match(covar_onco_eas$IID,rownames(prs)),];prs[,1:29] = -prs[,1:29]
prs_onco_eas = apply(prs,2,scale)

covar_FLCCA = fread("../PRS_yxzhang/validation/FLCCA/pheno_8807.txt")
load("trait_PRS_data/FLCCA_asian.RData")
prs = prs[match(covar_FLCCA$IID,rownames(prs)),];prs[,1:29] = -prs[,1:29]
prs_FLCCA = apply(prs,2,scale)

## AOU
load("AOU/multi_trait_score.RData");load("AOU/LungCa_V3_eas.RData")
score = cbind(score,LC_ref_multi_v5 = prs$score)
load("AOU/dta_final.RData");covar_aou = dta
prs_aou = score[match(covar_aou$person_id,rownames(score)),]
temp = colnames(score)
temp = gsub("_AFR","_ref_afr",temp);temp = gsub("_EAS","_ref_eas",temp);temp = gsub("_EUR","_ref_eur",temp)
temp = gsub("Asthma","Asthma_ref_multi",temp);temp = gsub("COPD","COPD_ref_multi",temp);temp = gsub("IPF","IPF_ref_multi",temp)
colnames(prs_aou) = temp

prs_aou = prs_aou[,match(colnames(prs_onco_eas),colnames(prs_aou))]

# tf = dta$ancestry_pred=="eas"
# summary(glm(lung_ca_combine~prs_aou[tf,"LungCa_V2_eas"],family = binomial,data = dta[tf,]))

tf = dta$ancestry_pred=="eas"
prs_aou1 = prs_aou[tf,]
prs_aou1 = scale(prs_aou1)
covar_aou1 = covar_aou[tf,]


load("TRICL_multiethnic.RData")
covar_TRICL = pheno
covar_TRICL$Lungcancer = covar_TRICL$disease
load("trait_PRS_data/TRICL_white.RData")
prs = prs[match(covar_TRICL$individual_ID,rownames(prs)),];prs[,1:29] = -prs[,1:29]
temp = fread("LungCa_V3/ALL_TRICL_asian_prs.txt")
prs[,35] = temp$PRS[match(covar_TRICL$individual_ID,temp$IID)]
prs = apply(prs,2,scale)
tf = covar_TRICL$ethnicity1==6
covar_TRICL1 = covar_TRICL[tf,]
prs_TRICL1 = prs[tf,]


y = c(covar_onco_eas$Lungcancer,covar_FLCCA$Lungcancer,covar_TRICL1$Lungcancer,covar_aou1$lung_ca_combine)
x = rbind(prs_onco_eas,prs_FLCCA,prs_TRICL1,prs_aou1)
x1 = x[,1:32]
label = c(rep("onco",nrow(covar_onco_eas)),rep("FLCCA",nrow(covar_FLCCA)),rep("TRICL",nrow(prs_TRICL1)),rep("aou",nrow(prs_aou1)))
# x2 = x1[label %in% c("onco","TRICL"),]

# for(i in 1:ncol(x)){
#     print(colnames(x)[i])
#     print(summary(glm(y~x[,i],family = binomial)))
#     print(roc(covar_aou1$lung_ca_combine,prs_aou1[,i])$auc)
# }

ORSD = function(y,x,ci=F){
  model = glm(y~scale(x),family="binomial")
  rlt = exp(coef(model)[2])
  if(ci){
    cl = exp(coef(model)[2]- 1.96*coef(summary(model))[2,2])
    cu = exp(coef(model)[2]+ 1.96*coef(summary(model))[2,2])      
    rlt = c(cl,rlt,cu)
  }
  return(rlt)
}


# model = glm(y~x[,c(8,11,20,29,31)],family = binomial);summary(model)
# y1 = predict(model)


library(grpreg)
set.seed(123)
group <- c(1,1,1,2,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,13,14)


# Cross-validation: best lambda
cv_fit <- cv.grpreg(x1, y, group, penalty = "grLasso",family="binomial")
print(cv_fit$lambda.min)
beta <- coef(cv_fit)[-1]
print(beta);y1 = x1 %*% beta

# fit <- grpreg(x1, y, group, penalty = "grLasso",family="binomial",lambda = 0.0025)
# fit_se = select(fit, "AIC")
# beta = fit_se$beta[-1];print(beta)
# y1 = x1 %*% beta
x2 = x1[,beta!=0]

y0 = x[,35]

library(catboost)

id1 = "FLCCA";id2 = "TRICL"
train_pool <- catboost.load_pool(data = x2[label==id1,], label = y[label==id1])
test_pool <- catboost.load_pool(data = x2[label==id2,], label = y[label==id2])
all_pool <- catboost.load_pool(data = x2, label = y)
#### look for best parameter
param_space = NULL
for(iterations in  c(100, 500, 1000,5000)){
  for(learning_rate in  c(0.001, 0.01 , 0.1, 0.2)){
    for(depth in c(3, 5, 10)){
      for(l2_leaf_reg in c(1,3,8)){
        temp = data.frame(iterations,learning_rate,depth,l2_leaf_reg)
        param_space = rbind(param_space,temp)
      }}}}

rlt = NULL
for (i in 1:nrow(param_space)) {
  param_cur <- param_space[i,]
  
  params <- list(
    iterations = param_cur$iterations, 
    learning_rate = param_cur$learning_rate,  
    depth = param_cur$depth,  
    l2_leaf_reg = param_cur$l2_leaf_reg, 
    loss_function = 'Logloss', 
    eval_metric = 'AUC',  
    random_seed = 12      
  )
  
  model <- catboost.train(train_pool, test_pool, params)
  y1 <- catboost.predict(model, all_pool, prediction_type = "Probability")
  
  id = "onco";
  y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
  AUC1 = roc(y[label==id],y2)$auc;OR1 = ORSD(y[label==id],y2)
  id = "TRICL";
  y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
  AUC2 = roc(y[label==id],y2)$auc;OR2 = ORSD(y[label==id],y2)
  id = "aou";
  y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
  AUC3 = roc(y[label==id],y2)$auc;OR3 = ORSD(y[label==id],y2)
  
  temp = data.frame(i, iterations = param_cur$iterations, learning_rate = param_cur$learning_rate, depth = param_cur$depth,AUC1,AUC2,AUC3,OR1,OR2,OR3)
  rlt = rbind(rlt,temp)
}
rlt$AUC_sum = rlt$AUC1 + rlt$AUC2 + rlt$AUC3 
rlt = rlt[order(-rlt$AUC_sum),]
rlt$l2_leaf_reg = param_space$l2_leaf_reg[match(rlt$i,rownames(param_space))]
write.csv(rlt,"catboost_rlt_asian.csv",row.names=F,quote=F)


params <- list(
  iterations = 500,         
  learning_rate = 0.1,       
  depth = 5,                
  loss_function = 'Logloss',
  l2_leaf_reg = 3,
  eval_metric = 'AUC',  
  random_seed = 12            
)
model <- catboost.train(train_pool, test_pool, params)
y1 <- catboost.predict(model, all_pool, prediction_type = "Probability")

id = "onco";roc(y[label==id],y1[label==id])$auc;ORSD(y[label==id],y1[label==id])
id = "FLCCA";roc(y[label==id],y1[label==id])$auc;ORSD(y[label==id],y1[label==id])
id = "TRICL";roc(y[label==id],y1[label==id])$auc;ORSD(y[label==id],y1[label==id])
id = "aou";roc(y[label==id],y1[label==id])$auc;ORSD(y[label==id],y1[label==id])

id = "onco";
y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
roc(y[label==id],y2)$auc;ORSD(y[label==id],y2)
id = "FLCCA";
y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
roc(y[label==id],y2)$auc;ORSD(y[label==id],y2)
id = "TRICL";
y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
roc(y[label==id],y2)$auc;ORSD(y[label==id],y2)
id = "aou";
y2 = predict(glm(y[label==id]~y1[label==id]+y0[label==id],family = binomial))
roc(y[label==id],y2)$auc;ORSD(y[label==id],y2)


id = "onco";roc(y[label==id],y0[label==id])$auc;ORSD(y[label==id],y0[label==id])
id = "FLCCA";roc(y[label==id],y0[label==id])$auc;ORSD(y[label==id],y0[label==id])
id = "TRICL";roc(y[label==id],y0[label==id])$auc;ORSD(y[label==id],y0[label==id])
id = "aou";roc(y[label==id],y0[label==id])$auc;ORSD(y[label==id],y0[label==id])

save(y,y0,y1,label,x1,x2,file = "trait_PRS_data/final_PRS_eas.RData")

