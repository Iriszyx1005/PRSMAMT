library(survival) 
library(survminer)
library(coxphf)
library(ggplot2)
library(pROC)
library(data.table)
library(tidyr)
library(dplyr)
library(openxlsx)
setwd("/home/sshen/Disk_m2/PRS/UKB")
load("UKB.RData")
load("UKB_PRSMA.RData")


load("~/UKB/ukb20240117.RData")
ukb = ukb[match(prs$IID,ukb$id_shen),]
score = score[match(prs$IID,score$IID),]
pROC::roc(ukb$lung_ca,prs$PRSMA)

mt = as.matrix(prs[,2:33])

load("/home/sshen/Disk_m2/PRS/trait_PRS_data/final_PRS_eur.RData")
x1 = x1[,order(colnames(x1))]
mt = mt[,order(colnames(mt))]
x1 = x1[,c(1:12,17:26,13,27:32,14:16)]
colnames(x1) = colnames(mt)

library(catboost)
id1 = "PLCO"
train_pool <- catboost.load_pool(data = x1[label==id1,], label = y[label==id1])
test_pool <- catboost.load_pool(data = mt, label = ukb$lung_ca)

params <- list(
  iterations = 5000, 
  learning_rate = 0.001,  
  depth = 5,   
  l2_leaf_reg = 8,  
  loss_function = 'Logloss', 
  eval_metric = 'AUC',  
  random_seed = 0929   
)

model <- catboost.train(train_pool, test_pool, params)

prs$PRSMT<- catboost.predict(model, test_pool, prediction_type = "Probability")


prs1=merge(prs,ukb,by.x="IID",by.y="id_shen")
prs1$Lungcancer=prs1$lung_ca

save(prs1,file="ukb_prsmamt_0522.RData")

######evaluation of prsmamt
setwd("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit")
load("result/UKB/UKB_prscsx.RData")
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
save(prs1,file="result/UKB/ukb_prsmamt_1116.RData")

#####compare with PGS004691#####
load("result/UKB/ukb_prsmamt_1116.RData")
prs = subset(prs1,personyear > 0)
rm(prs1)

prs$MAMT <- (prs$MAMT - mean(prs$MAMT[prs$Lungcancer==0]))/sd(prs$MAMT[prs$Lungcancer==0])
mamt_quant <- quantile(prs$MAMT, probs = seq(0, 1, by = 0.2))
prs$prsmamt_cat <- cut(prs$MAMT, breaks = mamt_quant, include.lowest = TRUE,labels = paste0("Q", 1:5))
prs$prsmamt_cat <- relevel(prs$prsmamt_cat, ref = "Q3")


prs$PGS004691 <- (prs$PGS004691 - mean(prs$PGS004691[prs$Lungcancer==0]))/sd(prs$PGS004691[prs$Lungcancer==0])
pgs4691_quant <- quantile(prs$PGS004691, probs = seq(0, 1, by = 0.2))
prs$pgs4691_cat <- cut(prs$PGS004691, breaks = pgs4691_quant, include.lowest = TRUE,labels = paste0("Q", 1:5))
prs$pgs4691_cat <- relevel(prs$pgs4691_cat, ref = "Q3")

re=NULL
decile=c(1,2,4,5)
df=prs

model = coxph(Surv(personyear,lung_ca) ~ factor(prsmamt_cat), data = df);summary(model)
cl=exp(confint(model)[1:4,1]);effect=exp(coef(model)[1:4]);cu=exp(confint(model)[1:4,2]);beta=summary(model)$coefficients[,1];se=summary(model)$coefficients[,3]
re=rbind(re,data.frame(race="ALL",PRS="PRSMAMT",metric="HR",cl,effect,cu,beta,se,decile))


model = coxph(Surv(personyear,lung_ca) ~ factor(pgs4691_cat), data = df);summary(model)
cl=exp(confint(model)[1:4,1]);effect=exp(coef(model)[1:4]);cu=exp(confint(model)[1:4,2]);beta=summary(model)$coefficients[,1];se=summary(model)$coefficients[,3]
re=rbind(re,data.frame(race="ALL",PRS="PGS004691",metric="HR",cl,effect,cu,beta,se,decile))
re=rbind(re,data.frame(race="ALL",PRS="PRSMAMT",metric="HR",cl=1,effect=1,cu=1,beta=0,se=NA,decile=3))
re=rbind(re,data.frame(race="ALL",PRS="PGS004691",metric="HR",cl=1,effect=1,cu=1,beta=0,se=NA,decile=3))

re1=re
rlt=NULL
for(i in c(1,2,4,5)){
  d=re1[re1$decile==i,]
  beta1 <- d$beta[d$PRS=="PGS004691"]
  se1   <- d$se[d$PRS=="PGS004691"]
  beta2 <- d$beta[d$PRS=="PRSMAMT"]
  se2   <- d$se[d$PRS=="PRSMAMT"]
  
  Z = (beta1 - beta2) / sqrt(se1^2 + se2^2)
  Z_p = 2 * pnorm(-abs(Z))
  
  #Q_p = metagen(TE = beta, seTE = se, sm = "HR",data=x)$pval.Q
  tmp=data.frame(decile=i,Z_p)
  rlt=rbind(rlt,tmp)
}
rlt$sig = ifelse(rlt$Z_p < 0.001, "***",
                 ifelse(rlt$Z_p>= 0.001  & rlt$Z_p < 0.01, "**",
                        ifelse(rlt$Z_p >= 0.01 & rlt$Z_p < 0.05, "*","")))


re1 <- left_join(re1, rlt, by = "decile")

re1=re1[order(re1$PRS,re1$decile),]
re1$PRS <- factor(re1$PRS, levels = c("PGS004691","PRSMAMT"))
custom_colors <- c( "PRSMAMT" ="#B8474D","PGS004691" = "#4E6691")

re1$terms <- factor(re1$decile, levels = paste0(1:5), 
                    labels = c("<20%", "20-40%", "40-60%", "60-80%","80-100%"))
#re1=re1[order(re1$score,re1$decile),]

write.xlsx(re1,paste0("result/UKB/ukb_prs_hr.xlsx"),rowNames = F,quote=F,sep="\t")

bracket_df <- data.frame(
  x_start = c(4.8),
  x_end   = c(5.2),
  y       = rep(2.2, 1)  
)

legend_labels <- c(
  "PRSMAMT" = expression(paste(PRS["MAMT"])),
  "PGS004691" = "PGS004691")


#tiff(paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/ukb_prs_hr1.tiff"),width = 15,height = 10,units = "in",res = 300)

p=ggplot(re1, aes(x = terms, y = effect, group = PRS, color = PRS)) +
  #geom_ribbon(aes(ymin = CL, ymax = CU, fill = PRS), alpha = 0.2, color = NA) +
  geom_point(size=5,pch = 19,position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = cl, ymax = cu), width = 0.2, linewidth = 1,
                position = position_dodge(width = 0.6)) + 
  scale_y_continuous(breaks = seq(0, 2, by = 0.5), limits = c(0, 2.3)) +
  geom_text(aes(y=2.2,label = sig),
            color = "black", size = 10) +
  geom_segment(data = bracket_df,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               inherit.aes = FALSE,
               linewidth = 0.8, color = "black")+
  geom_hline(yintercept = 1, linetype = "dashed",color = "#B9636A", size = 1) +
  scale_color_manual(values = custom_colors,labels = legend_labels) + 
  labs(x = "PRS Quantile", y = "Hazard Ratio (HR)", title = "HR Trend by PRS Quantiles for All Ancestries") +
  theme(panel.background=element_blank(),
        axis.title = element_text(size=25,color="black"),
        plot.title = element_text(size=25,face = "bold.italic",hjust=0.5),
        axis.text = element_text(size=25,color="black"),
        legend.text = element_text(size = 23),  
        legend.title = element_text(size = 23),  
        axis.line = element_line(color="black" ,linewidth = 0.5),
        axis.ticks.length = unit(-1,"mm"),
        axis.text.x = element_text(hjust = 0.5))
p
ggsave(paste0("Figure/ukb/ukb_prs_hr.pdf"),p,width = 15,height = 10,units = "in",dpi = 300)

dev.off()


#####compare with plco#####
load("/Users/zhangyixin/Desktop/PRS/ukb/ukb_plcom2012_0522.RData")
roc(df$lung_ca,df$plco)
df=df[,c("IID","plco")]

load("result/UKB/ukb_prsmamt_1116.RData")
prs=merge(df,prs1,by="IID")
rm(df,prs1)
prs=prs[prs$personyear>0,]

roc(prs$Lungcancer,prs$plco)
roc(prs$Lungcancer,prs$MAMT)

prs$ma_z <- as.numeric(scale(prs$PRSMA))
prs$mt_z <- as.numeric(scale(prs$PRSMT))
prs$mamt_z <- as.numeric(scale(prs$MAMT))
prs$plco_z <- as.numeric(scale(prs$plco))

model = glm(lung_ca ~ ma_z+plco_z, data = prs);summary(model)
prs$plco_prsma = predict(model)

model = glm(lung_ca ~ mt_z+plco_z, data = prs);summary(model)
prs$plco_prsmt = predict(model)

model = glm(lung_ca ~ mamt_z+plco_z, data = prs);summary(model)
prs$final = predict(model)

plco_quant <- quantile(prs$plco, probs = seq(0, 1, by = 0.2))
prs$plco_cat <- cut(prs$plco, breaks = plco_quant, include.lowest = TRUE,labels = paste0("Q",1:5))
prs$plco_cat <- relevel(prs$plco_cat, ref = "Q3")

plco_prsma_quant <- quantile(prs$plco_prsma, probs = seq(0, 1, by = 0.2))
prs$plco_prsma_cat <- cut(prs$plco_prsma, breaks = plco_prsma_quant, include.lowest = TRUE,labels = paste0("Q",1:5))
prs$plco_prsma_cat <- relevel(prs$plco_prsma_cat, ref = "Q3")

plco_prsmt_quant <- quantile(prs$plco_prsmt, probs = seq(0, 1, by = 0.2))
prs$plco_prsmt_cat <- cut(prs$plco_prsmt, breaks = plco_prsmt_quant, include.lowest = TRUE,labels = paste0("Q",1:5))
prs$plco_prsmt_cat <- relevel(prs$plco_prsmt_cat, ref = "Q3")

final_quant <- quantile(prs$final, probs = seq(0, 1, by = 0.2))
prs$final_cat <- cut(prs$final, breaks = final_quant, include.lowest = TRUE,labels = paste0("Q",1:5))
prs$final_cat <- relevel(prs$final_cat, ref = "Q3")

mamt_quant <- quantile(prs$MAMT, probs = seq(0, 1, by = 0.1))

re=NULL
#terms <- paste0("Q", 2:5)
decile=c(1,2,4,5)
model = coxph(Surv(personyear,lung_ca) ~ plco_cat, data = prs);summary(model)
tmp=cbind(summary(model)$coefficients,exp(confint(model)))
colnames(tmp)=c("beta","HR","se","z","P","CL","CU")
re=rbind(re,data.frame(beta=0,HR=1,se=NA,z=NA,P=NA,CL=1,CU=1,score="PLCOm2012",decile=3))
re=rbind(re,data.frame(tmp,score="PLCOm2012",decile))

model = coxph(Surv(personyear,lung_ca) ~ factor(plco_prsma_cat), data = prs);summary(model)
tmp=cbind(summary(model)$coefficients,exp(confint(model)))
colnames(tmp)=c("beta","HR","se","z","P","CL","CU")
re=rbind(re,data.frame(beta=0,HR=1,se=NA,z=NA,P=NA,CL=1,CU=1,score="PLCOm2012+PRSMA",decile=3))
re=rbind(re,data.frame(tmp,score="PLCOm2012+PRSMA",decile))

model = coxph(Surv(personyear,lung_ca) ~ factor(plco_prsmt_cat), data = prs);summary(model)
tmp=cbind(summary(model)$coefficients,exp(confint(model)))
colnames(tmp)=c("beta","HR","se","z","P","CL","CU")
re=rbind(re,data.frame(beta=0,HR=1,se=NA,z=NA,P=NA,CL=1,CU=1,score="PLCOm2012+PRSMT",decile=3))
re=rbind(re,data.frame(tmp,score="PLCOm2012+PRSMT",decile))

model = coxph(Surv(personyear,lung_ca) ~ final_cat, data = prs);summary(model)
tmp=cbind(summary(model)$coefficients,exp(confint(model)))
colnames(tmp)=c("beta","HR","se","z","P","CL","CU")
re=rbind(re,data.frame(beta=0,HR=1,se=NA,z=NA,P=NA,CL=1,CU=1,score="PLCOm2012+PRSMAMT",decile=3))
re=rbind(re,data.frame(tmp,score="PLCOm2012+PRSMAMT",decile))

re1=re[re$score=="PLCOm2012" |re$score=="PLCOm2012+PRSMAMT",]
re1=re

rlt=NULL
for(i in 1:5){
  d=re1[re1$decile==i,]
  beta1 <- d$beta[d$score=="PLCOm2012"]
  se1   <- d$se[d$score=="PLCOm2012"]
  beta2 <- d$beta[d$score=="PLCOm2012+PRSMAMT"]
  se2   <- d$se[d$score=="PLCOm2012+PRSMAMT"]
  
  Z = (beta1 - beta2) / sqrt(se1^2 + se2^2)
  Z_p = 2 * pnorm(-abs(Z))
  
  #Q_p = metagen(TE = beta, seTE = se, sm = "HR",data=x)$pval.Q
  tmp=data.frame(decile=i,Z_p)
  rlt=rbind(rlt,tmp)
}
rlt$sig = ifelse(rlt$Z_p < 0.001, "***",
                 ifelse(rlt$Z_p>= 0.001  & rlt$Z_p < 0.01, "**",
                        ifelse(rlt$Z_p >= 0.01 & rlt$Z_p < 0.05, "*","")))


re1 <- left_join(re1, rlt, by = "decile")


re1$score <- factor(re1$score, levels = c("PLCOm2012","PLCOm2012+PRSMAMT"))
custom_colors <- c( "PLCOm2012+PRSMAMT" ="#B8474D","PLCOm2012" = "#4E6691")

re1$score <- factor(re1$score, levels = c("PLCOm2012","PLCOm2012+PRSMA","PLCOm2012+PRSMT","PLCOm2012+PRSMAMT"))
custom_colors <- c( "PLCOm2012+PRSMAMT" ="#B8474D","PLCOm2012" = "#4E6691",
                    "PLCOm2012+PRSMA" ="#C39EBB","PLCOm2012+PRSMT" ="#A0522D90")

re1$terms <- factor(re1$decile, levels = paste0(1:5), 
                    labels = c("<20%", "20-40%", "40-60%", "60-80%","80-100%"))
re1=re1[order(re1$score,re1$decile),]

write.xlsx(re1,paste0("result/UKB/ukb_plco_prs3_5.xlsx"),rowNames = F,quote=F,sep="\t")

bracket_df <- data.frame(
  x_start = c(1.76,4.76),
  x_end   = c(2.2,5.2),
  y       = rep(9,2)  
)


legend_labels <- c(
  "PLCOm2012" = expression(paste(PLCO["m2012"])),
  "PLCOm2012+PRSMAMT" = expression(paste(PLCO["m2012"],"+",PRS["MAMT"])))

legend_labels <- c(
  "PLCOm2012" = expression(paste(PLCO["m2012"])),
  "PLCOm2012+PRSMA" = expression(paste(PLCO["m2012"],"+",PRS["MA"])),
  "PLCOm2012+PRSMT" = expression(paste(PLCO["m2012"],"+",PRS["MT"])),
  "PLCOm2012+PRSMAMT" = expression(paste(PLCO["m2012"],"+",PRS["MAMT"])))


#tiff(paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/ukb_plco_hr.tiff"),width = 15,height = 10,units = "in",res = 300)

p=ggplot(re1, aes(x = terms, y = HR, group = score, color = score)) +
  #geom_ribbon(aes(ymin = CL, ymax = CU, fill = PRS), alpha = 0.2, color = NA) +
  scale_y_continuous(breaks = seq(0, 9, by = 1), limits = c(0, 9)) +
  geom_point(size=5,pch = 19,position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.2, linewidth = 1,
                position = position_dodge(width = 0.6)) + 
  geom_text(aes(y=9,label = sig),color = "black", size = 10) +
  geom_hline(yintercept = 1, linetype = "dashed",color = "#B9636A", size = 1) +
  geom_segment(data = bracket_df,aes(x = x_start, xend = x_end, y = y, yend = y),inherit.aes = FALSE,linewidth = 0.8, color = "black")+
  scale_color_manual(values = custom_colors,labels = legend_labels) +  
  labs(x = "Quantile", y = "Hazard Ratio (HR)", title = "HR Trend by PRS Quantiles for All Ancestries") +
  theme(panel.background=element_blank(),
        axis.title = element_text(size=25,color="black"),
        plot.title = element_text(size=25,face = "bold.italic",hjust=0.5),
        axis.text = element_text(size=25,color="black"),
        legend.text = element_text(size = 23),  
        legend.title = element_text(size = 23),  
        axis.line = element_line(color="black" ,linewidth = 0.5),
        axis.ticks.length = unit(-1,"mm"),
        axis.text.x = element_text(hjust = 0.5))
p
ggsave(paste0("Figure/ukb/ukb_plco_prsmamt_5.pdf"),width = 15,height = 10,units = "in",dpi = 300)

dev.off()

#####Risk Stratification#####
plco_quant <- quantile(prs$plco, probs = seq(0, 1, by = 0.1))
mamt_quant <- quantile(prs$MAMT, probs = seq(0, 1, by = 0.1))

prs$score_d = rep(0,nrow(prs))
prs$score_d[prs$plco<plco_quant[9] & prs$MAMT<mamt_quant[9]] = 1
prs$score_d[prs$plco<plco_quant[9] & prs$MAMT>=mamt_quant[9]] = 2
prs$score_d[prs$plco>=plco_quant[9] & prs$MAMT<mamt_quant[9]] = 3
prs$score_d[prs$plco>=plco_quant[9] & prs$MAMT>=mamt_quant[9]] = 4

model = coxph(Surv(personyear,lung_ca) ~ factor(score_d), data = prs);summary(model)
#HR:1.81 6.57. 12.64
require(survminer)
fit=survfit(Surv(personyear,lung_ca)~ factor(score_d),data = prs)
legend_labels <- c(
  "factor(score_d)=1" = expression(paste("Low ",PLCO["m2012"], " at low genetic risk (reference)")),
  "factor(score_d)=2" = expression(paste("Low ",PLCO["m2012"], " at high genetic risk (HR = 1.82)")),
  "factor(score_d)=3" = expression(paste("High ",PLCO["m2012"], " at low genetic risk (HR = 6.57)")),
  "factor(score_d)=4" = expression(paste("High ",PLCO["m2012"], " at high genetic risk (HR = 12.64)")))

p <- ggsurvplot(
  fit,
  data = prs,
  conf.int = FALSE,
  ylim = c(0, 0.2),
  fun = "cumhaz",
  risk.table = FALSE,
  palette = "lancet",
  legend = c(0.4, 0.85),
  legend.title = "",
  xlab = "Follow-up years",
  ylab = "Cumulative risk",
  ggtheme = theme(panel.background = element_blank(),axis.line  = element_line(),
                  axis.title = element_text(size=20,color="black"),
                  plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
                  axis.text = element_text(size=20,color="black"),
                  legend.text = element_text(size = 20))
)
p$plot <- p$plot +
  scale_color_manual(values = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF"),
                     labels = legend_labels) +
  scale_linetype_manual(values = c(1,1,1,1),labels = legend_labels)
p
ggsave(paste0("Figure/ukb/Cumulative_risk.pdf"),width = 10,height =8,units = "in",dpi = 300)

dev.off()

cut_q = 0.7
score_plco = ifelse(prs$plco< quantile(prs$plco,cut_q),0,1)
score_mamt1 = ifelse(prs$MAMT<quantile(prs$MAMT,cut_q),0,1)

table(score_plco,score_mamt1,prs$lung_ca)
table(prs$lung_ca,score_mamt1)
table(prs$lung_ca,score_plco)

###NRI####
library(nricens)
library(Hmisc)

prs$logit_plco <- log(prs$plco / (1 - prs$plco))

fit <- glm(lung_ca ~ logit_plco + MAMT, data = prs,family = binomial)
prs$risk_mamt  <- predict(fit, type = "response")

fit1 <-  glm(lung_ca ~ plco, data = prs, family = binomial)
prs$pred1 <- predict(fit1, type = "response")

risk_cut <- c(0.01,0.05)

nri_result <- nribin(
  event   = prs$lung_ca,
  p.std   = prs$pred1,   # old model
  p.new  = prs$risk_mamt,    # new model
  cut     =risk_cut, 
  niter = 0,
  updown  = 'category' 
)

#####DCA####
library(rmda)
prs$logit_plco <- log(prs$plco / (1 - prs$plco))

fit <- glm(lung_ca ~ logit_plco + MAMT, data = prs,family = binomial)
#fit_mamt <- glm(lung_ca ~ mamt_z, data = prs,family = binomial)
prs$risk_mamt  <- predict(fit, type = "response")


fit1 <-  glm(lung_ca ~ plco, data = prs, family = binomial)

prs$pred1 <- predict(fit1, type = "response")

dca1 <- decision_curve(lung_ca ~ pred1, data = prs, family = binomial)
dca2 <- decision_curve(lung_ca ~ risk_mamt, data = prs, family = binomial)

tiff(paste0("Figure/ukb/DCA.tiff"),width = 6,height = 6,units = "in",res = 300)

plot_decision_curve(
  list(dca1, dca2),
  curve.names = c("PLCOm2012", "PLCOm2012+PRSMAMT"),
  col = c("#4E6691","#B8474D"),
  xlim =c(0.01,0.15),ylim=c(0,0.6),
  lwd = c(3, 2, 2, 1),
  legend.position = "topright"
)
dev.off()
