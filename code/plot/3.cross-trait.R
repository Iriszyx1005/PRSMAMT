require(data.table)
require(CMplot)
library(dplyr)
require(qqman) 
library(openxlsx)
library(ggplot2)
require(qqman) 
library(ggrepel)
library(ggbreak)
library(patchwork)

####cross-trait heatmap#####
require(data.table)
setwd("/home/sshen/Disk_m2/PRS")
require(openxlsx)
heat = read.xlsx("plot/AUC_OR_P.xlsx",sheet = "AUC1")
heat1 = as.matrix(heat[,-1])
for(i in 1:ncol(heat1)){
  temp = heat1[,i]
  temp = ifelse(temp<0.5,1-temp,temp)
  heat1[,i] = temp
}
rownames(heat1) = heat$Trait

library(pheatmap);library(grid)
colnames(heat1) <- gsub("_white", "_European", colnames(heat1))
colnames(heat1) <- gsub("_asian", "_Asian", colnames(heat1))
colnames(heat1) <- gsub("_black", "_African", colnames(heat1))
colnames(heat1) <- gsub("Onco_", "OncoArray_", colnames(heat1))
tiff('plot/pheatmap_multitrait1.tiff',width = 11, height =20 ,units = "in",res = 300)
pheatmap(heat1,fontsize=15, fontsize_row=23,fontsize_col=23,cluster_cols = F,cluster_rows = F,show_colnames = T,color = colorRampPalette(c("white","red"))(100),breaks = seq(0.5,0.6,length=100), border_color = "grey10",display_numbers = T,angle_col =90)
dev.off()


label = c("Age_of_smoking_afr","Age_of_smoking_eas","Age_of_smoking_eur","Cigarettes_per_day_afr","Cigarettes_per_day_eas","Cigarettes_per_day_eur","CRP_afr","CRP_eas","CRP_eur","Fev1_afr","Fev1_eas","Fev1_eur","Fev1fvc_afr","Fev1fvc_eas","Fev1fvc_eur","Fvc_afr","Fvc_eas","Fvc_eur","Smoking_cessation_afr","Smoking_cessation_eas","Smoking_cessation_eur","Smoking_initiation_afr","Smoking_initiation_eas","Smoking_initiation_eur","Tobacco_disorder_afr","Tobacco_disorder_eas","Tobacco_disorder_eur")
temp = apply(heat1[,9:12],1,mean)
temp = temp[names(temp) %in% label]
mean(temp[grep("afr",names(temp))])
mean(temp[grep("eas",names(temp))])
summary(temp[grep("eur",names(temp))])

dta_plot = data.frame(afr = temp[grep("afr",names(temp))],eas = temp[grep("eas",names(temp))],eur = temp[grep("eur",names(temp))])
dta_plot = reshape2::melt(dta_plot)
library(ggplot2)
ggplot(dta_plot,aes(variable,value))+
  stat_summary(mapping=aes(fill = variable),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7)+
  stat_summary(mapping=aes(fill = variable),fun=mean,geom = "point",fun.args = list(mult=1),size = 4)+
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+
  geom_jitter(aes(fill = variable),position = position_jitter(0.2),shape=21, size = 2,alpha=0.9)+
  labs(x = "Ancestry",y = "AUC")+
  scale_y_continuous(limits = c(0.5,0.55))+
  scale_fill_manual(values = c(
    "eur" = "#559FCD",
    "eas" = "#7FBC41",
    "afr" = "#F4C700"
  )) +
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25), 
        axis.line=element_line(colour="black",size=0.25), 
        axis.title=element_text(size=13,color="black"), 
        axis.text = element_text(size=12,color="black"), 
        legend.position="none"
  )
ggsave(paste0("plot/cross_trait_afr1.tiff"),width = 3,height =4,units = "in",dpi = 300)
