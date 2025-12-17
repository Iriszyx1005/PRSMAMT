library(ggplot2)
require(data.table)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(ggcorrplot)
library(corrplot)
library(ggbreak)
library(openxlsx)
library(patchwork)
####Circle Diagram
require(CMplot)
setwd("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/prsma系数")
AoU_ref_eur = fread("EUR_AoU_TRICL_OncoArray.txt")
AoU_ref_eur=AoU_ref_eur[,c("V1","V3")]
colnames(AoU_ref_eur)=c("SNP","AoU_ref_eur")

AoU_ref_eas = fread("EAS_AoU_TRICL_OncoArray.txt")
AoU_ref_eas=AoU_ref_eas[,c("V1","V3")]
colnames(AoU_ref_eas)=c("SNP","AoU_OncoArray_TRICL_ref_eas")

AoU_ref_afr = fread("AFR_AoU_TRICL_OncoArray.txt")
AoU_ref_afr=AoU_ref_afr[,c("V1","V3")]
colnames(AoU_ref_afr)=c("SNP","AoU_OncoArray_TRICL_ref_afr")

FLCCA_ref_eas = fread("EAS_FLCCA.txt")
FLCCA_ref_eas=FLCCA_ref_eas[,c("V1","V3")]
colnames(FLCCA_ref_eas)=c("SNP","FLCCA_ref_eas")

NCI_ref_afr = fread("AFR_NCI.txt")
NCI_ref_afr=NCI_ref_afr[,c("V1","V3")]
colnames(NCI_ref_afr)=c("SNP","NCI_ref_afr")

OncoArray_TRICL_ref_eur = fread("EUR_TRICL_OncoArray.txt")
OncoArray_TRICL_ref_eur=OncoArray_TRICL_ref_eur[,c("V1","V3")]
colnames(OncoArray_TRICL_ref_eur)=c("SNP","OncoArray_TRICL_ref_eur")

PLCO_ref_eur = fread("EUR_PLCO.txt")
PLCO_ref_eur=PLCO_ref_eur[,c("V1","V3")]
colnames(PLCO_ref_eur)=c("SNP","PLCO_ref_eur")


prs = merge(AoU_ref_eur,OncoArray_TRICL_ref_eur,by="SNP",all=T)
prs = merge(prs,PLCO_ref_eur,by="SNP",all=T)
prs = merge(prs,AoU_ref_eas,by="SNP",all=T)
prs = merge(prs,FLCCA_ref_eas,by="SNP",all=T)
prs = merge(prs,AoU_ref_afr,by="SNP",all=T)
prs = merge(prs,NCI_ref_afr,by="SNP",all=T)



prs$Chromosome = sapply(prs$SNP, function(x) strsplit(x,":")[[1]][1])
prs$Chromosome = as.numeric(gsub("chr","",prs$Chromosome))
prs$Position = as.numeric(sapply(prs$SNP, function(x) strsplit(x,":")[[1]][2]))

prs[, 2:8] <- lapply(prs[, 2:8], function(x) ifelse(is.na(x), 0, x))

prs = prs[,c(1,9,10,2:8)]
prs = data.frame(prs)
colnames(prs)=c("SNP","Chromosome","Position","AoU(eur)","OncoArray / TRICL(eur)","PLCO(eur)",
                "AoU / OncoArray / TRICL(eas)", "FLCCA(eas)","AoU / OncoArray / TRICL(afr)",
                "NCI(afr)")
#write.xlsx(prs,"prs_effect_7situation1.xlsx",rowNames = F,quote=F,sep="\t")
save(prs,file="prs_effect_7situation.RData")
write.table(prs, "prs_effect_7situation.txt", sep = "\t", row.names = FALSE, quote = FALSE)

load("prs_effect_7situation.RData")
for(i in 4:10) {
  prs[,i] = abs(prs[,i])
  prs[,i] = 10^(-prs[,i])
}

setwd("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit")
tiff(paste0("Figure/prs_Circular1.tiff"),width = 12,height = 12,units = "in",res = 300)

CMplot(prs,plot.type="c",r=1,cex = 0.6, col=matrix(c("#51999F","#7BC0CD","#BFDFD2","#DBCB92","#ECB66C","#EA9E58","#F4C700"),nrow=7),
       threshold.col = c("#F46D43","#559FCD"),cir.chr.h = 2,
       signal.cex = c(0.8,0.8), ylim=list(c(0,0.01),c(0,0.01),c(0,0.01),c(0,0.01),c(0,0.01),c(0,0.01),c(0,0.01)), signal.col=c("#F46D43","#4D9221"),file.output = F)
legend(
  "topright",
  legend = c("AoU(eur)","OncoArray / TRICL(eur)","PLCO(eur)",
             "AoU / OncoArray / TRICL(eas)", "FLCCA(eas)","AoU / OncoArray / TRICL(afr)",
             "NCI(afr)"),
  col = c("#51999F","#7BC0CD","#BFDFD2","#DBCB92","#ECB66C","#EA9E58","#F4C700"),
  pch = 19,
  pt.cex = 1.5,
  bty = "n",
  ncol = 1
)

dev.off()
