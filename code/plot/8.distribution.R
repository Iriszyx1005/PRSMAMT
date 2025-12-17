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
library(tidyr)

########Figure S1: PRS testing set#####
rows <- c("OncoArray_European","OncoArray_Asian","OncoArray_African", 
          "PLCO_European", "FLCCA_Asian", "NCI_African", 
          "TRICL_European","TRICL_Asian", "TRICL_African",
          "AoU_European","AoU_Asian","AoU_African")
cols <- c("INTEGRAL-ILCCO", "MVP", "FinnGen", "UKB", "PLCO", "NCI", "FLCCA", "CKB", "BBJ")

df <- expand.grid(PRS = rows, Column = cols, stringsAsFactors = FALSE) %>%
  mutate(value = ifelse(PRS == Column, 0, 1))

df$value[df$PRS=="OncoArray_European" & df$Column=="INTEGRAL-ILCCO"]=0
df$value[df$PRS=="OncoArray_Asian" & df$Column=="INTEGRAL-ILCCO"]=0
df$value[df$PRS=="OncoArray_African" & df$Column=="INTEGRAL-ILCCO"]=0
df$value[df$PRS=="TRICL_European" & df$Column=="INTEGRAL-ILCCO"]=0
df$value[df$PRS=="TRICL_Asian" & df$Column=="INTEGRAL-ILCCO"]=0
df$value[df$PRS=="TRICL_African" & df$Column=="INTEGRAL-ILCCO"]=0
df$value[df$PRS=="PLCO_European" & df$Column=="PLCO"]=0
df$value[df$PRS=="FLCCA_Asian" & df$Column=="FLCCA"]=0
df$value[df$PRS=="NCI_African" & df$Column=="NCI"]=0

display.brewer.all()  
all_colors <- brewer.pal.info["Blues", "maxcolors"]
selected_colors <- brewer.pal(all_colors, "Set1") 
subset_colors <- selected_colors[c(1:2)] 
colorRampPalette(subset_colors)(20)

tiff(paste0("/Users/zhangyixin/Desktop/PRS/z_script/Figure/meta_scource.tiff"),width = 10.05,height = 10.05,units = "in",res = 300)
ggplot(df, aes(x = Column, y = PRS, fill = factor(value))) +
  geom_tile(color = "black",linewidth=0.5) +
  scale_fill_manual(values = c("0" ="#BF2F3C80" , "1" = "white"),
                    name = "Datasets",
                    labels = c("Hold-out testing","PRS training")) +
  labs(x = "GWAS meta-analysis scources", y = "PRS testing sets") +
  theme_minimal(base_size = 12) +
  theme( panel.background = element_rect(fill = "transparent", color = NA),
         plot.background = element_rect(fill = "transparent", color = NA),
         axis.text.x = element_text(angle = 45, hjust = 1, size = 20, face = "bold"),
         axis.text.y = element_text(size = 20, face = "bold"),
         axis.title = element_text(size = 22, face = "bold"),
         legend.text = element_text(size = 15, face = "bold"),
         legend.title = element_text(size = 15, face = "bold"),
         panel.grid = element_blank(),
         axis.ticks = element_blank(),
         legend.position = "top",
         legend.direction = "horizontal" 
  )
#ggsave("/Users/zhangyixin/Desktop/PRS/z_result/cor/meta scource.png", plot = p, width = 9,height = 10, bg = "transparent")
dev.off()

###Figure S4: distribution_prscsx####OncoArray TRICL AoU PLCO FLCCA NCI
study="PLCO"
k="eur"
j="European"
load(paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/PRSMA_MT/final_",study,"_",k,"_prscsx.RData"))
data=prs1
data$PRSMA <- (data$PRSMA - mean(data$PRSMA[data$Lungcancer==0]))/sd(data$PRSMA[data$Lungcancer==0])
p <-ggplot(data,aes(x=PRSMA,fill = as.factor(Lungcancer), color = as.factor(Lungcancer))) +  
  geom_density(alpha = 0.6) + 
  theme(panel.background = element_blank(),axis.line  = element_line(),
        axis.title = element_text(size=20,color="black"),
        plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        axis.text = element_text(size=20,color="black"),
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20), 
        legend.position.inside=c(0.2,0.8),  plot.margin = margin(5, 5, 5, 30)) +
  xlab(bquote(PRS[MA] ~ "for" ~ .(paste0(study, "_", j)))) +
  ylab("Density")+
  scale_linetype_discrete(name=NULL,labels=c("Controls","Cases"))+
  scale_fill_manual(name = NULL, values = c("#2171B5", "#D02020"), labels = c("Controls", "Cases")) +
  scale_color_manual(name = NULL, values = c("#2171B5", "#D02020"), labels = c("Controls", "Cases"))+
  coord_cartesian(xlim = c(-6, 6))
p
ggsave(paste0("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/Figure/dis/Figure S_PRSMA_",study,"_",j,".png"), plot = p, width = 8, height = 6, units = "in", dpi = 300)
dev.off()
