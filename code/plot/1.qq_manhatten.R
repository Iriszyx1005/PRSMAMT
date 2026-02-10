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
####QQ plot#####
lc_meta = fread("../data/01META_All.txt.gz")
#lc_meta = subset(lc_meta,Freq1>0.01 & Freq1<0.99)
dta_plot = data.frame(SNP = lc_meta$MarkerName,Chromosome =lc_meta$chr,Position = lc_meta$pos,trait1 = lc_meta$PvalueARE)

p_values=dta_plot$trait1
z = qnorm(p_values/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3) 


setwd("/home/sshen/Disk_m2/PRS/plot")  
setnames(dta_plot,"trait1","P value")
CMplot(dta_plot,plot.type="q",box=FALSE,file="tiff",file.name="",dpi=300,
       conf.int=TRUE,conf.int.col=NULL,threshold.col="#4b8bad",threshold.lty=2,
       file.output=TRUE,verbose=TRUE,width=6,height=6)

dev.off()

####manhatten#####
new1=read.xlsx("/home/sshen/Disk_m2/PRS_yxzhang/indepent_snp/new/new_region_ref_ng_catalog_1204.xlsx")
new1$ID=paste0(new1$Chr,":",new1$Position)
dta_plot$Cytoband <- new1$Cytoband[match(dta_plot$SNP, new1$ID)]

#a=dta_plot[1:10000,]
#dta_plot=dta_plot[!is.na(dta_plot$Cytoband), ]
#dta_plot=rbind(a,dta_plot)

#===================================================================
# 2. Calculate -log10(P)
#===================================================================
dta_plot$logp <- -log10(dta_plot$trait1)

#===================================================================
#3. Calculate the length and cumulative position of each chromosome
#===================================================================
chr_len <- dta_plot %>%
  group_by(Chromosome) %>%
  summarise(chr_len = max(Position))

chr_pos <- chr_len %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

# Add the chromosome offset to the original position
dta_plot <- dta_plot %>%
  left_join(chr_pos, by = "Chromosome") %>%
  arrange(Chromosome, Position) %>%
  mutate(poscum = Position + total)

#===================================================================
# 4.X-axis label position
#===================================================================
axis_set <- dta_plot %>%
  group_by(Chromosome) %>%
  summarize(center = (max(poscum) + min(poscum)) / 2)

dta_plot$is_highlight <- ifelse(!is.na(dta_plot$Cytoband), TRUE, FALSE)

tiff(paste0("manhatten_new_cytoband.tiff"),width = 12.5,height = 5,units = "in",res = 300)
p <- ggplot(dta_plot, aes(x = poscum, y = logp)) +
  geom_point(aes(color = as.factor(Chromosome)), size = 1.3, alpha = 0.8) +
  scale_color_manual(values = rep(c("#4b8bad", "#83bfbd"), length(unique(dta_plot$Chromosome)))) +
  ylim(0,115) +
  geom_point(
    data = subset(dta_plot, is_highlight),
    color = "#BE281B", size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_label_repel(
    data = subset(dta_plot, is_highlight),
    aes(label = Cytoband),
    size = 4,box.padding = 0.8,
    point.padding = 0.6, nudge_y = 10,  
    color = "#cf365a") +
  scale_x_continuous(label = axis_set$Chromosome, breaks = axis_set$center) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme( panel.background=element_blank(),
         legend.position="none",
         axis.line = element_line(colour = "black"),
         axis.title.x = element_text(size = 15),   
         axis.title.y = element_text(size = 15),  
         axis.text.x  = element_text(size = 13, color = "black", angle = 0, vjust = 1),
         axis.text.y  = element_text(size = 13, color = "black"))

print(p)
dev.off()
