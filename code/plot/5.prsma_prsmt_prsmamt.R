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

###prsma prsmt prsmamt####
afr=read.xlsx("../data/final_perform_afr.xlsx")
afr$race="African"

eas=read.xlsx("../data/final_perform_eas.xlsx")
eas$race="Asian"

eur=read.xlsx("../data/final_perform_eur.xlsx")
eur$race="European"


afr$Effect=round(afr$Effect,3)
afr$Cl=round(afr$Cl,3)
afr$Cu=round(afr$Cu,3)
eas$Effect=round(eas$Effect,3)
eas$Cl=round(eas$Cl,3)
eas$Cu=round(eas$Cu,3)
eur$Effect=round(eur$Effect,3)
eur$Cl=round(eur$Cl,3)
eur$Cu=round(eur$Cu,3)

x=rbind(eur,eas)
x=rbind(x,afr)

data=x[x$Metric=="AUC",]
data=data[order(data$Study),]
data$PRS=ifelse(data$PRS=="PRSMA","y0",ifelse(data$PRS=="PRSMT","y1","y2"))

data$PRS_race <- interaction(data$PRS,data$race) 

final_colors=c("#C1DEF2","#88B7DC","#3780BD",
               "#C0DAAE","#90B17D","#578144",
               "#FDEDA1", "#F7D76B", "#E1B844")
names(final_colors) <- unique(data$PRS_race) 

#final_colors=c("#EFEDF5","#DADAEB","#BCBDDC")
#names(final_colors) <-unique(data$y)
legend_labels <- c(
  "y0.Asian" = expression(paste(PRS["MA"], "_Asian")),
  "y1.Asian" = expression(paste(PRS["MT"], "_Asian")),
  "y2.Asian" = expression(paste(PRS["MAMT"], "_Asian")),
  "y0.African" = expression(paste(PRS["MA"], "_African")),
  "y1.African" = expression(paste(PRS["MT"], "_African")),
  "y2.African" = expression(paste(PRS["MAMT"], "_African")),
  "y0.European" = expression(paste(PRS["MA"], "_European")),
  "y1.European" = expression(paste(PRS["MT"], "_European")),
  "y2.European" = expression(paste(PRS["MAMT"], "_European"))
)


p=ggplot(data, aes(x=Study, y=Effect,color=PRS_race,fill=PRS_race)) + 
  geom_bar(stat = "identity", width = 0.7,position = position_dodge(width = 0.7)) + 
  scale_fill_manual(values=final_colors, labels = legend_labels) +
  scale_color_manual(values = final_colors, labels = legend_labels) +  
  geom_errorbar(aes(ymin = Cl, ymax = Cu),position = position_dodge(width = 0.7), 
                width = 0.2, color = "black", linewidth = 0.5) +
  # scale_x_discrete(labels = custom_labels)+
  #geom_hline(yintercept = 1,color = "#C9CACA", linetype = 2)+
  #labs(x="Study") +
  scale_y_continuous("AUC",expand =c(0,0),limits =c(0,0.75),n.breaks = 8) +
  theme(panel.background=element_blank(),
        axis.title = element_text(size=20,color="black"),
        #  plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        axis.text = element_text(size=20,color="black"),
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20),  
        axis.line = element_line(color="black" ,linewidth = 0.5),
        axis.ticks.length = unit(-1,"mm"),
        strip.text = element_text(size = 20),
        #     axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5, 5, 5, 30),panel.spacing = unit(1, "lines")
  ) + coord_cartesian(clip = "off") + scale_y_break(c(0,0.4))+facet_wrap(~race,ncol=3)
p
ggsave(paste0("../results/bar_y3_auc.png"),p,width = 20,height = 10,units = "in",dpi = 300)
dev.off()


data=x[x$Metric=="OR/SD",]
data=data[order(data$Study),]
data$PRS=ifelse(data$PRS=="PRSMA","y0",ifelse(data$PRS=="PRSMT","y1","y2"))

data$PRS_race <- interaction(data$PRS,data$race) 
final_colors=c("#C1DEF2","#88B7DC","#3780BD",
               "#C0DAAE","#90B17D","#578144",
               "#FDEDA1", "#F7D76B", "#E1B844")

names(final_colors) <- unique(data$PRS_race) 
legend_labels <- c(
  "y0.Asian" = expression(paste(PRS["MA"], "_Asian")),
  "y1.Asian" = expression(paste(PRS["MT"], "_Asian")),
  "y2.Asian" = expression(paste(PRS["MAMT"], "_Asian")),
  "y0.African" = expression(paste(PRS["MA"], "_African")),
  "y1.African" = expression(paste(PRS["MT"], "_African")),
  "y2.African" = expression(paste(PRS["MAMT"], "_African")),
  "y0.European" = expression(paste(PRS["MA"], "_European")),
  "y1.European" = expression(paste(PRS["MT"], "_European")),
  "y2.European" = expression(paste(PRS["MAMT"], "_European"))
)


p=ggplot(data, aes(x=Study, y=Effect,color=PRS_race,fill=PRS_race)) + 
  geom_bar(stat = "identity", width = 0.7,position = position_dodge(width = 0.7)) + 
  scale_fill_manual(values=final_colors, labels = legend_labels) +
  scale_color_manual(values = final_colors, labels = legend_labels) +  
  geom_errorbar(aes(ymin = Cl, ymax = Cu),position = position_dodge(width = 0.7), 
                width = 0.2, color = "black", linewidth = 0.5) +
  # scale_x_discrete(labels = custom_labels)+
  geom_hline(yintercept = 1, color = "#F46D43", linetype = 2)+
  # labs(x="Study") +
  scale_y_continuous("OR",expand =c(0,0),limits =c(0,3),n.breaks = 3) +
  theme(panel.background=element_blank(),
        axis.title = element_text(size=20,color="black"),
        # plot.title = element_text(size=20,face = "bold.italic",hjust=0.5),
        axis.text = element_text(size=20,color="black"),
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20),  
        axis.line = element_line(color="black" ,linewidth = 0.5),
        axis.ticks.length = unit(-1,"mm"),
        strip.text = element_text(size = 20),
        #     axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5, 5, 5, 30),panel.spacing = unit(1, "lines")
  )+coord_cartesian(clip = "off") +facet_wrap(~race,ncol=3)
p
ggsave(paste0("../results/bar_y3_or.png"),width = 20,height = 10,units = "in",dpi = 300)
dev.off()
