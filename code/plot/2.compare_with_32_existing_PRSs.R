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

###compare with 32 existing PRS models#######
#####OR#####
x=fread("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/auc_or_32+1+1.txt")
x <- x[x$prs != "PRS", ]
x$prs[x$prs=="PRS_prscsx"]="PRSMA"

x$auc_cilower=ifelse(x$AUC<0.5,1-x$auc_ciupper,x$auc_cilower)
x$auc_ciupper=ifelse(x$AUC<0.5,1-x$auc_cilower,x$auc_ciupper)
x$AUC=ifelse(x$AUC<0.5,1-x$AUC,x$AUC)
x$or_cilower=ifelse(x$OR<1,1/x$or_ciupper,x$or_cilower)
x$or_ciupper=ifelse(x$OR<1,1/x$or_cilower,x$or_ciupper)
x$OR=ifelse(x$OR<1,1/x$OR,x$OR)


all_colors <- brewer.pal.info["Blues", "maxcolors"]
selected_colors <- brewer.pal(all_colors, "Blues") 
subset_colors <- selected_colors[c(3:7)] 
colors_32 <- colorRampPalette(subset_colors)(32)

plot_OR <- function(data, ancestry,label, show_y = TRUE) {
  highlight_data <- data[data$OR > 4, ]
  p <- ggplot(data, aes(x=prs, y=OR,fill = prs)) + 
    geom_bar(stat = "identity", width = 0.7) + 
    scale_fill_manual(values=final_colors) +
    geom_errorbar(aes(ymin = or_cilower, ymax = or_ciupper), 
                  width = 0.2, color = "black", linewidth = 0.5) +
    scale_x_discrete(labels = custom_labels)+
    geom_hline(yintercept = 1, color = "#F46D43", linetype = 2)+
    labs(title=ancestry,x="PRS") +
    scale_y_continuous("OR",expand =c(0,0),n.breaks = 4) +
    theme(panel.background=element_blank(),legend.position = "none",
          axis.title = element_text(size=15,color="black"),
          plot.title = element_text(size=15,face = "bold.italic",hjust=0.5),
          axis.text.x = element_text(size=8,color="black",angle = 45, hjust = 1),
          axis.text.y = if (show_y) element_text(size=15,color="black") else element_blank(),
          axis.line = element_line(color="black" ,linewidth = 0.5),
          axis.ticks.length = unit(-1,"mm"),
          axis.title.y = if (show_y) element_text() else element_blank(),
          plot.margin = margin(5, 5, 5, 30) 
    )+coord_cartesian(clip = "off",ylim = c(0, 4)) + 
    geom_text(data = highlight_data, 
              aes(x = prs, y = 4, label = paste0("↑ ", round(OR, 2))), 
              vjust = -0.3, size = 3.5) +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = 3, vjust = -1, 
             fontface = "bold", size = 5)
  
  return(p)
}


plots_OR <- list()
#plots_auc <- list()
labels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")

i <- 1
for (pop in unique(x$ancestry)) {
  data=x[x$ancestry==pop,]
  data$prs[data$prs=="PRS"]="PRSMA"
  data <- data[order(data$prs != "PRSMA", data$prs), ]
  
  data$prs <- factor(data$prs, levels = c("PRSMA","PGS000070","PGS000078","PGS000156","PGS000388",
                                          "PGS000389","PGS000390","PGS000391","PGS000392","PGS000393",
                                          "PGS000394","PGS000395","PGS000396","PGS000397","PGS000721",
                                          "PGS000740","PGS000789","PGS000880","PGS002270","PGS002808",
                                          "PGS003391","PGS003392","PGS003393","PGS004164","PGS004165",
                                          "PGS004246","PGS004325","PGS004442","PGS004512","PGS004691",
                                          "PGS004860","PGS004884","PGS004955"))
  
  
  custom_labels= c(expression("PRS"["MA"]),"PGS000070","PGS000078","PGS000156","PGS000388",
                   "PGS000389","PGS000390","PGS000391","PGS000392","PGS000393",
                   "PGS000394","PGS000395","PGS000396","PGS000397","PGS000721",
                   "PGS000740","PGS000789","PGS000880","PGS002270","PGS002808",
                   "PGS003391","PGS003392","PGS003393","PGS004164","PGS004165",
                   "PGS004246","PGS004325","PGS004442","PGS004512","PGS004691",
                   "PGS004860","PGS004884","PGS004955")
  
  auc_colors=rep("#559FCD",33)
  names(auc_colors) <- data$prs
  auc_colors["PRSMA"] <- "#BE281B"
  
  final_colors=c("#BE281B",colors_32)
  names(final_colors) <- data$prs
  if (pop %in% c("OncoArray_European")) {
    final_colors["PGS004860"] <-  "#D9D9D9"
    final_colors["PGS004691"] <-  "#D9D9D9"
    final_colors["PGS004246"] <- "#D9D9D9"
    final_colors["PGS002270"] <-  "#D9D9D9"
    final_colors["PGS000880"] <-  "#D9D9D9"
    final_colors["PGS000740"] <-  "#D9D9D9"
    final_colors["PGS000078"] <-  "#D9D9D9"
    auc_colors["PGS004860"] <-  "#D9D9D9"
    auc_colors["PGS004691"] <-  "#D9D9D9"
    auc_colors["PGS004246"] <-  "#D9D9D9"
    auc_colors["PGS002270"] <-  "#D9D9D9"
    auc_colors["PGS000880"] <-  "#D9D9D9"
    auc_colors["PGS000740"] <-  "#D9D9D9"
    auc_colors["PGS000078"] <-  "#D9D9D9"
  }
  if (pop %in% c("TRICL_European")) {
    final_colors["PGS004860"] <- "#C9CACA"
    final_colors["PGS004691"] <- "#C9CACA"
    final_colors["PGS004246"] <- "#C9CACA"
    final_colors["PGS003393"] <- "#C9CACA"
    final_colors["PGS003392"] <- "#C9CACA"
    final_colors["PGS003391"] <- "#C9CACA"
    final_colors["PGS002270"] <- "#C9CACA"
    final_colors["PGS000880"] <- "#C9CACA"
    final_colors["PGS000740"] <- "#C9CACA"
    final_colors["PGS000078"] <- "#C9CACA"
    auc_colors["PGS004860"] <- "#C9CACA"
    auc_colors["PGS004691"] <- "#C9CACA"
    auc_colors["PGS004246"] <- "#C9CACA"
    auc_colors["PGS003393"] <- "#C9CACA"
    auc_colors["PGS003392"] <- "#C9CACA"
    auc_colors["PGS003391"] <- "#C9CACA"
    auc_colors["PGS002270"] <- "#C9CACA"
    auc_colors["PGS000880"] <- "#C9CACA"
    auc_colors["PGS000740"] <- "#C9CACA"
    auc_colors["PGS000078"] <- "#C9CACA"
  }
  if (pop %in% c("OncoArray_Asian")) {
    final_colors["PGS004164"] <- "#C9CACA"
    final_colors["PGS004165"] <- "#C9CACA"
    auc_colors["PGS004164"] <- "#C9CACA"
    auc_colors["PGS004165"] <- "#C9CACA"
  }
  if (pop %in% c("FLCCA_Asian")) {
    final_colors["PGS004164"] <- "#C9CACA"
    final_colors["PGS004165"] <- "#C9CACA"
    final_colors["PGS000078"] <- "#C9CACA"
    auc_colors["PGS004164"] <- "#C9CACA"
    auc_colors["PGS004165"] <- "#C9CACA"
    auc_colors["PGS000078"] <- "#C9CACA"
  }
  if (pop %in% c("NCI_African")) {
    final_colors["PGS000078"] <- "#C9CACA"
    auc_colors["PGS000078"] <- "#C9CACA"
  }
  
  plots_OR[[i]] <- plot_OR(data, pop,labels[i], show_y =TRUE) 
  i <- i + 1
  
}


final_plot_OR <- plot_grid(plotlist = plots_OR, ncol = 3, align = "hv")

p=plot_grid(final_plot_OR, ncol = 1, rel_heights = c(1, 0.05))
p
ggsave("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/Figure/combine_or.png", p, width =15, height = 10, dpi = 300)

#####AUC#####
x$group=sapply(x$ancestry,function(x) strsplit(x,"_")[[1]][1])
x[which(x$group  %in% c("PLCO","FLCCA","NCI")),]$group="PLCO/FLCCA/NCI"

plot_AUC <- function(data, group,label, show_y = TRUE) {
  highlight_data <- data[data$AUC > 0.8, ]
  p  <- ggplot(data, aes(y=prs, x=AUC,color = color_group, fill = color_group)) +
    geom_point(size=5,pch = 19,position = position_dodge(width = 0.6)) + 
    geom_errorbar(aes(xmin = auc_cilower, xmax = auc_ciupper), width = 0.5, linewidth = 0.7,
                  position = position_dodge(width = 0.6)) +  
    # scale_x_discrete(labels = custom_labels)+
    scale_fill_manual(values = auc_colors ) + 
    scale_color_manual(values = auc_colors) + 
    geom_vline(xintercept = 0.5, color = "#C9CACA", linetype = 2)+
    labs(title=group,fill = NULL, color = NULL) +
    scale_x_continuous("AUC",expand =c(0,0),n.breaks = 5) +
    labs(y = "PRSs")+
    theme(panel.background=element_blank(),
          axis.title = element_text(size=25,color="black"),
          plot.title = element_text(size=25,face="bold.italic",hjust=0.5),
          axis.text = element_text(size=20,color="black"),
          axis.line = element_line(color="black" ,linewidth = 0.5),
          axis.ticks.length = unit(-1,"mm"),
          #  axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          axis.text.y = if (show_y) element_text() else element_blank(),
          axis.title.y = if (show_y) element_text() else element_blank(),panel.spacing = unit(0, "lines"),
          plot.margin = margin(20,32,6,15))+
    coord_cartesian(clip = "off",xlim = c(0.3,0.8)) + 
    geom_text(data = highlight_data, 
              aes(x = 0.8, y = prs, label = paste0(round(AUC, 3))), 
              vjust = -0.3, size = 5) +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = 3, vjust = -1.5,
             fontface = "bold", size = 9) 
  p
  
  return(p)
}

plots_auc <- list()
labels <- c("A", "B", "C","D")
order_group=c("OncoArray","TRICL","AoU","PLCO/FLCCA/NCI")
i <- 1
for (pop in order_group) {
  data=x[x$group==pop,]
  
  data$color_group <- data$race
  
  if (pop %in% c("OncoArray")) {
    data[data$prs=="PGS004860" & data$race=="European",]$color_group <-  "Overlap"
    data[data$prs=="PGS004691" & data$race=="European",]$color_group <-  "Overlap"
    data[data$prs=="PGS004246" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS002270" & data$race=="European",]$color_group <-  "Overlap"
    data[data$prs=="PGS000880" & data$race=="European",]$color_group <-  "Overlap"
    data[data$prs=="PGS000740" & data$race=="European",]$color_group <-  "Overlap"
    data[data$prs=="PGS000078" & data$race=="European",]$color_group <-  "Overlap"
    data[data$prs=="PGS004164" & data$race=="Asian",]$color_group <- "Overlap"
    data[data$prs=="PGS004165" & data$race=="Asian",]$color_group <- "Overlap"
  }
  if (pop %in% c("TRICL")) {
    data[data$prs=="PGS004860" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS004691" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS004246" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS003393" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS003392" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS003391" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS002270" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS000880" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS000740" & data$race=="European",]$color_group <- "Overlap"
    data[data$prs=="PGS000078" & data$race=="European",]$color_group <- "Overlap"
  }
  
  if (pop %in% c("PLCO/FLCCA/NCI")) {
    data[data$prs=="PGS004164" & data$race=="Asian",]$color_group <- "Overlap"
    data[data$prs=="PGS004165" & data$race=="Asian",]$color_group <- "Overlap"
    data[data$prs=="PGS000078" & data$race=="Asian",]$color_group <- "Overlap"
    data[data$prs=="PGS000078" & data$race=="African",]$color_group <- "Overlap"
  }
  
  #auc_colors= c("#559FCD","#7FBC41","#F4C700")
  #names(auc_colors) <- unique(data$race)
  auc_colors <- c("grey", "#559FCD", "#7FBC41", "#F4C700")
  names(auc_colors) <- c("Overlap","European","Asian","African" )
  
  plots_auc[[i]] <- plot_AUC(data, pop,labels[i], show_y = (i %% 4 == 1)) 
  i <- i + 1
}
final_plot_auc <- plot_grid(plotlist = plots_auc, ncol = 4, align = "hv")


com_p=plot_grid(final_plot_auc, ncol = 1, rel_heights = c(0.5, 0.01),rel_widths = c(2, 2, 2,2))
com_p

ggsave("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/Figure/combine_auc.pdf", com_p, width =23, height = 15, dpi = 300)
dev.off()

#####rank####
x=fread("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/result/auc_or_32+1+1.txt")
x$auc_cilower=ifelse(x$AUC<0.5,1-x$auc_ciupper,x$auc_cilower)
x$auc_ciupper=ifelse(x$AUC<0.5,1-x$auc_cilower,x$auc_ciupper)
x$AUC=ifelse(x$AUC<0.5,1-x$AUC,x$AUC)
x$or_cilower=ifelse(x$OR<1,1/x$or_ciupper,x$or_cilower)
x$or_ciupper=ifelse(x$OR<1,1/x$or_cilower,x$or_ciupper)
x$OR=ifelse(x$OR<1,1/x$OR,x$OR)
unique(x$ancestry)
x <- x[x$prs != "PRS", ]
x$prs[x$prs=="PRS_prscsx"]="PRSMA"

x=x[-which(x$prs=="PGS004860" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS004860" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS004691" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS004691" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS004246" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS004246" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS004164" & x$ancestry=="Onco_eas"),]
x=x[-which(x$prs=="PGS004164" & x$ancestry=="FLCCA_eas"),]

x=x[-which(x$prs=="PGS004165" & x$ancestry=="Onco_eas"),]
x=x[-which(x$prs=="PGS004165" & x$ancestry=="FLCCA_eas"),]

x=x[-which(x$prs=="PGS003393" & x$ancestry=="TRICL_eur"),]
x=x[-which(x$prs=="PGS003392" & x$ancestry=="TRICL_eur"),]
x=x[-which(x$prs=="PGS003391" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS002270" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS002270" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS000880" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS000880" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS000740" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS000740" & x$ancestry=="TRICL_eur"),]

x=x[-which(x$prs=="PGS000078" & x$ancestry=="Onco_eur"),]
x=x[-which(x$prs=="PGS000078" & x$ancestry=="TRICL_eur"),]
x=x[-which(x$prs=="PGS000078" & x$ancestry=="FLCCA_eas"),]
x=x[-which(x$prs=="PGS000078" & x$ancestry=="NCI_afr"),]


plot_rank <- function(data, group,label, show_y = TRUE) {
  
  p <- ggplot(data, aes(x=prs, y=rank_auc,fill = prs)) + 
    geom_bar(stat = "identity", width = 0.7) + 
    scale_fill_manual(values=colors_33) +
    scale_x_discrete(labels = custom_labels)+
    labs(x="PRSs",title=paste0(group)) +
    scale_y_continuous("Average rank",expand =c(0,0),limits = c(0,33),n.breaks = 6) +
    theme(panel.background=element_blank(),
          axis.title = element_text(size=15,color="black"),
          plot.title = element_text(size=15,face = "bold.italic",hjust=0.5),
          axis.text.x = element_text(size=12,color="black",angle = 45, hjust = 1),
          axis.line = element_line(color="black" ,linewidth = 0.5),
          axis.ticks.length = unit(-1,"mm"),
          legend.position = "none",
          axis.text.y = if (show_y) element_text(size=15,color="black") else element_blank(),
          axis.title.y = if (show_y) element_text() else element_blank(),panel.spacing = unit(0, "lines"),
          plot.margin = margin(5,5,5,30))+coord_cartesian(clip = "off") +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = 3, vjust = 0,
             fontface = "bold", size = 7) 
  return(p)
}

plots_rank <- list()
labels <- c("A", "B", "C","D")
group=c("Overall","European","Asian","African")
i <- 1
for (pop in group) {
  if (pop =="Overall"){
    df=x
  }else{
    df=x[x$race==pop,]
  }
  
  #!!OR/AUC
  df <- df %>%
    group_by(ancestry) %>%
    mutate(rank_auc = rank(-AUC)) 
  df_auc_avgrank <- aggregate(rank_auc ~ prs, data = df, FUN = mean)
  a=df_auc_avgrank$rank_auc[df_auc_avgrank$prs=="PRSMA"]
  df_auc_avgrank$prs[df_auc_avgrank$rank_auc<a]
  df_auc_avgrank$rank_auc=round(df_auc_avgrank$rank_auc,0)
  
  if (pop =="Overall" | pop =="Asian"){
    df_auc_avgrank <- df_auc_avgrank[order(df_auc_avgrank$rank_auc), ]
    
    custom_labels= c(expression("PRS"["MA"]),df_auc_avgrank$prs[-1])
    
  }else{
    df_auc_avgrank <- df_auc_avgrank[order(df_auc_avgrank$prs != "PRSMA", df_auc_avgrank$rank_auc), ]
    custom_labels= c(expression("PRS"["MA"]),df_auc_avgrank$prs[-1]) #black or/auc
    
  }
  
  df_auc_avgrank$prs <- factor(df_auc_avgrank$prs,levels = c(df_auc_avgrank$prs))
  
  all_colors <- brewer.pal.info["Blues", "maxcolors"]
  selected_colors <- brewer.pal(all_colors, "Blues") 
  subset_colors <- selected_colors[c(3:7)] 
  colors_33 <- colorRampPalette(subset_colors)(33)
  #colors_33 <- colorRampPalette(subset_colors)(34)
  names(colors_33) <- df_auc_avgrank$prs
  colors_33["PRSMA"] <- "#BE281B"
  
  plots_rank[[i]] = plot_rank(df_auc_avgrank, pop,labels[i],  show_y = TRUE) 
  
  i <- i + 1
}


final_plot_rank <- plot_grid(plotlist = plots_rank, ncol = 2, align = "hv")


com_p=plot_grid(final_plot_rank, ncol = 1, rel_heights = c(0.5, 0.01),rel_widths = c(2, 2, 2,2))
com_p
ggsave("/Users/zhangyixin/Desktop/03Nature_Communications/z_resubmit/Figure/rank_or.pdf", com_p,  width =20, height = 10, dpi = 300)


