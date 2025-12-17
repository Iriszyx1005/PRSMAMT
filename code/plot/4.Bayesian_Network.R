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
library(Rgraphviz)
library(bnlearn)

load('BN_eur.Rdata')
load('BN_eas.RData')
load('BN_afr.RData')
ls()
head(dag)
# rm(list = ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install(c("graph", "Rgraphviz"))
# install.packages("bnlearn")

gaussian.test = data.frame(nodes,LC=y)
# learn the network structure.
blacklist = tiers2blacklist(list(colnames(nodes),"LC"))
blacklist = rbind(blacklist,c("COPD","Smoking_initiation"))
dag = mmhc(gaussian.test,blacklist = blacklist)
# estimate the parameters of the Bayesian network
fitted = bn.fit(dag, gaussian.test, method = "mle-g")
gplot = graphviz.plot(fitted, shape = "ellipse", layout = "fdp")
plot(gplot)
save(gplot,gaussian.test,dag,fitted,file = "/home/sshen/Disk_m2/PRS/plot/BN_eur.RData")


#####
load("/Users/zhangyixin/Desktop/PRS/Bayesian_Network/BN_eur.RData")
library(igraph)
library(ggraph)
library(dplyr)
library(grid)
plot(ig)
# 1. Convert from bnlearn object to igraph object
# Here we use dag, directly building igraph using arcs information
edges <- arcs(fitted)
edges_df <- data.frame(edges)
edges_df <- edges_df %>%
  mutate(
    from = recode(from,
                  fev1 = "Fev1",
                  fvc = "Fvc",
                  fev1fvc = "Fev1fvc",
                  F17 = "Tobacco_disorder",
                  F17_2 = "Nicotine_dependence",
                  smoking_initiation = "Smoking_initiation",
                  age_of_smoking = "Age_of_smoking",
                  smoking_cessation = "Smoking_cessation",
                  emphysema = "Emphysema",
                  y = "LC",
                  .default = from),
    to = recode(to,
                fev1 = "Fev1",
                fvc = "Fvc",
                fev1fvc = "Fev1fvc",
                F17 = "Tobacco_disorder",
                F17_2 = "Nicotine_dependence",
                smoking_initiation = "Smoking_initiation",
                age_of_smoking = "Age_of_smoking",
                smoking_cessation = "Smoking_cessation",
                emphysema = "Emphysema",
                y = "LC",
                .default = to)
  )

ig <- graph_from_data_frame(edges_df, directed = TRUE)

# ig <- graph_from_data_frame(edges, directed = TRUE)

#2.Node Attributes: Grouping and Color
node_groups <- list(
  "LC" = "LC",
  "Smoking" = c("Age_of_smoking", "Smoking_initiation", "Smoking_cessation",
                "Cigarettes_per_day", "Nicotine_dependence", "Tobacco_disorder"),
  "Lung_Function" = c("Fev1", "Fvc", "Fev1fvc"),
  "CRP" = "CRP",
  "Respiratory_Disease" = c("Emphysema", "COPD", "Asthma", "IPF")
)

colors <- c(
  "LC" = "#F6B352",
  "Smoking" = "#d3a9c5",
  "Lung_Function" = "#7babdd",
  "CRP" = "#BDD088",
  "Respiratory_Disease" = "#8cd3d6"
)

# Set color attribute for igraph nodes
V(ig)$group <- sapply(V(ig)$name, function(n) {
  grp <- NA
  for (g in names(node_groups)) {
    if (n %in% node_groups[[g]]) {
      grp <- g
      break
    }
  }
  if (is.na(grp)) grp <- "Other"
  return(grp)
})

V(ig)$color <- colors[V(ig)$group]

# 3. Set the layout so that LC is in the center, with others surrounding it.
#Place LC at the center (0,0), and arrange the other nodes in a circular pattern around it.
center_node <- "LC"
others <- setdiff(V(ig)$name, center_node)
n_others <- length(others)
radius <- 1

#group_order already has an order, first sort by group, then others
group_order <- names(node_groups)
node_info <- data.frame(name = V(ig)$name, group = V(ig)$group, stringsAsFactors = FALSE)

# Sort by group order and name
ordered_others <- node_info %>%
  filter(name != center_node) %>%
  arrange(factor(group, levels = group_order), name) %>%
  pull(name)

# Distribute angles evenly on the circumference for points excluding LC
angle <- seq(0, 2 * pi, length.out = n_others + 1)[-(n_others + 1)]

layout_df <- data.frame(
  name = c(center_node, ordered_others),
  angle = c(0, angle),
  stringsAsFactors = FALSE
)

# Calculating Coordinates and Text Alignment Properties
layout_df <- layout_df %>%
  mutate(
    x = ifelse(name == center_node, 0, radius * cos(angle)),
    y = ifelse(name == center_node, 0, radius * sin(angle)),
    size = 13,
    hjust = ifelse(name == center_node, 0.5, ifelse(x >= 0, 0, 1)),
    nudge_x = ifelse(name == center_node, 0, ifelse(x >= 0, 0.1, -0.1))
  )

# Assign values back to igraph node attributes
V(ig)$x <- layout_df$x[match(V(ig)$name, layout_df$name)]
V(ig)$y <- layout_df$y[match(V(ig)$name, layout_df$name)]
V(ig)$hjust <- layout_df$hjust[match(V(ig)$name, layout_df$name)]
V(ig)$nudge_x <- layout_df$nudge_x[match(V(ig)$name, layout_df$name)]
V(ig)$size <- layout_df$size[match(V(ig)$name, layout_df$name)]
V(ig)$size[V(ig)$name == "LC"] <- 16.5  
V(ig)$fontface <- ifelse(V(ig)$name == "LC", "bold", "plain")


E(ig)$width <- 0.3
E(ig)$width[head_of(ig, E(ig))$name == center_node] <- 0.5
E(ig)$color <- "#777777"
E(ig)$color[head_of(ig, E(ig))$name == center_node] <- "#BE281B"


p <- ggraph(ig, layout = "manual", x = V(ig)$x, y = V(ig)$y) +
  geom_edge_arc(
    aes(edge_width = width, color = color),
    arrow = arrow(length = unit(1.26, 'mm'), type = 'closed',angle = 22),
    end_cap = circle(6.8, 'mm'), start_cap = circle(6, 'mm'),
    curvature = 0.1
  ) +
  geom_node_point(
    aes(fill = group, size = size),
    shape = 21,
    color = "white",  
    stroke = 0.5,
    show.legend = FALSE
  ) +
  geom_node_text(
    aes(label = gsub("_", " ", name), hjust = hjust,fontface = fontface),
    nudge_x = V(ig)$nudge_x,
    size = 7,
    show.legend = FALSE
  ) +
  coord_cartesian(
    xlim = c(min(V(ig)$x) - 0.5, max(V(ig)$x) + 0.5)
  ) +
  scale_fill_manual(values = colors) +   
  scale_edge_width_identity() +
  scale_edge_color_identity() +
  scale_size_identity(guide = "none") +
  theme_void() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
ggsave("/Users/zhangyixin/Desktop/PRS/z_script/Figure/network_eur1.tiff", plot = p, width = 12.5, height =7.5, dpi = 600)

