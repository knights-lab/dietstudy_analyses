
require(tidyverse)
require(reshape2)
require(ape)
require(vegan)
require(gridExtra)
require(cowplot)
 
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/")

set.seed(7)



source("lib/results_scripts/4a_Microbiome and food pair longitudinally/1_Figure3/load_data_lists.R")

### TODO: reduce to just .txt files in procrustes data


indiv_pro1 <- as.list(NULL)
pval.list.1 <- as.list(NULL)

for (n in names(mb.pcoa.1day)) {
  
  pcoa_f <- diet.pcoa.1day[[n]]
  pcoa_t <- mb.pcoa.1day[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_f, pcoa_t)
  pro_test <- protest(pcoa_f, pcoa_t, perm = 9999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Microbiome Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  
  indiv_pro1[[n]] <- ggplot(plot) +
    geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, fill = type), pch=21) +
    scale_fill_manual(values = c("#5dd047", "#E71D36")) +
    theme_bw() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "black", alpha = 0.7) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(.5, .5, .5, .5, "pt"),
          plot.title = element_text(size = 6, hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size=9),
          legend.position = 'none',
          axis.text = element_text(size=4),
          axis.title = element_text(size=4),
          aspect.ratio = 1) +
    #annotate("text", x = -0.17, y = -0.22, label = paste0("p-value=",pval), size = 3) +
    ggtitle(n) +
    xlab(paste0("PC 1 [",percent_var[1],"%]")) +
    ylab(paste0("PC 2 [",percent_var[2],"%]")) 
  
  pval.list.1[[n]] <- pval
  
}




indiv_pro2 <- as.list(NULL)
pval.list.2 <- as.list(NULL)

for (n in names(mb.pcoa.2day)) {
  
  pcoa_f <- diet.pcoa.2day[[n]]
  pcoa_t <- mb.pcoa.2day[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_f, pcoa_t)
  pro_test <- protest(pcoa_f, pcoa_t, perm = 9999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Microbiome Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  
  indiv_pro2[[n]] <- ggplot(plot) +
    geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, fill = type), pch=21) +
    scale_fill_manual(values = c("#5dd047", "#E71D36")) +
    theme_bw() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "black", alpha = 0.7) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(.5, .5, .5, .5, "pt"),
          plot.title = element_text(size = 6, hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size=9),
          legend.position = 'none',
          axis.text = element_text(size=4),
          axis.title = element_text(size=4),
          aspect.ratio = 1) +
    #annotate("text", x = -0.17, y = -0.22, label = paste0("p-value=",pval), size = 3) +
    ggtitle(n) +
    xlab(paste0("PC 1 [",percent_var[1],"%]")) +
    ylab(paste0("PC 2 [",percent_var[2],"%]")) 
  
  pval.list.2[[n]] <- pval
  
}



indiv_pro3 <- as.list(NULL)
pval.list.3 <- as.list(NULL)

for (n in names(mb.pcoa.3day)) {
  
  pcoa_f <- diet.pcoa.3day[[n]]
  pcoa_t <- mb.pcoa.3day[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_f, pcoa_t)
  pro_test <- protest(pcoa_f, pcoa_t, perm = 9999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Microbiome Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  
  indiv_pro3[[n]] <- ggplot(plot) +
    geom_point(size = 3, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
    scale_color_manual(values = c("#5a2071", "#5f86b7")) +
    theme_bw() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(.01, .01, .01, .01, "in"),
          plot.title = element_text(size = 9, hjust = 0.5),
          legend.position = 'none',
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1) +
   # ggtitle(gsub("MCTs", "", n)) +
  NULL
  
  pval.list.3[[n]] <- pval
  
}





#### keep running the procrustes with more days of diet ###
indiv_pro4 <- as.list(NULL)
pval.list.4 <- as.list(NULL)

for (n in names(mb.pcoa.4day)) {
  
  pcoa_f <- diet.pcoa.4day[[n]]
  pcoa_t <- mb.pcoa.4day[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_f, pcoa_t)
  pro_test <- protest(pcoa_f, pcoa_t, perm = 9999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Microbiome Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  
  indiv_pro4[[n]] <- ggplot(plot) +
    geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, fill = type), pch=21) +
    scale_fill_manual(values = c("#5dd047", "#E71D36")) +
    theme_bw() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "black", alpha = 0.7) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(.5, .5, .5, .5, "pt"),
          plot.title = element_text(size = 6, hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size=99),
          legend.position = 'none',
          axis.text = element_text(size=4),
          axis.title = element_text(size=4),
          aspect.ratio = 1) +
    #annotate("text", x = -0.17, y = -0.22, label = paste0("p-value=",pval), size = 3) +
    ggtitle(n) +
    xlab(paste0("PC 1 [",percent_var[1],"%]")) +
    ylab(paste0("PC 2 [",percent_var[2],"%]")) 
  
  pval.list.4[[n]] <- pval
  
}




indiv_pro5 <- as.list(NULL)
pval.list.5 <- as.list(NULL)

for (n in names(mb.pcoa.5day)) {
  
  pcoa_f <- diet.pcoa.5day[[n]]
  pcoa_t <- mb.pcoa.5day[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_f, pcoa_t)
  pro_test <- protest(pcoa_f, pcoa_t, perm = 9999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Microbiome Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  
  indiv_pro5[[n]] <- ggplot(plot) +
    geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, fill = type), pch=21) +
    scale_fill_manual(values = c("#5dd047", "#E71D36")) +
    theme_bw() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "black", alpha = 0.7) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(.5, .5, .5, .5, "pt"),
          plot.title = element_text(size = 6, hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size = unit(0.2, "in"),
          legend.text = element_text(size=9),
          legend.position = 'none',
          axis.text = element_text(size=4),
          axis.title = element_text(size=4),
          aspect.ratio = 1) +
    #annotate("text", x = -0.17, y = -0.22, label = paste0("p-value=",pval), size = 3) +
    ggtitle(n) +
    xlab(paste0("PC 1 [",percent_var[1],"%]")) +
    ylab(paste0("PC 2 [",percent_var[2],"%]")) 
  
  pval.list.5[[n]] <- pval
  
}





indiv_prodecay <- as.list(NULL)
pval.list.decay <- as.list(NULL)

for (n in names(mb.pcoa.decay)) {
  
  pcoa_f <- diet.pcoa.decay[[n]]
  pcoa_t <- mb.pcoa.decay[[n]]
  
  # procrustes
  pro <- procrustes(pcoa_f, pcoa_t)
  pro_test <- protest(pcoa_f, pcoa_t, perm = 9999) # may want to skip this for plotting?
  
  eigen <- sqrt(pro$svd$d)
  percent_var <- signif(eigen/sum(eigen), 4)*100
  
  beta_pro <- data.frame(pro$X)
  trans_pro <- data.frame(pro$Yrot)
  beta_pro$UserName <- rownames(beta_pro)
  beta_pro$type <- "Food Distance (Unweighted Unifrac)"
  trans_pro$UserName <- rownames(trans_pro)
  trans_pro$type <- "Microbiome Distance (Aitchison's)"
  
  colnames(trans_pro) <- colnames(beta_pro)
  
  pval <- pro_test$signif
  
  plot <- rbind(beta_pro, trans_pro)
  plot$StudyDay <- rep(1:dim(beta_pro)[1],2)
  
  indiv_prodecay[[n]] <- ggplot(plot) +
    geom_point(size = 2, alpha=0.75, aes(y = Axis.1, x = StudyDay, color = type)) +
    scale_color_manual(values = c("#5a2071", "#5f86b7")) +
    theme_bw() +
    geom_line(aes(y= Axis.1, x=StudyDay, group=UserName), col = "darkgrey", alpha = 0.6) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(.01, .01, .01, .01, "in"),
          plot.title = element_text(size = 9, hjust = 0.5),
          legend.position = 'none',
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1) +
    # ggtitle(gsub("MCTs", "", n)) +
    NULL
  
  pval.list.decay[[n]] <- pval
  
}


# sort for plotting A #
new_order <- names(sort(unlist(pval.list.decay[pval.list.decay < 0.05])))

order_indiv_prodecay <- indiv_prodecay[new_order]

new_order

names(order_indiv_prodecay)



### PLOT 3A  fig.height=2.7, fig.width=7 ####

grid.arrange(grobs = order_indiv_prodecay, nrow=3)

fig3A <- arrangeGrob(grobs = order_indiv_prodecay, nrow=3)

ggsave("../../../output/Figure3/Figure3A.pdf", fig3A, width = 7, height = 2.7)




### look at pvalues ###
plot <- as.data.frame(unlist(pval.list.1))
plot$pval.list.2 <- unlist(pval.list.2)
plot$pval.list.3 <- unlist(pval.list.3)
plot$pval.list.4 <- unlist(pval.list.4)
plot$pval.list.5 <- unlist(pval.list.5)
plot$pval.list.decay <- unlist(pval.list.decay)

colnames(plot) <- c("day1Monte Carlo p-value", "day2Monte Carlo p-value", "day3Monte Carlo p-value", "day4Monte Carlo p-value", "day5Monte Carlo p-value", "decayMonte Carlo p-value")
plot$UserName <- rownames(plot)




plot <- plot[order(plot$`day2Monte Carlo p-value`),]     # reorder by variable of choice
plot$UserName <- factor(plot$UserName, levels = plot$UserName)

plot_new <- plot %>% select(UserName, 
                        `day1Monte Carlo p-value`,
                        `day2Monte Carlo p-value`)

plot_melt <- melt(plot_new)
plot_melt <- na.omit(plot_melt)

# Make scatterplot of best v. 1 day, ordered by best
ggplot(plot_melt, aes(x = UserName, y = log10(value), group = variable)) + 
  geom_point(aes(color = variable), size = 3) +
  geom_point(shape = 1, size =3, color = "black") +
  scale_color_manual(values = c("#64baaa", "#5f86b7"), labels=c("1 day","2 day")) +
  geom_hline(yintercept = log10(0.05), linetype = 2, color = "darkgrey") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4),
        legend.position = c(0.88,0.13),
        legend.title = element_blank())+
  xlab(label = "Subject") +
  ylab(label = "log(Monte Carlo p-value)") +
  annotate("text", label = "p = 0.05", size = 3, x = 3, y = log10(0.04)) +
  scale_x_discrete(labels = gsub("MCTs", "", plot_melt$UserName))




plot2 <- plot[order(plot$`day3Monte Carlo p-value`),]     # reorder by variable of choice
plot2$UserName <- factor(plot2$UserName, levels = plot2$UserName)

plot2_new <- plot2 %>% select(UserName, 
                        `day1Monte Carlo p-value`,
                        `day3Monte Carlo p-value`)

plot2_melt <- melt(plot2_new)
plot2_melt <- na.omit(plot2_melt)

# Make scatterplot of best v. 1 day, ordered by best
ggplot(plot2_melt, aes(x = UserName, y = log10(value), group = variable)) + 
  geom_point(aes(color = variable), size = 3) +
  geom_point(shape = 1, size =3, color = "black") +
  scale_color_manual(values = c("#64baaa", "#7b3294"), labels=c("1 day","3 day")) +
  #scale_color_manual(values = brewer.pal(5, "Dark2")) +
  geom_hline(yintercept = log10(0.05), linetype = 2, color = "darkgrey") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4),
        legend.position = c(0.88,0.13),
        legend.title = element_blank())+
  xlab(label = "Subject") +
  ylab(label = "log(Monte Carlo p-value)") +
  annotate("text", label = "p = 0.05", size = 3, x = 3, y = log10(0.04)) +
  scale_x_discrete(labels = gsub("MCTs", "", plot2_melt$UserName))




### Make boxplot for Figure 3B ###


colMeans(plot[1:6])
apply(plot[1:6], 2, function(x) median(x) )
table(apply(plot[1:6], 1, function(x) which.min(x)))

plot_melt <- melt(plot)
plot_melt <- na.omit(plot_melt)

source("../../../lib/colors/UserNameColors.R")

names(UserNameColors)<- gsub("MCTs", "", names(UserNameColors))

ggplot(plot_melt, aes(x = variable, y = -log10(value))) +
  geom_boxplot(outlier.shape = NA) + 
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "All with\ndecaying\ninfluence")) +
  geom_jitter(shape = 21, size = 3, fill = "light grey", alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "darkgrey") +
  theme_classic() +
  xlab(label = "Number of days combined") +
  ylab(label = "Monte-carlo p-value")

plot_melt$UserName <- gsub("MCTs", "", plot_melt$UserName)

myplot <- ggplot(plot_melt, aes(x = variable, y = -log10(value))) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "All with\ndecaying\ninfluence")) +
  geom_line(aes(group = UserName), color = "grey", alpha = 0.5) +
  geom_jitter(aes(color = UserName), size = 3, alpha = 0.7, width = 0.1) +
  scale_color_manual(values = UserNameColors) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "darkgrey") +
  geom_boxplot(outlier.shape = NA, fill = NA) + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "right") +
  guides(color = guide_legend(nrow = 5, title = "Subject", title.position = "top")) +
  xlab(label = "Number of days combined") +
  ylab(label = "-log(Monte Carlo p-value)")


myplot_leg <- get_legend(myplot) 

# and replot suppressing the legend
myplot_1 <- myplot + theme(legend.position='none')

myplot_1

ggsave("../../../output/Figure3/Figure3B.pdf", myplot_1, height = 2.5, width = 4)


#pvals_boxplot_legend, fig.height = 1.5, fig.width=4

fig3_legend <- myplot_leg


ggsave("../../../output/Figure3/Figure3_legend.pdf", fig3_legend, height = 1.5, width = 4)



### MAKE 3C venn, echo = F, include = T, fig.height=4, fig.width=4 ##

require(VennDiagram)
require(gridExtra)

venn_plot <- plot

venn_plot[1:6] <- apply(venn_plot[1:6], 2, function(x) ifelse(x < 0.05, as.character(venn_plot$UserName), NA))


xlist <- apply(venn_plot[c(2,3,4,6)], 2, function(x) x[!is.na(x)])

venn.plot <- venn.diagram(xlist , filename = NULL, fill=c("red", "orange", "#5a2071", "pink"), alpha=c(0.4,0.4,0.4,0.4), 
                          cex = 2, cat.fontface=4, category.names=c("2 days", "3 days", "4 days", "All days"))
 
# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

ggsave("../../../output/Figure3/Figure3C.pdf", venn.plot, height = 4, width = 4)
