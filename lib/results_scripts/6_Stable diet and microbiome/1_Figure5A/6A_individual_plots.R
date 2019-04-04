
require(rmarkdown)
require(knitr)
require(tidyverse)
require(RColorBrewer) 
require(cowplot)
require(reshape2)
require(ggdendro)
require(vegan)
require(gridExtra)

# set the path for root directory
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/")

#### load the maps for downstream use ####
# full sample map - every possible sample
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# map for each subject
map_username <- read.table("data/maps/UserName_map.txt", sep = "\t", header = T, comment = "")

# map of just the smaples with taxonomy information after cleaning
tax_map <- read.table("data/maps/taxonomy_norm_map.txt", sep = "\t", header = T, comment = "")

# set up some of the colors palettes (http://phrogz.net/css/distinct-colors.html)
pal <- rev(c("#ff40f2", "#ff0000", "#008c4b", "#00138c", "#8c235b", "#ffbfbf", "#8c7723", "#468c75", "#8091ff", "#ff80c4", "#8c3123", "#fff2bf", "#40fff2", "#69698c", "#ff0044", "#ff9180", "#e5ff80", "#bffbff", "#5940ff", "#8c696e", "#8c7369", "#858c69", "#40d9ff", "#c480ff", "#ff8c40", "#4b8c00", "#23698c", "#69238c", "#8c4b00", "#bfffbf", "#004b8c", "#eabfff", "#ffc480", "#40ff59", "#80c4ff", "#ffd940" ))

pal4 <- c("#bf0000", "#f29d3d", "#a3d9b1", "#bfd9ff", "#f780ff", "#ff0000", "#ffe1bf", "#00d957", "#0020f2", "#e60099", "#730f00", "#7f6600", "#336655", "#293aa6", "#a6538a", "#8c4f46", "#e5c339", "#00ffcc", "#333a66", "#40202d", "#f29979", "#fbffbf", "#00b3a7", "#8091ff", "#cc335c", "#594943", "#a3d936", "#003033", "#300059", "#400009", "#cc5200", "#354020", "#39c3e6", "#cbace6", "#d9a3aa", "#7f5940", "#6a8040", "#006fa6", "#cc00ff", "#402910", "#0a4d00", "#566573", "#83008c")


# Determining clustering for order
# import average microbime per person for clustering
# here we use the normalized values for determing clustering from Bray Curtis distances
# These are used over the compositionally corrected relative abundances because they make more sense visually
microbes <- read.delim("data/microbiome/processed_average/UN_taxonomy_norm_s.txt", row = 1)

# create normalized sqrt relative abundance
# cluster by sqrt becuase it takes better into account lower abundance microbes
# also results in nice ordering
microbes.sq <- sweep(sqrt(microbes), 2, colSums(sqrt(microbes)), "/")


# create distance matrix for clustering
# using bray becuase it works well for grouping the prevotella dominant people on one side of the dendrogram
mydist <- vegdist(t(microbes.sq), method = "bray") 
myclust <- hclust(mydist, method = "average")
order <- myclust$order
microbes <- microbes[,order]
ord_factor <- as.character(colnames(microbes))

###### MICROBES aim for width 7.5, height 2.6########

# import microbiome data
taxa <- read.table("data/microbiome/processed_sample/taxonomy_norm_s.txt", header = T, row.names = 1, sep = "\t", comment = "")

#readjust for plotting
#taxa = sweep(sqrt(taxa),2,colSums(sqrt(taxa)),'/')
taxa = sweep(taxa, 2, colSums(taxa),'/')

# for plotting, limit to top ~20 bugs
taxa = taxa[rowMeans(taxa) >= 0.0047,]

# sort by highest average relative abundance
taxa <- taxa[order(rowMeans(taxa), decreasing = F),]

# make tax ordering factor
tax_ord_factor <- as.character(rownames(taxa))

plot <- as.data.frame(t(taxa))
plot <- rownames_to_column(plot, var = "X.SampleID")
plot <- melt(plot, id = "X.SampleID", variable.name = "Species")

# Check there are no duplicates/sum all of the same species for the same person
plot <- plot %>% group_by(X.SampleID, Species) %>% dplyr::summarise(newvalue = sum(value))

# create other category for plotting
other <- plot %>% group_by(X.SampleID) %>% dplyr::summarise(newvalue = 1- sum(newvalue))
other$Species <- "Other"
other <- other %>% select(X.SampleID, Species, newvalue)

# before we can join the two dfs with rbind must change from tbl dataframe caused by dplyr
plot <- as.data.frame(plot)
other <- as.data.frame(other)

# join for plotting
plot <- as.data.frame(rbind(plot,other))

# fix labeling for plotting
plot$Species <- gsub("?.*s__", "", plot$Species)
plot$Species <- gsub("?.*g__", "Uncl. Genus ", plot$Species)
plot$Species <- gsub("?.*f__", "Uncl. Family ", plot$Species)
plot$Species <- gsub("?.*o__", "Uncl. Order ", plot$Species)
plot$Species <- gsub(";NA", "", plot$Species)

# merge with map to get day variable
plot <- merge(plot, map, by = "X.SampleID")
plot$StudyDayNo <- as.numeric(plot$StudyDayNo)

# get right number of colors for plotting
no_cols <- length(unique(plot$Species))
colors_micro <- c("lightgrey", "#afd35a", "#ec9bfa", "#3c3e00", "#aaa4e1", "#f4f4f4", "#a717d9","#f4e909",
                  "#ffaa00", "#ff0000", "#69ef7b", "#d20000", "#ff5084", "#b04fbc", "#5bc07f",
                  "#e242a8", "#ffaae8", "#df97ac", "#ffede7", "#0024ff", "#ffa787", "#71619c",
                  "#668e57", "#3d7600", "#1f72d4", "#12982d", "#dd6c43", "#c25945", "#b6601e",
                  "#01ceff", "#02531d", "#d7dd4c", "#8f0000", "#4b03a9", "#A0447F")

# order participants by their microbiome diversity
# reorder by ordered factor
plot$UserName <- gsub("MCTs", "", plot$UserName)
ord_factor <- gsub("MCTs", "", ord_factor)
plot$UserName <- factor(plot$UserName, levels = ord_factor)

# reorder food species and clean up naming for the ordering
tax_ord_factor <- gsub("?.*s__", "", tax_ord_factor)
tax_ord_factor <- gsub("?.*g__", "Uncl. Genus ", tax_ord_factor)
tax_ord_factor <- gsub("?.*f__", "Uncl. Family ", tax_ord_factor)
tax_ord_factor <- gsub("?.*o__", "Uncl. Order ", tax_ord_factor)
tax_ord_factor <- gsub(";NA", "", tax_ord_factor)
tax_ord_factor <- c("Other", tax_ord_factor)

plot$Species <- as.factor(plot$Species)
plot$Species <-factor(plot$Species, levels = tax_ord_factor)

# plot as area plot
micro_plot <- ggplot(data = plot, aes(x=StudyDayNo, y = newvalue, fill=Species)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_micro) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.x = element_text(angle = 0, size = 6, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_rect(color = "grey"), 
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        panel.spacing.x=unit(0.05, "lines")) +
  guides(fill = guide_legend(reverse = FALSE,
                             keywidth = 0.5,
                             keyheight = 0.5,
                             ncol = 1)) +
  #nrow = 5)) + # for full page
  #nrow = 3)) + # for slide
  ylab("Relative Abundance") 



####### DIET width 7.5 height 2.6##########

# import diet data
food <- read.table("data/diet/processed_food/dhydrt.txt", header = T, sep = "\t", comment = "", row.names = "taxonomy")
food <- food[, !(colnames(food) == "X.FoodID")]

# collapse at level 2
# Summarizing at different levels - makes changes to everything downstream
split <- strsplit(rownames(food),";")             # Split and rejoin on lv7 to get species level
foodStrings <- sapply(split,function(x) paste(x[1:1],collapse=";"))
food<- rowsum(food,foodStrings)              # Collapse by taxonomy name

# sqrt relative abundance
food.n <- sweep(sqrt(food),2,colSums(sqrt(food)),'/')

# drop to top 20 ish foods
food <- food[rowMeans(food.n) > 0.015,]

# how many unique foods?
# length(unique(rownames(food)))

# relative abundance
#food <- sweep(sqrt(food), 2, colSums(sqrt(food)), '/')
food <- sweep(food, 2, colSums(food), '/')

# sort by highest average relative abundance
food <- food[order(rowMeans(food), decreasing = F),]

# make food ordering factor
food_ord_factor <- as.character(rownames(food))

plot3 <- as.data.frame(t(food))
plot3 <- rownames_to_column(plot3, var = "X.SampleID")
plot3 <- melt(plot3, id = "X.SampleID", variable.name = "Food")

# combine all "<x% abundance" foods into one for plotting
plot3 <- plot3 %>% group_by(X.SampleID, Food) %>% dplyr::summarise(newvalue = sum(value))

# fix labeling for plotting
plot3$Food <- gsub(".*L1_", " ", plot3$Food)
plot3$Food <- gsub("_", " ", plot3$Food)

# merge with map to get day variable
plot3 <- merge(plot3, map, by = "X.SampleID")

# set seed to get nice colors
set.seed(3)

# get right number of colors for plotting
no_cols <- length(unique(plot3$Food))


##TODO fix the color of meats/milks
colors_food <- rev(c("#fe9700",
                     "#c91acb",
                     "#d43f1f",
                     "#00a2f2",
                     "#5dd047",
                     "#ffff59",
                     "#662a00",
                     "#a8863a",
                     "#737373"))

# order participants by their microbiome diversity
# reorder by ordered factor and rename
plot3$UserName <- gsub("MCTs", "", plot3$UserName)
ord_factor <- gsub("MCTs", "", ord_factor)
plot3$UserName <- factor(plot3$UserName, levels = ord_factor)

# reorder food levels
food_ord_factor <- gsub(".*L1_", " ", food_ord_factor)
food_ord_factor <- gsub("_", " ", food_ord_factor)
plot3$Food <- as.factor(plot3$Food)
plot3$Food <-factor(plot3$Food, levels = food_ord_factor)

# make the plot
food_plot<-ggplot(data = plot3, aes(x=StudyDayNo, y = newvalue, fill=Food)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_food) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_blank(), 
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        panel.spacing.x=unit(0.05, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, 
                             keywidth = 0.5, 
                             keyheight = 0.5, 
                             ncol = 1)) +
  #nrow = 5)) + # for full page figure
  #nrow = 1)) + # for slide figure
  ylab("Relative Abundance") 


# Expanded plots on a per-person basis
# single_microbes, fig.width=5, fig.height= 3.5
 
single_plots <- as.list(NULL)

plot <- plot[order(plot$UserName),]

for (i in unique(plot$UserName)) {
  plot_un <- plot[plot$UserName == i,]
  
  single_plots[[i]] <- ggplot(data = plot_un, aes(x=StudyDayNo, y = newvalue, fill=Species)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_micro) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.x = element_text(angle = 0, size = 6, face = "italic"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "grey"), 
        legend.position = "none")

}



#single_foods, fig.width=5, fig.height= 3.5}

single_foods <- as.list(NULL)

plot3 <- plot3[order(plot3$UserName),]

for (i in unique(plot3$UserName)) {
  plot3_un <- plot3[plot3$UserName == i,]
  
  single_foods[[i]] <- ggplot(data = plot3_un, aes(x=StudyDayNo, y = newvalue, fill=Food)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_food) +
  scale_x_discrete(drop = TRUE) +
  theme_classic() +
  theme(strip.text.x = element_text(angle = 0, size = 6, face = "italic"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(color = "grey"), 
        legend.position = "none") 
}



# Area plots subset for figure 5A
#soylent_10, echo=F, fig.height=4, fig.width=5}

mygrobs <- c(single_plots[c(22,23,18)], single_foods[c(22,23,18)])

fig6A <- arrangeGrob(grobs = mygrobs, nrow=2)

ggsave("output/Figure6/Figure6A.pdf", fig6A, width = 5, height = 4)




micro_plot_leg <- get_legend(micro_plot) 
food_plot_leg <- get_legend(food_plot) 

smallplotlegend <- plot_grid(micro_plot_leg, food_plot_leg, ncol =1, align = "v")

save_plot("output/Figure6/Figure6A_legend.pdf", smallplotlegend, base_width = 3, base_height = 8)
