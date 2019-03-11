---
title: "Figure 1"
author: "Abigail Johnson"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output: 
  pdf_document:
    dev: png
---




```{r setup, include=FALSE, message = FALSE}

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
opts_knit$set(root.dir = "/Users/abby/Documents/Projects/dietstudy/")

# create a directory for the figure and set it's resolution/format
opts_chunk$set(echo = TRUE, fig.path = "Figs_Jan2018/", dev = c("pdf"), dpi = 300)
```

```{r load long map, echo = FALSE}
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


```

```{r dend_plot, echo = F, fig.width = 10, fig.height=4, message=F, include = F}

# Determining clustering for order
# import average microbime per person for clustering
# here we use the normalized values for determing clustering from Bray Curtis distances
# These are used over the compositionally corrected relative abundances because they make more sense visually
microbes <- read.delim("data/processed_UN_tax/UN_taxonomy_norm_s.txt", row = 1)

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


```



```{r plot microbime, include =F, fig.width = 7.5, fig.height=2.6, message = F}
# import microbiome data
taxa <- read.table("data/processed_tax/taxonomy_norm_s.txt", header = T, row.names = 1, sep = "\t", comment = "")

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

#names(colors_micro) <- levels(plot$Species) 


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




```




```{r plot_function, include =F, fig.width = 7.5, fig.height=2.6, message = F}
#path_map <- read.table("raw/functional/Kegg_ID_Map_plus.txt", sep = "\t", header = T, comment = "")

# import microbiome data for function
modules <- read.table("data/processed_KEGG/preprocessed_KEGG_modules.txt_filtered", header = T, row.names = 1, sep = "\t", comment = "")

#readjust for plotting (don't know how to do with CLR data)
#modules = sweep(sqrt(modules),2,colSums(sqrt(modules)),'/')
modules = sweep(modules,2,colSums(modules), '/')

# sort by highest average relative abundance
modules <- modules[order(rowMeans(modules), decreasing = F),]

# make tax ordering factor
mod_ord_factor <- as.character(rownames(modules))

plot2 <- as.data.frame(t(modules))
plot2 <- rownames_to_column(plot2, var = "X.SampleID")
plot2 <- melt(plot2, id = "X.SampleID", variable.name = "Modules")

# combine all "<x% abundance" foods into one for ploting
plot2 <- plot2 %>% group_by(X.SampleID, Modules) %>% dplyr::summarise(newvalue = sum(value))


# merge with map to get day variable
plot2 <- merge(plot2, map, by = "X.SampleID")

set.seed(1)

# get right number of colors for plotting
no_cols <- length(unique(plot2$Modules))
#colors_func <- sample(cols_func(no_cols))
colors_func <- sample(pal, no_cols)

# order participants by their microbiome diversity
# reorder by ordered factor
plot2$UserName <- gsub("MCTs", "", plot2$UserName)
plot2$UserName <- factor(plot2$UserName, levels = ord_factor)

# reorder modules
plot2$Modules <- as.factor(plot2$Modules)
plot2$Modules <-factor(plot2$Modules, levels = mod_ord_factor)

# plot2 <- merge(plot2, path_map, by.x = "Modules", by.y = "Shogun_ID")
# plot2$KEGG_Label <- gsub(" \\[.*$", "", plot2$KEGG_Label)
# plot2$KEGG_Label <- gsub("M..... ", "", plot2$KEGG_Label)

soylents2 <- subset(plot2,UserName == c('11', '12'))
soylents2 <- soylents2[!duplicated(soylents2$UserName),]

# make the plot2
func_plot <- ggplot(data = plot2, aes(x=StudyDayNo, y = newvalue, fill=Modules)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_func) +
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
                             #nrow = 5)) +
  ylab("Relative Abundance") 



```

```{r plot diet, include =F, fig.width = 7.5, fig.height=2.6, message = F}
# import diet data
food <- read.table("data/processed_food/dhydrt.txt", header = T, sep = "\t", comment = "", row.names = "taxonomy")
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
"#00a2f2",
"#d43f1f",
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

soylents3 <- subset(plot3,UserName == c('11', '12'))
soylents3 <- soylents3[!duplicated(soylents3$UserName),]


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


```


```{r plot nutr, include =F, fig.width = 7.5, fig.height=2.6, message = F}

# load nutrition intake for each person
nutr <- as.data.frame(t(read.delim("data/processed_nutr/nutr_65.txt", row.names = 1, sep = "\t", comment = "")))


# this table already has everything in grams so it's ready to use to plot normalized grams for macro and micronutrients

plot4 <- nutr %>% select(-c(KCAL, MOIS,  
                            SFAT, S040, S060, S080, S100, S120, S140, S160, S180,
                            MFAT, M161, M181, M201, M221, 
                            PFAT, P182, P183, P184, P204, P205, P225, P226))
plot4$CARB <- plot4$CARB-plot4$SUGR
plot4 <- t(plot4)


# normalize
plot4 <- sweep(sqrt(plot4),2,colSums(sqrt(plot4)),'/')
#plot4 <- sweep(plot4,2,colSums(plot4),'/')

# load data dictionary to change names
dict <- read.delim("raw/Uncleaned_diet/Data_Dictionary/DataDictionary_ITEMS.txt", stringsAsFactors = F)
dict <- dict[colnames(dict) == c("Field.Name", "Description")]
dict$Field.Name <- gsub(" ", "", dict$Field.Name )
dict <- dict[(dict$Field.Name) %in% rownames(plot4),]

rownames(plot4) <- dict$Description
rownames(plot4) <- gsub(" \\(mcg\\) ", "", rownames(plot4))
rownames(plot4) <- gsub(" \\(mg\\) ", "", rownames(plot4))
rownames(plot4) <- gsub(" \\(g\\) ", "", rownames(plot4))
rownames(plot4) <- gsub(" \\(mg\\)", "", rownames(plot4))
rownames(plot4) <- gsub(" \\(mcg\\)", "", rownames(plot4))
rownames(plot4) <- gsub(" \\(mcg_RAE\\) ", "", rownames(plot4))
rownames(plot4) <- gsub(" \\(mcg_DFE\\)", "", rownames(plot4))
rownames(plot4) <- gsub(", DFE", "", rownames(plot4))
rownames(plot4) <- gsub(", RAE", "", rownames(plot4))

# sort by highest average relative abundance
plot4 <- plot4[order(rowMeans(plot4), decreasing = F),]

plot4 <- as.data.frame(t(plot4))



plot4 <- rownames_to_column(plot4, var = "X.SampleID")

plot4 <- melt(plot4, id = "X.SampleID", variable.name = "Nutrient")


# combine all "<x% abundance" foods into one for plotting
plot4 <- plot4 %>% group_by(X.SampleID, Nutrient) %>% dplyr::summarise(newvalue = sum(value))


# merge with map to get day variable
plot4 <- merge(plot4, map, by = "X.SampleID")


# order participants by their microbiome diversity
# reorder by ordered factor and rename
plot4$UserName <- gsub("MCTs", "", plot4$UserName)
ord_factor <- gsub("MCTs", "", ord_factor)
plot4$UserName <- factor(plot4$UserName, levels = ord_factor)

soylents4 <- subset(plot4,UserName == c('11', '12'))
soylents4 <- soylents4[!duplicated(soylents4$UserName),]


no_cols <- length(unique(plot4$Nutrient))

set.seed(2)
colors_nutr <-sample(pal4)[1:no_cols]

# plot kcal filled with % protein, fat, and carb
nutr_plot <- ggplot(data = plot4, aes(x = StudyDayNo, y = newvalue , fill = Nutrient)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values= colors_nutr) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_blank(),
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
  ylab("sqrt(Relative Abundance)")




```


```{r figure_1, echo = F, fig.width=8, fig.height =8.4}
# combine into one big plot

micro_plot_leg <- get_legend(micro_plot) 
func_plot_leg <- get_legend(func_plot) 
food_plot_leg <- get_legend(food_plot) 
nutr_plot_leg <- get_legend(nutr_plot) 

# and replot suppressing the legend
micro_plot_1 <- micro_plot + theme(legend.position='none')
func_plot_1 <- func_plot + theme(legend.position='none')
food_plot_1 <- food_plot + theme(legend.position='none')
nutr_plot_1 <- nutr_plot + theme(legend.position='none')


plot_grid(micro_plot_1, func_plot_1, food_plot_1, nutr_plot_1,  
          ncol = 1, 
          scale = c(1, 1, 1, 1), 
          rel_heights = c(1, 1, 1, 1),
          labels = "AUTO",
          align = "h")

```





```{r legend, echo = F, fig.width = 13, fig.height = 5}

plot_grid(micro_plot_leg, func_plot_leg, food_plot_leg, nutr_plot_leg, nrow =1, align = "h")


```

# Expanded plots on a per-person basis
```{r single_microbes, include = F, fig.width=5, fig.height= 3.5}
 
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

```


```{r multiplot_microbes, echo=F, fig.height=11, fig.width=8.5}

grid.arrange(grobs = single_plots[1:34], nrow=7)

```

```{r single_func, include = F, fig.width=5, fig.height= 3.5}

single_funcs <- as.list(NULL)

plot2 <- plot2[order(plot2$UserName),]

for (i in unique(plot2$UserName)) {
  plot2_un <- plot2[plot2$UserName == i,]
  
  single_funcs[[i]] <- ggplot(data = plot2_un, aes(x=StudyDayNo, y = newvalue, fill=Modules)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_func) +
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

```

```{r multiplot_func, echo=F, fig.height=11, fig.width=8.5}

grid.arrange(grobs = single_funcs[1:34], nrow=7)

```


```{r single_nutr, include = F, fig.width=5, fig.height= 3.5}

single_nutr <- as.list(NULL)

plot4 <- plot4[order(plot4$UserName),]

for (i in unique(plot4$UserName)) {
  plot4_un <- plot4[plot4$UserName == i,]
  
  single_nutr[[i]] <- ggplot(data = plot4_un, aes(x=StudyDayNo, y = newvalue, fill=Nutrient)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_nutr) +
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

```

```{r multiplot_nutr, echo=F, fig.height=11, fig.width=8.5}

grid.arrange(grobs = single_nutr[1:34], nrow=7)

```



```{r single_foods, include = F, fig.width=5, fig.height= 3.5}

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

```

```{r multiplot_foods, echo=F, fig.height=11, fig.width=8.5}

grid.arrange(grobs = single_foods[1:34], nrow=7)

```



```{r soylentfoods_11_12, echo = F, fig.width=7, fig.height= 2, include = F}
# For response to reviewers
#### single food but level 6 ####
# import diet data
food <- read.table("data/processed_food/dhydrt.txt", header = T, sep = "\t", comment = "", row.names = "taxonomy")
food <- food[, !(colnames(food) == "X.FoodID")]

mcts <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
food <- food[, colnames(food) %in% mcts$X.SampleID]


# relative abundance
food <- sweep(food, 2, colSums(food), '/')

# drop to top 30 ish foods
food <- food[rowMeans(food) > 0.00001,]

# fix names
rownames(food) <- gsub(".*;", "", rownames(food))

# sort by highest average relative abundance
food <- food[order(rowMeans(food), decreasing = F),]

plot6 <- as.data.frame(t(food))
plot6 <- rownames_to_column(plot6, var = "X.SampleID")
plot6 <- melt(plot6, id = "X.SampleID", variable.name = "Food")


# merge with map to get day variable
plot6 <- merge(plot6, map, by = "X.SampleID")


# set seed to get nice colors
set.seed(3)

# get right number of colors for plotting
colors_food_soylent <- c("#775e42",
                 "#fe9700",
                 "#961ACB",
                 "#c91acb",
                 "#ff93ff",
                 "#a100aa",
                 "#dcabda")



single_food_plot_6<-ggplot(data = plot6, aes(x=StudyDayNo, y = value, fill=Food)) +
  geom_area(stat = "identity") +
  facet_grid(.~UserName, scales = "free") +
  scale_fill_manual(values = colors_food_soylent) +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        panel.spacing.x=unit(0.05, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, 
                             keywidth = 0.5, 
                             keyheight = 0.5, 
                             ncol = 1)) +
                             #nrow = 5)) + # for full page figure
                             #nrow = 1)) + # for slide figure
  ylab("Relative Abundance") +
 xlab("")


single_food_plot_6



```
```{r soylent_nutr, echo=F, fig.height=2, fig.width=3, include = F}

mygrobs <- c(single_nutr[c(22,23)])

grid.arrange(grobs =mygrobs, nrow=1)

```

# Area plots subset for figure 5A
```{r soylent_10, echo=F, fig.height=4, fig.width=5}

mygrobs <- c(single_plots[c(22,23,18)], single_foods[c(22,23,18)])

grid.arrange(grobs =mygrobs, nrow=2)

```


