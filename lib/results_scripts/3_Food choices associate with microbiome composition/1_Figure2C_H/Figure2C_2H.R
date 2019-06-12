
require(rmarkdown)
require(knitr)
require(tidyverse)
require(RColorBrewer) 
require(cowplot)
require(reshape2)
require(ggdendro)
require(vegan)
require(ape)

# set the path for root directory
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/")

# set seed
set.seed(42)


#### load the maps for downstream use ####
map_sample <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
map_username <- read.table("data/maps/UserName_map.txt", sep = "\t", header = T, comment = "")
map_sample$StudyDate <- as.Date.factor(map_sample$StudyDate, format = "%m/%d/%y")

#load other maps
food_map <- read.table("data/maps/food_map.txt", sep = "\t", header = T, comment = "")
tax_map <- read.table("data/maps/taxonomy_norm_map.txt", sep = "\t", header = T, comment = "")


##########food##############
# load the food distance matrix, unweighted unifrac
food_un <- read.delim("data/diet/processed_food/dhydrt_smry_no_soy_beta/unweighted_unifrac_dhydrt.smry.no.soy.txt", row = 1) # weighted in not significant
food_dist <- as.dist(food_un)

# make non-tree distance matrix for food (for supplemental)
food <- read.delim("data/diet/processed_food/dhydrt.smry.no.soy.txt", row = 1)
food <- food[,colnames(food) %in% colnames(food_un)]
no_tree_dist <- dist(t(food))


##############taxonomy############
# load taxonomy collapsed for each person
tax <- read.delim("data/microbiome/processed_average/UN_tax_CLR_mean_norm_s.txt", row = 1) # updates from reviewer comments
# drop soylents
tax <- tax[,colnames(tax) %in% colnames(food_un)]

# make taxonomy distance matrix
tax_dist <- dist(t(tax))


###############nutrients############
# load nutrition data
nutr <- read.delim("data/diet/processed_nutr/nutr_65_smry_no_soy.txt", row = 1)

# normalize nutrition data across features (rows)
nutr_n <- sweep(nutr, 1, rowSums(nutr), "/")

# make nutrition distance matrix (euclidean)
nutr_dist <- dist(t(nutr_n))

# identify soylent samples
soylent <- map_sample[map_sample$UserName %in% c("MCTs11", "MCTs12"),]
soylent <- as.character(soylent$X.SampleID)



source(file = "lib/colors/UserNameColors.R")

tax_all <- read.delim("data/microbiome/processed_sample/taxonomy_clr_s.txt", row = 1)

# drop soylent 
tax_all_ns <- tax_all[, !colnames(tax_all) %in% soylent]

#make dist matrix
tax_all_dist <- dist(t(tax_all_ns))

# make pcoa
tax_all_pcoa <- data.frame(pcoa(tax_all_dist)$vectors)

eigen <- pcoa(tax_all_dist)$values$Eigenvalues

percent_var <- signif(eigen/sum(eigen), 4)*100

# move rownames to a column
tax_all_pcoa <- rownames_to_column(tax_all_pcoa, var = "X.SampleID")

#merge map and PCOA by SampleID - select both lines and run together - won't work otherwise
tax_all_pcoa <- inner_join(tax_all_pcoa, map_sample, by = 'X.SampleID')

## PERMANOVA for within subject grouping

mytest <- adonis(tax_all_dist~tax_all_pcoa$UserName, permutations = 999)
plot_p <- mytest$aov.tab$`Pr(>F)`[1]

mytest


#### MAKE 2C ######


tax_all_plot_one_color <- ggplot(tax_all_pcoa, aes(x = Axis.1, y = Axis.2, fill = UserName)) +
  geom_point(color =  "#5f86b7", alpha=0.75, size = 2, alpha = 0.6) +
  stat_ellipse(color = "dark grey", type = "norm", linetype = 2, level = 0.95, alpha = 0.5) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) +
  theme_classic() +                                                                      
    theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 4),
        axis.title = element_text(size=9)) +
  guides(color = guide_legend(ncol = 4))
  #annotate("text", x = 15, y = -30, label = paste0("p-value = ",plot_p))

tax_all_plot_one_color + theme(legend.position = 'none')



food_all <- read.delim("data/diet/processed_food/dhydrt_beta/unweighted_unifrac_dhydrt.txt", row = 1)


# drop soylent 
food_all_ns <- food_all[!rownames(food_all) %in% soylent, !colnames(food_all) %in% soylent]
#food_all_dist <- as.dist(food_all_ns) #for adonis testing
food_all_pcoa <- data.frame(pcoa(food_all_ns)$vectors)

eigen <- pcoa(food_all_ns)$values$Eigenvalues

percent_var <- signif(eigen/sum(eigen), 4)*100

# move rownames to a column
food_all_pcoa <- rownames_to_column(food_all_pcoa, var = "X.SampleID")

#merge map and PCOA by SampleID
food_all_pcoa <- inner_join(food_all_pcoa, map_sample, by = 'X.SampleID')

## PERMANOVA for within subject grouping
#check order
#colnames(tax_all_ns) == tax_all_pcoa$X.SampleID

#mytest <- adonis(food_all_dist~food_all_pcoa$UserName, permutations = 999)
#plot_p <- mytest$aov.tab$`Pr(>F)`[1]

#mytest

### MAKE 2D ######

food_all_plot_one_color <- ggplot(food_all_pcoa, aes(x = Axis.1, y = Axis.2, fill = UserName)) +
  geom_point(color = "#5a2071", alpha=0.75, size = 2, alpha = 0.6) +
  stat_ellipse(color = "dark grey", type = "norm", linetype = 2, level = 0.95, alpha = 0.5) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) +
  theme_classic() +                                                                      
    theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 4),
        axis.title = element_text(size=9)) +
  guides(color = guide_legend(ncol = 1))
  #annotate("text", x = 15, y = -30, label = paste0("p-value = ",plot_p))

food_all_plot_one_color + theme(legend.position = 'none')

#### MAKE 2E #####

# make pcoas 
pcoa_f <- as.data.frame(pcoa(food_dist)$vectors)
pcoa_t <- as.data.frame(pcoa(tax_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_f, pcoa_t)
pro_test <- protest(pcoa_f, pcoa_t, perm = 999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Food (Unweighted Unifrac)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome (Aitchison's)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 1)

plot <- rbind(beta_pro, trans_pro)

food_micro <- ggplot(plot) +
  geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#5a2071", "#5f86b7")) +
    theme_classic() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
         aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.3, y = -0.27, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 
  

food_micro_leg <- get_legend(food_micro)


food_micro + theme(legend.position = "none")


## check for gender diff in average microbiome

map_username_sub <- map_username[map_username$UserName %in% rownames(pcoa_t),]
gentest <- adonis(tax_dist~map_username_sub$Gender, permutations = 999)
gen_p <- gentest$aov.tab$`Pr(>F)`[1]

gentest

## same for foods
foodtest <- adonis(food_dist~map_username_sub$Gender, permutations = 999)
foodtest

#### MAKE 2F #####
# make pcoas 
pcoa_n <- as.data.frame(pcoa(nutr_dist)$vectors)
pcoa_t <- as.data.frame(pcoa(tax_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_n, pcoa_t)
pro_test <- protest(pcoa_n, pcoa_t, perm = 999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Nutrient (Euclidean)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome (Aitchison's)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 1)

plot <- rbind(beta_pro, trans_pro)

nutr_micro <- ggplot(plot) +
  geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#5f86b7", "#2bbaa7")) +
    theme_classic() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
         aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.01, y = -0.19, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 

nutr_micro + theme(legend.position = "none")

nutr_micro_leg <- get_legend(nutr_micro)



#### MAKE 2G #####


# Import the beta diversity table from grains_beta (weighted or unweighted)
food_beta <- read.table(file="data/diet/fiber/grains_data/grains_beta/unweighted_unifrac_grains_fiber.txt")

food_dist <- as.dist(food_beta)
# refer this code to make a pcoa you will have to fix naming and play with the data structure of the dataframe called plot to get this to work.


# Import the beta diversity table for taxonomy
taxa_beta <- read.table(file="data/microbiome/processed_average/UN_tax_beta/euclidean_UN_taxonomy_clr_s.txt")

taxa_beta <- taxa_beta[rownames(taxa_beta) %in% rownames(food_beta),colnames(taxa_beta) %in% colnames(food_beta)]

tax_dist <- as.dist(taxa_beta)




# make pcoas 
pcoa_f <- as.data.frame(pcoa(food_dist)$vectors)
pcoa_t <- as.data.frame(pcoa(tax_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_f, pcoa_t)
pro_test <- protest(pcoa_f, pcoa_t, perm = 999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Grain Fiber (Unweighted Unifrac)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome (Aitchison's)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif)

plot <- rbind(beta_pro, trans_pro)

grain_micro <- ggplot(plot) +
  geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#fe9700", "#5f86b7")) +
    theme_classic() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
         aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.05, y = -0.27, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 

grain_micro + theme(legend.position = "none")

grain_leg <- get_legend(grain_micro)


### MAKE 2H #####
# Import the beta diversity table from fruit_beta (weighted or unweighted)
food_beta <- read.table(file="data/diet/fiber/fruit_data/fruit_beta/unweighted_unifrac_fruit_fiber.txt")

food_dist <- as.dist(food_beta)



# make pcoas 
pcoa_f <- as.data.frame(pcoa(food_dist)$vectors)
pcoa_t <- as.data.frame(pcoa(tax_dist)$vectors)

# procrustes
pro <- procrustes(pcoa_f, pcoa_t)
pro_test <- protest(pcoa_f, pcoa_t, perm = 999)

eigen <- sqrt(pro$svd$d)
percent_var <- signif(eigen/sum(eigen), 4)*100

beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
beta_pro$UserName <- rownames(beta_pro)
beta_pro$type <- "Fruit Fiber (Unweighted Unifrac)"
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Microbiome (CLR-Euclidean)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 3)

plot <- rbind(beta_pro, trans_pro)

fruit_micro <- ggplot(plot) +
    geom_point(size = 2, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
  scale_color_manual(values = c("#CBD13E", "#5f86b7")) +
    theme_classic() +
    geom_line(aes(x= Axis.1, y=Axis.2, group=UserName), col = "darkgrey", alpha = 0.6) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = 'bottom',
        axis.text = element_text(size=4),
        axis.title = element_text(size=9),
         aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  annotate("text", x = 0.33, y = -0.33, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 

fruit_micro + theme(legend.position = "none")

fruit_leg <- get_legend(fruit_micro)



plotC_H <- plot_grid(tax_all_plot_one_color + theme(legend.position = "none"), 
          food_all_plot_one_color + theme(legend.position = "none"), 
          food_micro + theme(legend.position = "none"), 
          nutr_micro + theme(legend.position = "none"), 
          grain_micro + theme(legend.position = "none"), 
          fruit_micro + theme(legend.position = "none"), 
          ncol = 2, 
          align = "h", 
          labels = c("C", "D", "E", "F", "G", "H"))

save_plot("output/Figure2/Figure2C_H.pdf", plotC_H, base_width = 3.5, base_height = 5.5)


plotC_H_leg <- plot_grid(food_micro_leg, nutr_micro_leg, grain_leg, fruit_leg, ncol = 1, align = "h")

save_plot("output/Figure2/Figure2C_H_legend.pdf", plotC_H_leg, base_width = 3, base_height = 2)
