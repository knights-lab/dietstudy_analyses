# use phILR
setwd("~/Documents/Projects/dietstudy/")
# source("https://bioconductor.org/biocLite.R")
# biocLite("philr")

library(philr)
library(phyloseq)
library(ape)
library(ggplot2)

## so to repeat with my data, need to force my data to be a phyloseq object
# convert map to phyloseq obect
map <- read.table("data/maps/UserName_map.txt", sep = "\t", header = T, comment = "", row.names = 1)
map <- map[!rownames(map) %in% c("MCTs11", "MCTs12"),]

map_seq <- sample_data(map)

# rename OTUs as OTU`X` to allow us to match separately with Taxonomy
otu <- read.delim("data/greengenes/UN_strain.processed.txt", row.names = 1)
otu <- otu[,!colnames(otu)%in% c("MCTs11", "MCTs12"),]

otu <- as.matrix(otu)

# make phyloseq otu file
otu_seq <- otu_table(otu, taxa_are_rows = TRUE)

# add tree
tree <- read_tree_greengenes("data/greengenes/gg97.tre") # exported from iTOL
#tree <- read.tree('data/greengenes/gg97.tre.cleaned') #from phyloseq package
pruned_tree <- prune_taxa(taxa = rownames(otu), x = tree) # matches class of the argument 

pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')

## Make phyloseq object
datP <- merge_phyloseq(otu_seq, map_seq, pruned_tree)

#minimal process phILR based on example: https://bioconductor.org/packages/release/bioc/vignettes/philr/inst/doc/philr-intro.html
datP <-  filter_taxa(datP, function(x) sum(x > 2) > (0.2*length(x)), TRUE) # drop taxa without at least 2 counts in at least 20% of people
datP <-  filter_taxa(datP, function(x) sd(x)/mean(x) > 2.0, TRUE) # drop taxa with low variation
datP <- transform_sample_counts(datP, function(x) x+1)

# check number of taxa remaining
datP

# store otu names for later
otu_to_keep <- rownames(otu_table(datP))

# check tree
is.rooted(phy_tree(datP))
is.binary.tree(phy_tree(datP))

# look at how many taxa remain



#transpose for phILR
otu.table <- t(otu_table(datP))
tree <- phy_tree(datP)
metadata <- sample_data(datP)



dat.philr <- philr(otu.table, tree, 
                  part.weights='enorm.x.gm.counts', 
                  ilr.weights='blw.sqrt')


dat.philr[1:5, 1:5] # lots fo columns with NA only, so remove them?

dat.philr <- dat.philr[,colSums(is.na(dat.philr))<nrow(dat.philr)]

dat.dist <- dist(dat.philr, method="euclidean")
dat.pcoa <- pcoa(dat.dist)

# load the food distance matrix, unweighted unifrac
food_un <- read.delim("data/processed_food/dhydrt_smry_no_soy_beta/unweighted_unifrac_dhydrt.smry.no.soy.txt", row = 1) # weighted in not significant
food_dist <- as.dist(food_un)

# make non-tree distance matrix for food (for supplemental)
food <- read.delim("data/processed_food/dhydrt.smry.no.soy.txt", row = 1)
food <- food[,colnames(food) %in% colnames(food_un)]
no_tree_dist <- dist(t(food))


# now see if these disances pair with the food distances
set.seed(7)

# make pcoas 
pcoa_f <- as.data.frame(pcoa(food_dist)$vectors)
pcoa_t <- dat.pcoa$vectors

require(vegan)
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
trans_pro$type <- "Microbiome (PhILR, gg97)"

colnames(trans_pro) <- colnames(beta_pro)

pval <- signif(pro_test$signif, 3)

plot <- rbind(beta_pro, trans_pro)

food_micro <- ggplot(plot) +
  geom_point(size = 3, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = type)) +
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
  annotate("text", x = 0.27, y = -0.27, label = paste0("p-value=",pval), size = 2) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) 



food_micro_pro <- food_micro + theme(legend.position = "none")

##### repeat for the big data set #####

## so to repeat with my data, need to force my data to be a phyloseq object
# convert map to phyloseq obect
map <- read.table("data/maps/taxonomy_norm_map.txt", sep = "\t", header = T, comment = "", row.names = 1)

map_seq <- sample_data(map)

# rename OTUs as OTU`X` to allow us to match separately with Taxonomy
otu <- read.delim("data/greengenes/strain.minimal.processed.txt", row.names = 1)

otu <- as.matrix(otu)

# filter to same otus as used for the collapsed data
otu <- otu[otu_to_keep,]

# make phyloseq otu file
otu_seq <- otu_table(otu, taxa_are_rows = TRUE)

# add tree
tree <- read_tree_greengenes("data/greengenes/gg97.tre") # exported from iTOL
#tree <- read.tree('data/greengenes/gg97.tre.cleaned') #from phyloseq package
pruned_tree <- prune_taxa(taxa = rownames(otu), x = tree) # matches class of the argument 

pruned_tree <- makeNodeLabel(pruned_tree, method="number", prefix='n')

## Make phyloseq object
datP <- merge_phyloseq(otu_seq, map_seq, pruned_tree)

#minimal process for philr

datP <- transform_sample_counts(datP, function(x) x+1)

datP

#transpose for phILR
otu.table <- t(otu_table(datP))
tree <- phy_tree(datP)
metadata <- sample_data(datP)



dat.philr <- philr(otu.table, tree, 
                   part.weights='enorm.x.gm.counts', 
                   ilr.weights='blw.sqrt')


dat.philr[1:5, 1:5] # lots fo columns with NA only, so remove them?

dat.philr <- dat.philr[,colSums(is.na(dat.philr))<nrow(dat.philr)]

dat.dist <- dist(dat.philr, method="euclidean")
dat.pcoa <- as.data.frame(pcoa(dat.dist)$vectors)

plot <- merge(dat.pcoa, map, by = 0)

eigen <- pcoa(dat.dist)$values$Eigenvalues

percent_var <- signif(eigen/sum(eigen), 4)*100

source(file = "lib/UserNameColors.R")

tax_all_plot <- ggplot(plot, aes(x = Axis.2, y = Axis.1, fill = UserName)) +
  geom_point(aes(color = UserName), alpha=0.75, size = 2.5, alpha = 0.7) +
  stat_ellipse(aes(color = UserName), type = "norm", linetype = 2, level = 0.95, alpha = 0.5) +
  scale_fill_manual(name = "Subject", values = UserNameColors) + 
  scale_color_manual(values = UserNameColors) +
  xlab(paste0("PC 1 [",percent_var[1],"%]")) +
  ylab(paste0("PC 2 [",percent_var[2],"%]")) +
  theme_classic() +                                                                      
  theme(legend.position = "none",
        legend.title.align = 0.5,
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 4),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = F, fill = guide_legend(ncol = 4)) 
#annotate("text", x = 15, y = -30, label = paste0("p-value = ",plot_p))
#theme(legend.title = element_text("Subject"))

tax_all_plot

require("cowplot")

myleg <- get_legend(food_micro)
# plot <- plot_grid(tax_all_plot + theme(legend.position = "none"), 
#                   food_micro + theme(legend.position = "none"), 
#                   ncol = 2, 
#                   align = "h")
# 
# save_plot("analysis/greengeens/Figure_phILR.pdf", plot, base_width = 4, base_height = 2.5)



ggsave("analysis/greengeens/FigureS3_gg1.pdf", tax_all_plot, height = 1.8, width = 1.8)
ggsave("analysis/greengeens/FigureS3_gg2.pdf", food_micro_pro, height = 1.8, width = 1.8)
ggsave("analysis/greengeens/FigureS3_gg2_leg.pdf", myleg, height =1, width = 2)
