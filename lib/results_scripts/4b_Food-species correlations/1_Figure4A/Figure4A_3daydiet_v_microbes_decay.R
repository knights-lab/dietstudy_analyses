# plot food v microbe per person
require(reshape2)

setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_decay")

# load datlist (out put from running correlations)
load(file = "../../../data/food_v_microbes_per_person_decay.RData")

# write out a file showing all the significant correlations from each person
#sigs <- lapply(datlist, function(x) subset(x, fdr_pval <= 0.9))
sigs <- lapply(datlist, function(x) subset(x, fdr_pval <= 0.1))
allsigs <- do.call("rbind", sigs)

# make a column with L1
allsigs$FoodL1 <- gsub(";L2_.*", "", allsigs$Food)
allsigs$FoodL1 <- gsub("L1_", "", allsigs$FoodL1)


fix.names.family <- function(x) {
  x <- gsub(";g__.*", "", x)
  x <- gsub("?.*f__", "", x)
  x <- gsub("?.*o__", "Uncl. Order", x)
  x <- gsub("?.*p__", "Uncl. Phylum ", x)
  x <- gsub("?.*k__", "Uncl. Kingdom ", x)
  x <- gsub(";NA", "", x)
  x <- gsub("\\[", "", x)
  x <- gsub("\\]", "", x)
  x <- gsub("-", "_", x)
  x <- gsub("\\/", "_", x)
  return(x)
}


fix.names.species <- function(x) {
  x <- gsub("?.*s__", "", x)
  x <- gsub("?.*g__", "Uncl. Genus ", x)
  x <- gsub("?.*f__", "Uncl. Family ", x)
  x <- gsub("?.*o__", "Uncl. Order ", x)
  x <- gsub("?.*p__", "Uncl. Phylum ", x)
  x <- gsub(";NA", "", x)
  return(x)
}

allsigs$Family<- fix.names.family(allsigs$Taxa)
allsigs$Family <- ifelse(allsigs$Family == "", allsigs$Taxa, allsigs$Family)
allsigs$Family <- gsub(";f.*", "", allsigs$Family)
allsigs$Family <- gsub("?.*o__", "Uncl. Order ", allsigs$Family)

allsigs$Species <- fix.names.species(allsigs$Taxa)

allsigs$bin <- ifelse(allsigs$Correlation < -0.65, "Negative (+/-)", 
                      ifelse(allsigs$Correlation > 0.65, "Positive (+/+ or -/-)", "NA"))
allsigs$bin <- factor(allsigs$bin, levels = c("Positive (+/+ or -/-)", "Negative (+/-)"))

allsigs$FoodL1 <- gsub("Dry_Beans_Peas_Other_Legumes_Nuts_and_Seeds", "Legumes", allsigs$FoodL1)
allsigs$FoodL1 <- gsub("Fats_Oils_and_Salad_Dressings", "Fats", allsigs$FoodL1)
allsigs$FoodL1 <- gsub("Grain_Product", "Grains", allsigs$FoodL1)
allsigs$FoodL1 <- gsub("Milk_and_Milk_Products", "Milks", allsigs$FoodL1)
allsigs$FoodL1 <- gsub("Meat_Poultry_Fish_and_Mixtures", "Meats", allsigs$FoodL1)
allsigs$FoodL1 <- gsub("Sugars_Sweets_and_Beverages", "Sweets and Beverages", allsigs$FoodL1)


# load colors
source(file = "../../../lib/colors/UserNameColors.R")

## get information to label key values repeated in more than one person for labeling
allsigs_names <- allsigs[colnames(allsigs) %in% c("Food", "Taxa")]
allsigs_names_pairs <- paste(allsigs_names$Food, allsigs_names$Taxa)

names_repeated<-as.data.frame(table(allsigs_names_pairs))
names_repeated<- subset(names_repeated, names_repeated$Freq >2)
table(names_repeated$Freq)
names_repeated <- colsplit(names_repeated$allsigs_names_pairs, " ", c("Food", "Taxa"))
names_repeated$label <- "*"




allsigs <- merge(allsigs, names_repeated, all.x = T)

allsigs$ID <- gsub("MCTs", "", allsigs$ID)

table(allsigs$ID)

length(table(allsigs$ID))


names(UserNameColors) <- gsub("MCTs", "", names(UserNameColors))

require(ggplot2)
myplot <- ggplot(data = allsigs, aes(x = Correlation, y = Family, size = -log10(fdr_pval), color = ID)) +
  geom_point(alpha = 0.8) +
  #geom_point(alpha = 0.8, color = "darkgrey", pch = 21) +
  facet_grid(FoodL1~bin, scales = "free", space = "free_y")+
  scale_color_manual(values = UserNameColors) +
  theme_classic() +
  guides(color = guide_legend(nrow = 10, title = "Subject", title.position = "top"),
         size = guide_legend(title.position = "top", title = "-log(FDR p-value)")) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 5, color = "black"),
        axis.text.x = element_text(size = 6, color = "black"),
        panel.grid.major = element_line(colour = "lightgrey"),
        strip.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8)) +
  ylab("Family-level taxonomy classification") +
  xlab("Spearman correlation") +
  scale_x_continuous(trans = "reverse")


#myplot

# L1_Meat_Poultry_Fish_and_Mixtures		#d43f1f
# L1_Grain_Product		#fe9700
# L1_Vegetables		#5dd047
# L1_Milk_and_Milk_Products	#00a2f2
# L1_Sugars_Sweets_and_Beverages	#c91acb
# L1_Fruits	#ffff59
# L1_Dry_Beans_Peas_Other_Legumes_Nuts_and_Seeds #662a00
# L1_Eggs	#a8863a
# L1_Fats_Oils_and_Salad_Dressings	#737373
sort(unique(allsigs$FoodL1))

library(grid)
g <- ggplot_gtable(ggplot_build(myplot))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("white", 
           "white",
           "#a8863a", #eggs
           "#7d7d7d", #fats
           "#ffff59",# fruits
           "#fe9700",#grains
           "#a14200",#legumes
           "#d43f1f",#meats
           "#00a2f2",#milk
           "#c91acb",#sugar
           "#5dd047"#vege
           )
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
#grid.draw(g)

ggsave("../../../output/Figure4/Figure4A_3daydietpersonalizedcor_decay.pdf",g, device = "pdf", width = 5, height = 9.5, dpi = 500)



# Fix labeling for export table
export <- subset(allsigs, allsigs$label == "*")
export$Food <- gsub(".*L2_", "", export$Food)
export <- export[colnames(export) %in% c("Food", "Species", "Correlation", "Pvalue", "fdr_pval", "ID")]
export <- export[order(export$Food, export$Species),]
export <- export[,c("Food", "Species", "Correlation", "Pvalue", "fdr_pval", "ID")]
export$Correlation <- round(export$Correlation,2)
export$Pvalue <- ifelse(export$Pvalue < 0.001, "< 0.001", round(export$Pvalue,3))
export$fdr_pval <- signif(export$fdr_pval, 2)
export$ID <- gsub("MCTs", "", export$ID)



# write.table(export, file = "../../../output/TableS2/Table S2_decay.txt",
#             sep = "\t",
#             col.names = T,
#             row.names = F,
#             quote = F)


# how many significant values in each person?
sigs.person <- unlist(lapply(sigs, function(x) dim(x)[1]))
length(sigs.person)
table(sigs.person>=1)
round(table(sigs.person>=1)[2]/length(sigs.person),2)
mean(sigs.person)
sd(sigs.person)       
median(sigs.person)
quantile(sigs.person, probs = c(5,95)/100)



