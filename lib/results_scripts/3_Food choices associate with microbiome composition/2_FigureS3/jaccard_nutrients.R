homedir <- "/Users/abby/Documents/Projects/dietstudy/"
setwd(dir = homedir)
# assign 0/1 based on high or low value relative to the sample median

# get the median value for each nutrient from nhanes
require(RNHANES)
pop.nutrients <- nhanes_load_data("DR1TOT", "2011-2012")
nhanes.map <- nhanes_load_data("DEMO", "2011-2012")

# limit to just adults
nhanes.adult <- nhanes.map[nhanes.map$RIDAGEYR >= 18,]
pop.nutrients <- pop.nutrients[pop.nutrients$SEQN %in% nhanes.adult$SEQN,]


# fix names
colnames(pop.nutrients) <- gsub("DR1T", "", colnames(pop.nutrients))

# load per-person-day nutrient table
nutr <-read.delim(file = "data/processed_nutr/nutr_65.txt", row = 1)

#limit to values of interest
pop.nutrients <- pop.nutrients[,colnames(pop.nutrients) %in% rownames(nutr)]

# drop rows without information
pop.nutrients <- pop.nutrients[!is.na(rowSums(pop.nutrients)),]

pop.medians <- apply(pop.nutrients, 2, median)

# adjust so everything is in the same unit
# change most variables in the nutr_65 table into grams (i.e. mg to grams, etc)
# Prot, fat, carb, mois, alc, sugars,fiber, sfat, mfat, pfat
# and all the specific fatty acids are already in grams

# mg variables
mgvar <- c("CAFF", "THEO", "CALC", "IRON", "MAGN", "PHOS", "POTA", "SODI",
           "ZINC", "COPP", "VC", "VB1", "VB2", "NIAC", "VB6", "ATOC", "CHOLE",
           "CHOLN", "VITE_ADD")
mcgvar <- c("SELE", "FOLA", "FA", "FF", "FDFE", "VB12", "VARA", "RET", "BCAR", 
            "ACAR", "CRYP", "LYCO", "LZ", "VK", "VITD", "B12_ADD")

#mg values to grams
mgvals <- names(pop.medians) %in% c(mgvar)
pop.medians[mgvals] <- pop.medians[mgvals]/1000

#mcg to grams 
mcgvals <- names(pop.medians) %in% c(mcgvar)
pop.medians[mcgvals] <- pop.medians[mcgvals]/1000000



#limit nutrients table to just those that are in the population nutrients table
nutr <- nutr[rownames(nutr) %in% names(pop.medians),]

# reorder pop.medians 
pop.medians <- pop.medians[order(factor(names(pop.medians), levels = rownames(nutr)))]

nutr_test <- apply(nutr, 2, function(x) x > pop.medians)
nutr_test <- as.data.frame(nutr_test)
nutr_test <- sapply(nutr_test, as.numeric)  # this gives a matrix of 0s and 1s to use for jaccard distance.

#which(colSums(nutr_test) == 0) # no values with 0 colSums (good!)
# can't get meaningful distances for those days, so have to drop them
# since there aren't any, don't worry about this
#nutr_test <- nutr_test[,!colSums(nutr_test) == 0]


require(vegan)
jacc_dists <- vegdist(t(nutr_test), method = "jaccard")

require(ape)
jacc_pcoa <- as.data.frame(pcoa(jacc_dists)$vectors)    

require(ggplot2)

map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "", row.names = 1)
plot <- jacc_pcoa[1:2]
plot <- merge(plot, map, by = 0)

ggplot(data = plot, aes(x = Axis.1, y = Axis.2, color = UserName)) + geom_point() + stat_ellipse(level = 0.75)

#######################################

# repeat for summary data
#  load per-person summary nutrient table
summary_nutr <- read.delim(file = "data/processed_nutr/nutr_65_smry.txt", row = 1)


#limit nutrients table to just those that are in the population nutrients table
summary_nutr <- summary_nutr[rownames(summary_nutr) %in% names(pop.medians),]

nutr_summ_test <- apply(summary_nutr, 2, function(x) x > pop.medians)
nutr_summ_test <- as.data.frame(nutr_summ_test)
nutr_summ_test <- sapply(nutr_summ_test, as.numeric)  # this gives a matrix of 0s and 1s to use for jaccard distance.

# drop soylents
nutr_summ_test <- nutr_summ_test[,!colnames(nutr_summ_test) %in% c("MCTs11", "MCTs12")]


require(vegan)
jacc_summ_dists <- vegdist(t(nutr_summ_test), method = "jaccard")

require(ape)
jacc_summ_pcoa <- as.data.frame(pcoa(jacc_summ_dists)$vectors)    

require(ggplot2)

plot <- jacc_summ_pcoa[1:2]
plot$UserName <- rownames(plot)

ggplot(data = plot, aes(x = Axis.1, y = Axis.2, color = UserName)) + geom_point()
