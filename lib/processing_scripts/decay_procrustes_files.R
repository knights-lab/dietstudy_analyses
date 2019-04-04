

require(rmarkdown)
require(knitr)
require(tidyverse)
require(stringr)
require(dplyr)
require(reshape2)
require(zoo)
require(RColorBrewer)


setwd("/Users/abby/Documents/Projects/dietstudy_analyses/")


map_sample <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
map_username <- read.table("data/maps/UserName_map.txt", sep = "\t", header = T, comment = "")

map_sample$StudyDate <- as.Date.factor(map_sample$StudyDate, format = "%m/%d/%y")

#load other maps
food_map <- read.table("data/maps/food_map_no_soy.txt", sep = "\t", header = T, comment = "")
tax_map <- read.table("data/maps/taxonomy_norm_map.txt", sep = "\t", header = T, comment = "")


tax <- read.delim("data/microbiome/processed_sample/taxonomy_clr_s.txt", row = 1)

food <- read.table("data/diet/processed_food/dhydrt.txt", header = T, sep = "\t", comment = "", row = 1)
food_taxa <- read.table("data/diet/diet.taxonomy.txt", header = T, sep = "\t", comment = "", row = 3)
food_taxa <- food_taxa[,colnames(food_taxa) != "FoodID", drop=F]

# drop samples from the food_map that are not in the food file and vice versa to drop soylents
food_map <- food_map[food_map$X.SampleID %in% colnames(food),] # drop rows
food <- food[,colnames(food) %in% food_map$X.SampleID]

table(food_map$UserName) # almost all people have complete diet data (only 3 incomplete timelines)

table(tax_map$UserName) # many more people are missing microbiome samples


# remove blanks
map_sample <- map_sample %>% filter(Study.Status != "Blank")

# remove dropouts
map_sample <- map_sample %>% filter(Study.Status != "Dropped")

# remove soylent
map_sample <- map_sample %>% filter(UserName != "MCTs11", UserName != "MCTs12")


# make a table for tax, map, and food for each individual
#loop through and write each to a txt file

#make lists to fill
foodsub <- list()
taxsub <- list()
mapsub <- list()

# calculate decaying influence diet using the decaying version of rollapply
for (i in unique(map_sample$UserName)){
  subgroup <- map_sample[map_sample$UserName == i,]
  
  # foods
  foodsub[[i]] <- food[,colnames(food) %in% subgroup$X.SampleID] # subset to one username
  foodsub[[i]] <- foodsub[[i]][,order(names(foodsub[[i]]))] # ensure order is sequential by day

  newnames <- colnames(foodsub[[i]]) 
  # apply rollapply
  foodsub[[i]] <- t(foodsub[[i]])                                   # transpose for the rollapply function
  my_fun = function(x) {d = (length(x):1)-1; f = 2^-d; sum(f*x)}
  foodsub[[i]] <- rollapply(foodsub[[i]], width = 17, FUN = my_fun, by.column=T, align="right", partial = T) 
  foodsub[[i]] <- t(foodsub[[i]])                                   # transpose back
  foodsub[[i]] <- as.data.frame(foodsub[[i]])
  # fix names on foodsub df
  #rename the new summarized data starting with the first summarized day's name
  
  colnames(foodsub[[i]]) <- newnames 
  
  # pair tax and summarized food days
  taxsub[[i]] <- tax[,colnames(tax) %in% colnames(foodsub[[i]]), drop = F]  # subset to match the foodsub ids
  foodsub[[i]] <- foodsub[[i]][,colnames(foodsub[[i]]) %in% colnames(taxsub[[i]])] # do the same to match the taxsub ids
  
  # prep for export
  taxsub[[i]] <- rownames_to_column(taxsub[[i]], var = "#taxonomy")
  foodsub[[i]] <- merge(foodsub[[i]], food_taxa, by = 0)         # add taxonomy information back to each df
  #rownames(foodsub[[i]]) <- foodsub[[i]][,"Row.names"]           # format for export
  #foodsub[[i]] <- foodsub[[i]][,colnames(foodsub[[i]]) != "Row.names", drop = F]  # format for export
  colnames(foodsub[[i]])[1] <- "#FoodID"                         # add name for export
  
  
  # map
  mapsub[[i]] <- subgroup[subgroup$X.SampleID %in% colnames(foodsub[[i]]),] # subset to match the foodsub ids
  colnames(mapsub[[i]])[1] <- "#SampleID"

  }


for (df in names(mapsub)) {
  #write foods to table
  write.table(foodsub[[df]], file = paste0("data/procrustes/data_username_decay/", df, "_food", ".txt"), sep = "\t", quote = F,  row.names = F, col.names = T)
  
  # write taxa to table
  write.table(taxsub[[df]], file = paste0("data/procrustes/data_username_decay/", df, "_tax", ".txt"), sep = "\t", quote = F,  row.names = F, col.names = T)
  
  # write map to table
  write.table(mapsub[[df]], file = paste0("data/procrustes/data_username_decay/", df, "_map", ".txt"), sep = "\t", quote = F,  row.names = F, col.names = T)
  
}



# scripts for qiime are in lib/procrustes.sh

# or in lib/procrustes_offset.sh





