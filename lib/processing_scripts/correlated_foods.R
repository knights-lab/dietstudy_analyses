# Food v. food correlations within a person?
# What foods do people eat together as meals? Can we identify them from the data?


require(reshape2)
require(psych)
require(ggpubr)

setwd("/Users/abby/Documents/Projects/dietstudy/")

# load map
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID", "UserName","StudyDayNo","Diet.Week.Day")]

# load food table for everyone
food <- read.delim("data/processed_food/dhydrt.txt", row = "taxonomy")
food <- food[,!colnames(food) == "X.FoodID"]

# collapse at different levels?
split <- strsplit(rownames(food),";")             
foodStrings <- sapply(split,function(x) paste(x[1:3],collapse=";"))
food<- rowsum(food,foodStrings)


#limit map to the people with mb samples
map <- map[map$X.SampleID %in% colnames(food),]
map <- droplevels(map)

# listify distance matrices by person
list.food<- NULL

for (i in seq_along(unique(map$UserName))){
  a <- as.character(unique(map$UserName)[i])
  ids <- droplevels(map$X.SampleID[map$UserName == a])
  
  sub.food <- food[,colnames(food) %in% ids]
  list.food[[a]] <- sub.food
  
} 


# find repeated correlations across people
mycors.list <- NULL
for (i in seq_along(list.food)) {
  mycors <- cor(t(list.food[[i]]), method = "spearman")
  mycors[lower.tri(mycors, diag = T)] <- NA
  mycors <- melt(mycors, value.name = "correlation")
  
  mycors.p <- corr.test(t(list.food[[i]]), adjust = "fdr", use = "complete")$p
  mycors.p[lower.tri(mycors.p, diag = T)] <- NA
  mycors.p <- melt(mycors.p, value.name = "fdr_p")
  
  test <- merge(mycors, mycors.p)
  test <- test[!is.na(test$correlation),]
  test <- test[test$fdr_p < 0.01,]

  mycors.list[[i]] <- test
}


foods.cors <- unlist(lapply(mycors.list, function(x) dim(x)[[1]]))
names(foods.cors) <- names(list.food)
mean(foods.cors[foods.cors > 3])

ggqqplot(foods.cors) # approximately normal distributrion

mean(foods.cors)
sd(foods.cors)
