# Make edge table for all dietary data from study
# should be long format, need to specify which are source attributes and which are target attributes
# and which are edge attributes when you load into cytoscape

require(reshape2)
require(tidyverse)
require(tibble)

setwd(dir = "/Users/abby/Documents/Projects/dietstudy")

map <- read.table(file = "data/maps/UserName_map.txt", sep = "\t", header = T)
food <- read.table(file ="data/processed_food/food.smry.txt", sep = "\t", header = T, comment = "", row.names = "taxonomy")
food_taxonomy <- read.table(file = "raw/diet.taxonomy.txt", sep = "\t", header = T)

### limit food to L6
food <- food[,!colnames(food) == "X.FoodID"]
split <- strsplit(rownames(food), ";")
foodStrings <- sapply(split, function(x) paste(x[1:6], collapse= ";"))
food <- rowsum(food, foodStrings)

### make edge table ###
food <- rownames_to_column(food, "L6")
food <- melt(food)
food <- food %>% select(variable, L6, value)
food$L1 <- food$L6

edge.table <- food
edge.table$L1 <- gsub(";L2_.*$", "",edge.table$L1)

colnames(edge.table) <- c("Source", "Target", "amount", "taxonomy")
edge.table <- edge.table %>% filter(!(amount == 0))

edge.table <- inner_join(edge.table, map, by = c("Source"="UserName"))
edge.table <- edge.table %>% select(Source, Target, amount, taxonomy, Gender, Age, Supplement, Activity.Factor, mean_kcal)

write.table(edge.table, "figures/figure 3/food_network/food.edge.table.txt", sep = "\t", quote = F, row.names = F, col.names = T)


# ### make node table ###
# this will make a variable that distinguishes between food_node and user_node
node.table <- map %>% select(UserName)
node.table$ntype <- "user_node"
colnames(node.table)[1] <- "node_name"
node.table <- node.table %>% filter(node_name %in% unique(edge.table$Source))

 
node.table2 <- NULL
node.table2$node_name <- unique(edge.table$Target)
node.table2 <- as.data.frame(node.table2)
node.table2$ntype <- "food_node"

node.table <- rbind(node.table,node.table2)
 
write.table(node.table, "figures/figure 3/food_network/food.node.table.txt", sep = "\t", quote = F, row.names = F, col.names = T)
