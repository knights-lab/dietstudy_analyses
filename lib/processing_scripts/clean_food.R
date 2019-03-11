require(dplyr)
require(tibble)

homedir <- "/Users/abby/Documents/Projects/dietstudy/"
setwd(dir = homedir)

# load map
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# load food taxonomy
food_tax <- read.table("raw/diet.taxonomy.txt", sep = "\t", header = T, comment = "")
# drop foodID
food_tax <- food_tax[,!colnames(food_tax) == "FoodID"]
# rename main.food.description to be "#FoodID" for later
colnames(food_tax) <- c("taxonomy", "#FoodID")

# load food otus and nutrition data
dhydrt <- read.table("raw/diet.dhydrt.txt", header = T, sep = "\t", comment = "")
food <- read.table("raw/diet.food.txt", header = T, sep = "\t", comment = "") 
fiber <- read.table("raw/diet.fiber.txt", header = T, sep ="\t", comment = "")
nutr <- read.table("raw/nutrition_totals.txt", header = T, sep = "\t", comment = "")
food_items <- read.delim("raw/Uncleaned_diet/MCTs_23887_Items.csv", header = T, sep = ",", comment = "")

# determine dietary outliers from nutritional data table
# ASA24 recommends a few differnet methods for identifying outliers, 
# Based on the 5th and 95th percentile reported in NHANES:
nutrmap <- inner_join(nutr, map, by = "X.SampleID")

# we excluded dietary data for days when the particpant reported marconutrient intakes as follows:

# Kcal:
# Women <600 and >4400
# Men <650 and >5700
# Protein: 
# Women <10 and >180
# Men <25 and >240
# Fat:
# Women <15 >185
# Men <25 >230

macronutrient_outliers <- nutrmap %>% 
  filter((Gender == "F" & (KCAL < 600 | KCAL > 4400)) | 
           (Gender == "M" & (KCAL < 650 | KCAL > 5700)) |
           (Gender == "F" & (PROT < 10 | PROT > 180)) | 
           (Gender == "M" & (PROT < 25 | PROT > 240)) |
           (Gender == "F" & (TFAT < 15 | TFAT > 185)) | 
           (Gender == "M" & (TFAT < 25 | TFAT > 230)))

# 11 observations
# dietary outliers
outliers = data.frame(macronutrient_outliers$UserName.x, macronutrient_outliers$RecordDayNo)
colnames(outliers) <- c("UserName", "RecordDayNo")

# manual review of outliers using items file
outlier_items <- merge(food_items, outliers) 

# export for review in excell
write.table(outlier_items, file = "raw/outlier_items_for_review.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# after manual review, only one record from MCT02 (day 4) and MCTO5(day 8) appeared to be incomplete.
# others appeared to simply suggest normal fluctuations in dietary habits and seem complete.
# Participant MCT05 also provided partial intake data for MCT.f.0085, so that sample is also removed.

outliers_to_drop <- macronutrient_outliers[c(1,3,7),]
outliers_to_drop <- outliers_to_drop$X.SampleID
outliers_to_drop <- c(outliers_to_drop, "MCT.f.0085")

# drop dropouts 
dropouts =  as.character(map[map$UserName %in% c("MCTs02", "MCTs17", "MCTs30"),"X.SampleID"])

# drop outliers and dropouts from each food table
dhydrt <- dhydrt[,!(colnames(dhydrt) %in% c(outliers_to_drop, dropouts))]
food <- food[,!(colnames(food) %in% c(outliers_to_drop, dropouts))]
fiber <- fiber[,!(colnames(fiber) %in% c(outliers_to_drop, dropouts))]
nutr <- nutr[!nutr$X.SampleID %in% c(outliers_to_drop, dropouts),]


# normalize and round dhydrate for QIIME test
dhydrt.n <- dhydrt
l.cols <- ncol(dhydrt.n)-1
dhydrt.n[,2:l.cols] <- sweep(dhydrt.n[,2:l.cols], 2, colSums(dhydrt.n[,2:l.cols]), "/")

# multiply normalized table by median food weight
m.cols <- median(colSums(dhydrt[,2:l.cols]))
dhydrt.n[,2:l.cols] <- ceiling(dhydrt.n[,2:l.cols] * m.cols) # this makes a nice normalized counts table for qiime


# normalize and round fiber for QIIME test
fiber.n <- fiber
l.cols <- ncol(fiber.n)-1
fiber.n[,2:l.cols] <- sweep(fiber.n[,2:l.cols], 2, colSums(fiber.n[,2:l.cols]), "/")

# multiply normalized table by median food weight
m.cols <- median(colSums(fiber[,2:l.cols]))
fiber.n[,2:l.cols] <- ceiling(fiber.n[,2:l.cols] * m.cols) # this makes a nice normalized counts table for qiime


# fix nameing of first column to use with QIIME
colnames(food)[1] <- "#FoodID"
colnames(dhydrt)[1] <- "#FoodID"
colnames(fiber)[1] <- "#FoodID"
colnames(dhydrt.n)[1] <- "#FoodID"
colnames(fiber.n)[1] <- "#FoodID"

# make a map that matches the food tables (for QIIME)
food_map <- map[map$X.SampleID %in% colnames(food),] #568 food records + FoodID and taxonomy
colnames(food_map)[1] <- "#SampleID"

# filter to most commonly consumed foods (similar to filtering for taxonomy)
# Drop foods that don't appear in at least 
food_per_person <- list()

for (i in unique(food_map$UserName)) {
  submap <- food_map[food_map$UserName == i,]
  subset <- dhydrt[, colnames(dhydrt) %in% submap$`#SampleID`]
  subset <- subset[rowSums(subset > 0) > ncol(subset)/12, ] # these are the foods this person ate on at least 1/12 of the study days
  myfoods <- rownames(subset)
  myfoods <- dhydrt[rownames(dhydrt) %in% myfoods,]
  food_per_person[[i]] <- myfoods$taxonomy
}

n25 <- round(length(food_per_person) * 0.1) # 10% of people
nonuni <- unlist(food_per_person)
counts <- table(nonuni)
counts <- counts[counts > n25]
keep <- as.vector(names(counts))

# limit to foods we know many people ate during the study
dhydrt_limited <- dhydrt[dhydrt$taxonomy %in% keep,]  
fiber_limited <- fiber[fiber$taxonomy %in% keep,]

# make no soylent version
soylent_ids <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
dhydrt_no_soy <- dhydrt[!colnames(dhydrt) %in% soylent_ids$X.SampleID]
fiber.n_no_soy <- fiber.n[!colnames(fiber.n) %in% soylent_ids$X.SampleID]
food_map_no_soy <- food_map[!food_map$`#SampleID` %in% soylent_ids$X.SampleID,]

# make pre and post versions
pre_ids <- sort(na.omit(c(as.character(map[map$Timing == "Pre",]$X.SampleID), "#FoodID", "taxonomy")))
post_ids <- sort(na.omit(c(as.character(map[map$Timing == "Post",]$X.SampleID), "#FoodID", "taxonomy")))
dhydrt_pre <- dhydrt[colnames(dhydrt) %in% pre_ids]
dhydrt_post <- dhydrt[colnames(dhydrt) %in% post_ids]
fiber_pre <- fiber[colnames(fiber) %in% pre_ids]
fiber_post <- fiber[colnames(fiber) %in% post_ids]


# make base nutritional variables
nutr_65 <- nutr %>% select(-c(UserName, StudyDayNo, RecordDayNo))
nutr_65 <- remove_rownames(nutr_65)
nutr_65 <- column_to_rownames(nutr_65, var = "X.SampleID")
nutr_65 <- as.data.frame(t(nutr_65))
nutr_65 <- nutr_65[1:65,]
nutr_65 <- rownames_to_column(nutr_65, var = "#SampleID")


# change most variables in the nutr_65 table into grams (i.e. mg to grams, etc)
# Prot, fat, carb, mois, alc, sugars,fiber, sfat, mfat, pfat
# and all the specific fatty acids are already in grams

# mg variables
mgvar <- c("CAFF", "THEO", "CALC", "IRON", "MAGN", "PHOS", "POTA", "SODI",
           "ZINC", "COPP", "VC", "VB1", "VB2", "NIAC", "VB6", "ATOC", "CHOLE",
           # "SFAT", "S040", "S060", "S080", "S100", "S120", "S140", "S160", "S180",
           # "MFAT", "M161", "M181", "M201", "M221", 
           # "PFAT", "P182", "P183", "P184", "P204", "P205", "P225", "P226",
           "CHOLN", "VITE_ADD")
mcgvar <- c("SELE", "FOLA", "FA", "FF", "FDFE", "VB12", "VARA", "RET", "BCAR", 
            "ACAR", "CRYP", "LYCO", "LZ", "VK", "VITD", "B12_ADD")

#mg values to grams
mgrows <- which(nutr_65[,1] %in% c(mgvar))
nutr_65[mgrows,2:ncol(nutr_65)] <- nutr_65[mgrows,2:ncol(nutr_65)]/1000

#mcg to grams 
mcgrows <- which(nutr_65[,1] %in% c(mcgvar))
nutr_65[mcgrows,2:ncol(nutr_65)] <- nutr_65[mcgrows,2:ncol(nutr_65)]/1000000

# make hei nutrition table
nutr_hei <- nutr %>% select(-c(UserName, StudyDayNo, RecordDayNo))
nutr_hei <- remove_rownames(nutr_hei)
nutr_hei <- column_to_rownames(nutr_hei, var = "X.SampleID")
nutr_hei <- as.data.frame(t(nutr_hei))
nutr_hei <- nutr_hei[66:(nrow(nutr_hei)-1),]
nutr_hei <- rownames_to_column(nutr_hei, var = "#SampleID")

# write tables
write.table(food, file = "data/processed_food/food.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt, file = "data/processed_food/dhydrt.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_limited, file = "data/processed_food/dhydrt.limited.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_no_soy, file = "data/processed_food/dhydrt.no.soy.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt.n, file = "data/processed_food/dhydrt.norm.counts.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_pre, file = "data/processed_food/dhydrt.pre.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_post, file = "data/processed_food/dhydrt.post.txt", quote = F, sep = "\t", row.names = F, col.names = T)

write.table(fiber, file = "data/processed_food/fiber.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(fiber.n, file ="data/processed_food/fiber.norm.counts.txt", quote = F, sep = "\t",row.names = F, col.names = T)
write.table(fiber.n_no_soy, file ="data/processed_food/fiber.norm.counts.no.soy.txt", quote = F, sep = "\t",row.names = F, col.names = T)
write.table(fiber_limited, file = "data/processed_food/fiber.limited.txt", sep = "\t", row.names = F, col.names = T)
write.table(fiber_pre, file = "data/processed_food/fiber.pre.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(fiber_post, file = "data/processed_food/fiber.post.txt", quote = F, sep = "\t", row.names = F, col.names = T)

write.table(nutr, file = "data/processed_nutr/nutr_totals.txt", quote = F, sep ="\t", row.names = F, col.names = T)
write.table(nutr_65, file = "data/processed_nutr/nutr_65.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(nutr_hei, file = "data/processed_nutr/nutr_hei.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(food_map, file="data/maps/food_map.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(food_map_no_soy, file = "data/maps/food_map_no_soy.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# summarize by person

# set up the summary dataframes
food_smry <- data.frame(matrix(nrow=nrow(food), ncol=0))
dhydrt_smry <- data.frame(matrix(nrow=nrow(dhydrt), ncol=0))
dhydrt_limited_smry <- data.frame(matrix(nrow=nrow(dhydrt_limited), ncol=0))
dhydrt_pre_smry <- data.frame(matrix(nrow=nrow(dhydrt), ncol=0))
dhydrt_post_smry <- data.frame(matrix(nrow=nrow(dhydrt), ncol=0))
fiber_smry <- data.frame(matrix(nrow=nrow(fiber), ncol=0))
fiber_limited_smry <- data.frame(matrix(nrow=nrow(fiber_limited), ncol=0))
nutr_smry <- data.frame(matrix(nrow=nrow(nutr_65), ncol=0))
nutr_hei_smry <- data.frame(matrix(nrow=nrow(nutr_hei), ncol=0))

# add rownames to each
rownames(food_smry) <- food$taxonomy
rownames(dhydrt_smry) <- dhydrt$taxonomy
rownames(dhydrt_limited_smry) <- dhydrt_limited$taxonomy
rownames(dhydrt_pre_smry) <- dhydrt$taxonomy
rownames(dhydrt_post_smry) <- dhydrt$taxonomy
rownames(fiber_smry) <- fiber$taxonomy
rownames(fiber_limited_smry) <- fiber_limited$taxonomy
rownames(nutr_smry) <- nutr_65$`#SampleID`
rownames(nutr_hei_smry) <- nutr_hei$`#SampleID`

#loop through and summarize by username
for (i in unique(food_map$UserName)){
  subgroup <- food_map[food_map$UserName == i,]
  # food
  foodsub <- food[,colnames(food) %in% subgroup$`#SampleID`]
  tmp <- as.data.frame(rowMeans(foodsub))                    # need to choose if means or sums is better for analysis
  colnames(tmp) <- i
  food_smry <- cbind(food_smry,tmp)
  # dhydrt
  dhydrtsub <- dhydrt[,colnames(dhydrt) %in% subgroup$`#SampleID`]
  tmp2 <- as.data.frame(rowMeans(dhydrtsub))                    # need to choose if means or sums is better for analysis
  colnames(tmp2) <- i
  dhydrt_smry <- cbind(dhydrt_smry,tmp2)
  # dhydrt_limited
  dhydrt_limitedsub <- dhydrt_limited[,colnames(dhydrt_limited) %in% subgroup$`#SampleID`]
  tmp3 <- as.data.frame(rowMeans(dhydrt_limitedsub))                    # need to choose if means or sums is better for analysis
  colnames(tmp3) <- i
  dhydrt_limited_smry <- cbind(dhydrt_limited_smry,tmp3)
  #dhydrt_pre
  dhydrt_presub <- dhydrt_pre[,colnames(dhydrt_pre) %in% subgroup$`#SampleID`]
  tmp8 <- as.data.frame(rowMeans(dhydrt_presub))                    # need to choose if means or sums is better for analysis
  colnames(tmp8) <- i
  dhydrt_pre_smry <- cbind(dhydrt_pre_smry,tmp8)
  #dhydrt_post
  dhydrt_postsub <- dhydrt_post[,colnames(dhydrt_post) %in% subgroup$`#SampleID`]
  tmp9 <- as.data.frame(rowMeans(dhydrt_postsub))                    # need to choose if means or sums is better for analysis
  colnames(tmp9) <- i
  dhydrt_post_smry <- cbind(dhydrt_post_smry,tmp9)
  # fiber
  fibersub <- fiber[,colnames(fiber) %in% subgroup$`#SampleID`]
  tmp4 <- as.data.frame(rowMeans(fibersub))                    # need to choose if means or sums is better for analysis
  colnames(tmp4) <- i
  fiber_smry <- cbind(fiber_smry,tmp4)
  # fiber_limited
  fiber_limtedsub <- fiber_limited[,colnames(fiber_limited) %in% subgroup$`#SampleID`]
  tmp5 <- as.data.frame(rowMeans(fiber_limtedsub))
  colnames(tmp5) <- i
  fiber_limited_smry <- cbind(fiber_limited_smry,tmp5)
  # nutr summary
  nutr_smry_sub <- nutr_65[,colnames(nutr_65) %in% subgroup$`#SampleID`]
  tmp6 <- as.data.frame(rowMeans(nutr_smry_sub))
  colnames(tmp6) <- i
  nutr_smry <- cbind(nutr_smry, tmp6)
  # nutr healthy eating index summary
  nutr_hei_sub <- nutr_hei[,colnames(nutr_hei) %in% subgroup$`#SampleID`]
  tmp7 <- as.data.frame(rowMeans(nutr_hei_sub))
  colnames(tmp7) <- i
  nutr_hei_smry <- cbind(nutr_hei_smry, tmp7)
  
}

# move rownames to column and label as taxonomy
food_smry <- rownames_to_column(food_smry, "taxonomy")
dhydrt_smry <- rownames_to_column(dhydrt_smry, "taxonomy")
dhydrt_limited_smry <- rownames_to_column(dhydrt_limited_smry, "taxonomy")
dhydrt_pre_smry <- rownames_to_column(dhydrt_pre_smry, "taxonomy")
dhydrt_post_smry <- rownames_to_column(dhydrt_post_smry, "taxonomy")
fiber_smry <- rownames_to_column(fiber_smry, "taxonomy")
fiber_limited_smry <- rownames_to_column(fiber_limited_smry, "taxonomy")
nutr_smry <- rownames_to_column(nutr_smry, "nutr_value")
nutr_hei_smry <- rownames_to_column(nutr_hei_smry, "nutr_value")


# re-add the food id and move to the first column
food_smry <- merge(food_tax, food_smry)
dhydrt_smry <- merge(food_tax, dhydrt_smry)
dhydrt_limited_smry <- merge(food_tax, dhydrt_limited_smry)
dhydrt_pre_smry <- merge(food_tax, dhydrt_pre_smry)
dhydrt_post_smry <- merge(food_tax, dhydrt_post_smry)
fiber_smry <- merge(food_tax, fiber_smry)
fiber_limited_smry <- merge(food_tax, fiber_limited_smry)

# change order for saving
food_smry <- food_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)
dhydrt_smry <- dhydrt_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)
dhydrt_limited_smry <- dhydrt_limited_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)
dhydrt_pre_smry <- dhydrt_pre_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)
dhydrt_post_smry <- dhydrt_post_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)
fiber_smry <- fiber_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)
fiber_limited_smry <- fiber_limited_smry %>% select(-taxonomy,`#FoodID`, everything(), taxonomy)

# make no soylent version
dhydrt_smry_no_soy <- dhydrt_smry[,!colnames(dhydrt_smry) %in% c("MCTs11", "MCTs12")]
dhydrt_limited_smry_no_soy <- dhydrt_limited_smry[,!colnames(dhydrt_limited_smry) %in% c("MCTs11", "MCTs12")]
nutr_smry_no_soy <- nutr_smry[,!colnames(nutr_smry) %in% c("MCTs11", "MCTs12")]
nutr_hei_smry_no_soy <- nutr_hei_smry[,!colnames(nutr_hei_smry) %in% c("MCTs11", "MCTs12")]

# write summary tables
write.table(food_smry, file = "data/processed_food/food.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_smry, file = "data/processed_food/dhydrt.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_smry_no_soy, file = "data/processed_food/dhydrt.smry.no.soy.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_limited_smry, file = "data/processed_food/dhydrt.limited.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_limited_smry_no_soy, file = "data/processed_food/dhydrt.limited.smry.no.soy.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_pre_smry, file = "data/processed_food/dhydrt.pre.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(dhydrt_post_smry, file = "data/processed_food/dhydrt.post.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(fiber_smry, file = "data/processed_food/fiber.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(fiber_limited_smry, file = "data/processed_food/fiber.limited.smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(nutr_smry, file = "data/processed_nutr/nutr_65_smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(nutr_smry_no_soy, file = "data/processed_nutr/nutr_65_smry_no_soy.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(nutr_hei_smry, file = "data/processed_nutr/nutr_hei_smry.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(nutr_hei_smry_no_soy, file = "data/processed_nutr/nutr_hei_smry_no_soy.txt", quote = F, sep = "\t", row.names = F, col.names = T)

