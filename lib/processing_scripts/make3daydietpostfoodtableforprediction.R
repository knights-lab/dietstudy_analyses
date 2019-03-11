# Make one complete 3 day diet table for the entire study
require(zoo)
require(tidyr)

setwd(dir = "~/Documents/Projects/dietstudy_analyses/")

# read in the Food totals for the complete dataset
food <- read.delim(file = "data/diet/processed_food/dhydrt.txt", row.names = 1)

# read in the mapping data and strip blanks and dropouts
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID","UserName", "StudyDayNo","DietDayNo", "Study.Status")]
map <- map[map$Study.Status %in% c("Complete", "No blood draw"),]
map <- map[colnames(map) %in% c("X.SampleID","UserName", "StudyDayNo","DietDayNo")]
map <- droplevels(map)

# Goal is to sum sets of three days that are decoupled from microbiome days
# want to label these with the first day's date so they can be matched to microbiome samples prior to diet intake
# triplets of food days correspond to microbiome days (1,2,3 = 1; 2,3,4 = 2; 3,4,5 = 3, etc.)

master3day <- NULL

for (i in 1:length(unique(map$UserName))) {
  x <- unique(map$UserName)[i]
  
  submap <- map[map$UserName == x,]
  
  StudyDayNo <- c(1:17)
  StudyDayNo <- as.data.frame(StudyDayNo)
  
  index <- merge(StudyDayNo, submap, all.x = T)
  
  foodsub <- food[,colnames(food) %in% submap$X.SampleID]
  foodsub <- as.data.frame(t(foodsub))
  foodsub$X.SampleID <- rownames(foodsub)
  
  foodsub <- merge(index, foodsub, all.x = T)  

  #double check sorting
  foodsub <- foodsub[order(foodsub$StudyDayNo),]
  
  # drop diet day 0, since that is prior to the first sample day
  foodsub <- foodsub[!foodsub$DietDayNo %in% c(0,1),]

  # use rollapply to sum every three rows
  # limit the dataframe to just the summable columns
  food3day <- foodsub[,5:ncol(foodsub)]
  food3day <- as.data.frame(rollapply(food3day, 3, FUN = sum)) # note, this works so any groupings with missing days are summed to NA

  # add back information to map diet days to microbiome samples
  mynames <- submap$X.SampleID[1:dim(food3day)[1]]
  
  rownames(food3day) <- mynames
  
  tmp <- food3day
  
  # append each person to the same master dataframe
  master3day <- rbind(master3day,tmp)

  }


# now drop any with rowsums == NA to get the final dataframe

master3day <- master3day[!is.na(rowSums(master3day)),]

# for processing with qiime, this needs to be transposed and have the taxonomy added back in 

master3day <- as.data.frame(t(master3day))
master3day$taxonomy <- food$taxonomy
master3day <- rownames_to_column(master3day, var = "#SampleID")

write.table(master3day, file = "data/diet/processed_food/post3daydietforprediction.txt", sep = "\t", quote = F, row.names = F, col.names = T)
