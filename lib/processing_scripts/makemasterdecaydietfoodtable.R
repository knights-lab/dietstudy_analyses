## GOAL:
# Would a weighted average of the dietary profiles with 
# decaying weights 2^(-n) for day 'n' do better or worse?

# Make one decaying day diet table for the entire study
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




# New goal is to sum days with decaying weights so for a person with diet days from days prior
# for MB day 3 we would sum Diet day 0, 1, and 2 with a rate of decay of 0.5 
# so: diet from day 2 + 1/2(diet day 1) + 1/4(diet day 2)
# then continue with subsequent days

masterdecay <- NULL

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

  # use rollapply to sum every three rows
  # limit the dataframe to just the summable columns
  fooddecay <- foodsub[,5:ncol(foodsub)]
  my_fun = function(x) {d = (length(x):1)-1; f = 2^-d; sum(f*x)}
  fooddecay <- rollapply(fooddecay, width = 17, FUN = my_fun, by.column=T, align="right", partial = T) 

  # add back information to map diet days to microbiome samples
  mynames <- submap$X.SampleID
  
  rownames(fooddecay) <- mynames
  
  tmp <- fooddecay
  
  # append each person to the same master dataframe
  masterdecay <- rbind(masterdecay,tmp)

  }


# now drop any with rowsums == NA to get the final dataframe

masterdecay <- masterdecay[!is.na(rowSums(masterdecay)),]

# for processing with qiime, this needs to be transposed and have the taxonomy added back in 

masterdecay <- as.data.frame(t(masterdecay))
masterdecay$taxonomy <- food$taxonomy
masterdecay <- rownames_to_column(masterdecay, var = "#SampleID")

write.table(masterdecay, file = "data/diet/processed_food/masterdecaydiet.txt", sep = "\t", quote = F, row.names = F, col.names = T)
