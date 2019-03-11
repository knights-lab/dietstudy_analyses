### Reviewer comment #1 ###

#Figure S2 monte carlo test sounds strange for what the authors want to test (page 4 top). 
#If the idea is to test whether the day before was the most similar, 
#why not just test it directly and compare with distance to other days, 
#rather than this permutation test? Am I missing something? 
#And why is Bray-Curtis distance used here?

### Reviewer comment #2 ###

#The authors state that "each subject's microbiome samples were most similar to 
#their own sample from the previous day compared to all other timepoints from that subject".
#Was this done also with comparing to all other unrelated samples, or only to samples within
#the same subject? The proper analysis should be to compare as background to the whole cohort.


### Notes from meeting with Dan ###

#so the analysis should be redone using the Aitchison's distance. 
#Get each day's distance to the next day within in a person and compare to all 
# of the distances that are 2 or more days apart.

#Then test whether the mean difference is <0 with a t-test. 

# load map
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID", "UserName","StudyDayNo","Sample.Week.Day")]

# load mb table for everyone
mb <- read.delim("data/processed_tax/taxonomy_clr_s.txt", row = 1)

# get mb distances (euclidean from clr)
mbdist <- dist(t(mb))
mbdist <- as.matrix(mbdist)

#limit map to the people with mb samples
map <- map[map$X.SampleID %in% colnames(mb),]
map <- droplevels(map)

# listify distance matrices by person
list.mbdist<- NULL

for (i in seq_along(unique(map$UserName))){
  a <- as.character(unique(map$UserName)[i])
  ids <- droplevels(map$X.SampleID[map$UserName == a])
  
  sub.mbdist <- mbdist[rownames(mbdist) %in% ids,colnames(mbdist) %in% ids]
  list.mbdist[[i]] <- sub.mbdist
  
} 

# calculate distance to the next day and distance to days 2 days apart

### This is the meat of the problem, and where I think I'm probably messing up! ###

day1diff <- NULL
allotherdaysdiff <- NULL
for (i in 1:length(list.mbdist)) {
  # subset to just one person 
  x <- list.mbdist[[i]]
  n <- dim(x)[1]
  # pull the distances between sucessive days
  # create an indicator for the diagonals in the matrix
  d <- row(x) - col(x)
  # use split to group on these values by their diagonals 
  vals <- split(x, d)
  # get the mean of the diagonal that corresponds to 1 day diff
  day1diff[[i]] <- unlist(vals[n+1])
  # get the mean of all the other distances
  allotherdaysdiff[[i]] <- unlist(vals[(n+2):(2*n-1)])
  
}

x <- unlist(day1diff)
y <- unlist(allotherdaysdiff)

# test if distance to next day is less than distance to all other days
t.test(x, y, alternative = "less")

#plot

x <- as.data.frame(x)
colnames(x) <- "distance"
x$mb_dist_to <- "To Next Day"

y <- as.data.frame(y)
colnames(y) <- "distance"
y$mb_dist_to <- "To Other Days"

plot <- rbind(x,y)
myplot <- ggplot(data = plot, aes(x = mb_dist_to, y = distance)) + 
  geom_jitter(alpha = 0.2, width = 0.35, pch = 21, fill = "#64baaa") +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  annotate("text", x = 1.5, y = 30, label = "italic(p) <= 0.001", size = 4, parse = T) +
  theme_classic() +
  xlab("Microbiome Distance from Today") +
  ylab("Aitchison’s distance")
  
  
myplot


