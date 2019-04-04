require(rmarkdown)
require(knitr)
require(ggplot2)
require(reshape2)

# set the path for root directory
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/")

# load map
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID", "UserName","StudyDayNo","Sample.Week.Day")]

# load mb table for everyone
mb <- read.delim("data/microbiome/processed_sample/taxonomy_clr_s.txt", row = 1)

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
day1diffall <- NULL
weekdaysall <- NULL
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
  day1diff[[i]] <- mean(unlist(vals[n+1]))
  day1diffall[[i]] <- unlist(vals[n+1])
  # get the mean of all the other distances
  allotherdaysdiff[[i]] <- mean(unlist(vals[(n+2):(2*n-1)]))
  
  # get the days for plotting
  days <- colnames(x)[2:length(colnames(x))]
  w <- subset(map, map$X.SampleID %in% days)
  weekdaysall[[i]] <- as.vector(w[,"StudyDayNo"])
  
}

x <- unlist(day1diff)
y <- unlist(allotherdaysdiff)



wilcox.test(x, y, paired = T)



#plot fig.height=2.5, fig.width=2.5


x <- as.data.frame(x)
colnames(x) <- "distance"
x$mb_dist_to <- "To next day"

y <- as.data.frame(y)
colnames(y) <- "distance"
y$mb_dist_to <- "To other days"

plot <- rbind(x,y)
myplot <- ggplot(data = plot, aes(x = mb_dist_to, y = distance)) + 
  geom_jitter(alpha = 0.7, width = 0.35, color = "#5f86b7", size = 2) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  annotate("text", x = 1.5, y = 25, label = "italic(p) < 0.001", size = 3, parse = T) +
  theme_classic() +
  xlab("Distance from today") +
  ylab("Mean Aitchison's distance")+
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size =7))
  
  
myplot

ggsave("output/Figure4/Figure4B.pdf", myplot, width = 2.5, height = 2.5)
