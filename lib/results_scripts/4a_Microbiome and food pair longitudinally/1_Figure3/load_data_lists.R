


# load microbiome and diet flies as list for procrustes
###### Scale up for each person #########
# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_1day/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.1day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.1day <- lapply(mb.dist.1day, function(x) as.data.frame(pcoa(x)$vectors))

names(mb.pcoa.1day) <- gsub("_.*", "", temp)

##############################

# read in each dietary distance matrix and calculate pcoas for each person
# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.1day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.1day <- lapply(diet.dist.1day, function(x) as.data.frame(pcoa(x)$vectors))

names(diet.pcoa.1day) <- gsub("_.*", "", temp)

############################


# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_2day/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.2day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.2day <- lapply(mb.dist.2day, function(x) as.data.frame(pcoa(x)$vectors))

names(mb.pcoa.2day) <- gsub("_.*", "", temp)

##############################

# read in each dietary distance matrix and calculate pcoas for each person
# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.2day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.2day <- lapply(diet.dist.2day, function(x) as.data.frame(pcoa(x)$vectors))

names(diet.pcoa.2day) <- gsub("_.*", "", temp)

############################


# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_3day/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.3day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.3day <- lapply(mb.dist.3day, function(x) as.data.frame(pcoa(x)$vectors))

names(mb.pcoa.3day) <- gsub("_.*", "", temp)

##############################

# read in each dietary distance matrix and calculate pcoas for each person
# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.3day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.3day <- lapply(diet.dist.3day, function(x) as.data.frame(pcoa(x)$vectors))

names(diet.pcoa.3day) <- gsub("_.*", "", temp)

############################




############################


# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_4day/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.4day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.4day <- lapply(mb.dist.4day, function(x) as.data.frame(pcoa(x)$vectors))

names(mb.pcoa.4day) <- gsub("_.*", "", temp)

##############################

# read in each dietary distance matrix and calculate pcoas for each person
# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.4day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.4day <- lapply(diet.dist.4day, function(x) as.data.frame(pcoa(x)$vectors))

names(diet.pcoa.4day) <- gsub("_.*", "", temp)

############################




############################


# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_5day/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.5day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.5day <- lapply(mb.dist.5day, function(x) as.data.frame(pcoa(x)$vectors))

names(mb.pcoa.5day) <- gsub("_.*", "", temp)

##############################

# read in each dietary distance matrix and calculate pcoas for each person
# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.5day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.5day <- lapply(diet.dist.5day, function(x) as.data.frame(pcoa(x)$vectors))

names(diet.pcoa.5day) <- gsub("_.*", "", temp)

############################




# load microbiome and diet flies as list for procrustes
###### Scale up for each person #########
# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_decay/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# drop MCTs05, MCTs06, MCTs28 and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
mb.dist.decay <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.decay <- lapply(mb.dist.1day, function(x) as.data.frame(pcoa(x)$vectors))

names(mb.pcoa.decay) <- gsub("_.*", "", temp)

##############################

# read in each dietary distance matrix and calculate pcoas for each person
# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]

# drop MCTs05, MCTs06, MCTs28, and MCTs29 (less than 10 days of microbiome or food data points)
temp <- temp[grep("MCTs05", temp, invert = T)]
temp <- temp[grep("MCTs06", temp, invert = T)]
temp <- temp[grep("MCTs28", temp, invert = T)]
temp <- temp[grep("MCTs29", temp, invert = T)]

# create a list containing each of these
diet.dist.decay <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.decay <- lapply(diet.dist.decay, function(x) as.data.frame(pcoa(x)$vectors))

names(diet.pcoa.decay) <- gsub("_.*", "", temp)

############################


# clean up and rm the dist matrices, don't need them for downstream analysis
rm(list=ls(pattern = "dist"))
