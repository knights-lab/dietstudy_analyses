
setwd(dir = "/Users/abby/Documents/Projects/dietstudy_analyses/")
require(ape)
require(tidyverse)

# load map
map <- read.delim(file = "data/maps/SampleID_map.txt")
map <- map[colnames(map) %in% c("X.SampleID", "StudyDayNo","Sample.Week.Day","UserName")]

###### Scale up for each person #########
# read in each microbiome distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_decay/")

# get the file paths for the euclidean distance matrices for microbiome
temp <- list.files(pattern = "*_tax.txt", recursive = T)
temp <- temp[grep("euclidean", temp)]

# create a list containing each of these
mb.dist.1day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
mb.pcoa.1day <- lapply(mb.dist.1day, function(x) as.data.frame(pcoa(x)$vectors))

# identify people with < 5 time points
drops <- which(lapply(mb.pcoa.1day, function(x) dim(x)[1])<=6)

# drop the people with too few time points
mb.pcoa.1day <- mb.pcoa.1day[-drops]

# limit to the top 10 axis
mb.pcoa.1day <- lapply(mb.pcoa.1day, function (x) x[,1:5])
mb.pcoa.1day <- lapply(mb.pcoa.1day, setNames, paste0("Mb.Axis.", 1:5))


##############################

# read in each dietary distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_decay/")

# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_food.txt", recursive = T)
temp <- temp[grep("unweighted", temp)]


# create a list containing each of these
diet.dist.3day <- lapply(temp, function(x) read.delim(x, row.names = 1))

# make the PCOA from the distances
diet.pcoa.3day <- lapply(diet.dist.3day, function(x) as.data.frame(pcoa(x)$vectors))

# identify people with < 5 time points
drops <- which(lapply(diet.pcoa.3day, function(x) dim(x)[1])<=6)

# drop the people with too few time points
diet.pcoa.3day <- diet.pcoa.3day[-drops]


# limit to the top 5 axis
diet.pcoa.3day <- lapply(diet.pcoa.3day, function (x) x[,1:5])
diet.pcoa.3day <- lapply(diet.pcoa.3day, setNames, paste0("Dt.Axis.", 1:5))


# join the microbiome and diet dataframes and only keep days we have data for both
mylength <- length(mb.pcoa.1day)
predictors <- lapply(1:mylength, function(i) merge(mb.pcoa.1day[[i]], diet.pcoa.3day[[i]], by = 0))

############################

# load the microbiome relative abundance table for each person
# read in each dietary distance matrix and calculate pcoas for each person
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_decay/")

# get the file paths for the unweighted distance matrices for diet
temp <- list.files(pattern = "*_tax.txt")


# create a list containing each of these taxonomy tables
mb.ra <- lapply(temp, function(x) read.delim(x, row.names = 1))

# identify people with < 5 time points
drops <- which(lapply(mb.ra, function(x) dim(x)[2])<=6)

# drop the people with too few time points
mb.ra <- mb.ra[-drops]

# fix the naming of each species
mb.ra <- lapply(mb.ra, function(x) {
  rownames(x) <- gsub("?.*s__", "", rownames(x))
  rownames(x) <- gsub("?.*g__", "Uncl. Genus ", rownames(x))
  rownames(x) <- gsub("?.*f__", "Uncl. Family ", rownames(x))
  rownames(x) <- gsub("?.*o__", "Uncl. Order ", rownames(x))
  rownames(x) <- gsub("?.*p__", "Uncl. Phylum ", rownames(x))
  rownames(x) <- gsub("?.*k__", "Uncl. Kingdom ", rownames(x))
  rownames(x) <- gsub(";NA", "", rownames(x))
  rownames(x) <- gsub("\\[", "", rownames(x))
  rownames(x) <- gsub("\\]", "", rownames(x))
  rownames(x) <- gsub(" ", "_", rownames(x))
  rownames(x) <- gsub("-", "_", rownames(x))
  rownames(x) <- gsub("\\/", "_", rownames(x))
  return(x)
})

# sort each dataframe within the list by most abundant speices
mb.abund <- lapply(mb.ra, function(x) rowMeans(x))

# get names of top 100 most variable per person
top.abund <- lapply(mb.abund, function(x) names(head(sort(x, decreasing = T), n = 130)))

# # limit to just the top most abundant species
keepspecies <- names(which(sort(table(unlist(top.abund))) > 29))


# subset mb.ra to these
mb.ra <- lapply(1:mylength, function(i) mb.ra[[i]][keepspecies,])


#################################

# subset to just the samples for this person
mb.ra <- lapply(1:mylength, function(i) mb.ra[[i]][,colnames(mb.ra[[i]]) %in% predictors[[i]]$Row.names])


################################
# we care about the previous day ra for the model, so add these to the predictors
ra.predictors <- lapply(mb.ra, function(x) as.data.frame(t(x)))
ra.predictors <- lapply(ra.predictors, function(x) {
  x$Row.names <- rownames(x)
  return(x)
})

# we also care about these as the response variables
response <- lapply(mb.ra, function(x) as.data.frame(t(x)))
response <- lapply(response, function(x) {
  x$Row.names <- rownames(x)
  return(x)
})


# add the day variable to both predictors and response variables
predictors_day <- lapply(1:mylength, function(i) merge(map, predictors[[i]], by.x = "X.SampleID", by.y = "Row.names", all =F))
response_day <- lapply(1:mylength, function(i) merge(map, response[[i]], by.x = "X.SampleID", by.y = "Row.names", all = F))

# get the indicies to use for filtering
possible.days <- c(1:17) # vector to compare to 1:17 becuase there are 17 possible days
pred_days <- lapply(1:mylength, function(i) possible.days %in% predictors_day[[i]]$StudyDayNo) # this gives the list of days we values for

# make a list of the dataframes made from the vector
possible.days.df <- as.data.frame(possible.days, drop = F)

possible.days.list <- lapply(1:mylength, function(i) {cbind(possible.days.df, pred_days[[i]])})
possible.days.list <- lapply(possible.days.list, function(x){
  colnames(x) <- c("possible.days", "pred.days")
  return(x)
})
possible.days.list <- lapply(possible.days.list, function(x){
  x$pred_days <- ifelse(x$pred.days == TRUE, x$possible.days, NA)
  x$resp_days <- x$pred_days + c(diff(x$pred_days), NA)
  x$pred_days <- x$resp_days - 1
  keeps <- c("possible.days", "pred_days", "resp_days")
  x <- x[keeps]
  return(x)
})

# possible.days.list can now act as a framework to add the appropriate prediction and response days

# merge to get all the information in one table that we will use for testing
# either boosted regression or randomforests
dat <- lapply(1:mylength, function(i){merge(possible.days.list[[i]], predictors_day[[i]], by.x = "pred_days", by.y = "StudyDayNo")})
dat <- lapply(1:mylength, function(i){merge(dat[[i]], ra.predictors[[i]], by.x = "X.SampleID", by.y = "Row.names")})
dat <- lapply(1:mylength, function(i){merge(dat[[i]], response_day[[i]], by.x = "resp_days", by.y = "StudyDayNo")})


mynames <- as.character(temp)
mynames <- gsub("_.*", "", mynames)

# drop
mynames <- mynames[-drops]


dat <- setNames(dat, mynames)


# identify people with < 5 time points
drops <- which(lapply(dat, function(x) dim(x)[1])<=5)

# drop the people with too few time points
dat <- dat[!names(dat) %in% names(drops)]

save(dat, file = "../../../data/test_personal_diet_dat.Rdata")
save(keepspecies, file = "../../species.to.test.Rdata")
# ################### TEST DIET IMPACT ###########
# 
# 
# 
# # want to look at only same speices across people
# # so find the species in everyonen - can probably do this more eligantly later if needed 
# all_species <- unlist(most.abund)
# common_species <- names(which(table(all_species)==28)) # uses 38 common species found in everyone
# 
# 
# "test.diet.impact" <- function(){
#   
#   set.seed(42)
#   
#   Master.own.diet <- NULL
#   Master.others.diet <- NULL
#   Master.own.microbiome <- NULL
#   
#   for (k in 1:length(dat)) {
#     # person A index
#     a <- k
#     
#     # person B index
#     possible.b <- 1:length(dat)
#     possible.b <- possible.b[possible.b != a]
#     b <- sample(possible.b,1)
#     
#     # get the species names
#     today.species <- paste0(common_species, ".x")
#     tomorrow.species <- paste0(common_species, ".y")
#     
#     # create empty list for storing values
#     own.diet.predictions <- NULL
#     others.diet.predictions <- NULL
#     own.microbiome.predictions <- NULL
#     
#     
#     for (i in 1:length(tomorrow.species)) {
#       
#       tomorrow <- tomorrow.species[i]
#       today <- today.species[i]
#       
#       mb.features <- c(today, "Mb.Axis.1","Mb.Axis.2")
#       dt.features <- c("Dt.Axis.1","Dt.Axis.2")
#       
#       # predict person A's microbiome from microbiome alone
#       # regress tomorrow (yA) from today's species and microbiome
#       # predict tomorrow (yAhat)
#       yA.feats <- dat[[a]][,mb.features]
#       yA <- dat[[a]][,tomorrow]
#       yA.model <- lm(yA ~., yA.feats)
#       yAhat <- predict(yA.model)
#       
#       # get residuals
#       yAr <- yAhat - yA
#       
#       # regress the residuals on diet
#       modelA.feats <- dat[[a]][,dt.features]
#       modelA <- lm(yAr ~., modelA.feats)
#       yrAhat <- predict(modelA)
#       
#       # subtract predicted residuals (modelA) from yhat
#       # and correlate with y
#       own.diet <- cor(yAhat - yrAhat, yA)
#      
#       ##
#       # predict person B
#       # predict person B's species from B microbiome alone
#       yB.feats <- dat[[b]][,mb.features]
#       yB <- dat[[b]][,tomorrow]
#       yB.model <- lm(yB ~., yB.feats)
#       yBhat <- predict(yB.model)
#       # get residuals
#       yBr <- yBhat - yB
#       
#       #run predict on the model with person A, but passing in the diet data for person B 
#       #(you can pass in a new data matrix to predict) 
#       #to get their externally-predicted residuals yr.B.hat.
#       modelB.feats <- dat[[b]][,dt.features]
#       yrBhat <- predict(modelA, modelB.feats)
# 
#       #subtract yr.B.hat from yhat.B and correlate the final result with y.B.
#       others.diet <- cor(yBhat - yrBhat, yB)
#       
#       # compare yA and yAhat to 
#       own.microbiome <- cor(yAhat, yA)
#       
#       # store variables
#       
#       own.diet.predictions[[tomorrow]] <- own.diet
#       others.diet.predictions[[tomorrow]] <- others.diet
#       own.microbiome.predictions[[tomorrow]] <- own.microbiome
# 
#     }
#     Master.own.diet[[k]] <- own.diet.predictions
#     Master.others.diet[[k]] <- others.diet.predictions
#     Master.own.microbiome[[k]] <- own.microbiome.predictions
# 
#   }
#   
#   return(list(Master.own.diet = Master.own.diet, Master.others.diet = Master.others.diet, Master.own.microbiome = Master.own.microbiome))
# }
# 
# 
# myresults <- test.diet.impact()
# 
# # quick look at the means
# mean.own.diet <- unlist(lapply(myresults$Master.own.diet, function(x) mean(x)))
# mean.others.diet <- unlist(lapply(myresults$Master.others.diet, function(x) mean(x)))
# mean.own.microbiome <- unlist(lapply(myresults$Master.own.microbiome, function(x) mean(x)))
# 
# 
# # format myresults as dataframes for plotting
# own.diet <- data.frame(matrix(unlist(myresults$Master.own.diet), nrow = length(common_species), byrow = F))
# rownames(own.diet) <- common_species
# mynames <- gsub("_tax.txt", "", temp)
# mynames <- gsub("MCTs", "", mynames)
# colnames(own.diet) <- mynames
# 
# others.diet <- data.frame(matrix(unlist(myresults$Master.others.diet), nrow = length(common_species), byrow = F))
# rownames(others.diet) <- common_species
# colnames(others.diet) <- mynames
# 
# own.microbiome <- data.frame(matrix(unlist(myresults$Master.own.microbiome), nrow = length(common_species), byrow = F))
# rownames(own.microbiome) <- common_species
# colnames(own.microbiome) <- mynames
# 
# 
# 
# 
# # format for ggplot to make plots
# require(ggplot2)
# require(reshape2)
# 
# own.diet.melt <- rownames_to_column(own.diet, var = "Taxa")
# own.diet.melt <- melt(own.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.corr")
# 
# others.diet.melt <- rownames_to_column(others.diet, var = "Taxa")
# others.diet.melt <- melt(others.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "others.diet.corr")
# 
# own.microbiome.melt <-rownames_to_column(own.microbiome, var = "Taxa")
# own.microbiome.melt <- melt(own.microbiome.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.mb.corr")
# 
# master.df <- merge(own.diet.melt,others.diet.melt)
# master.df <- merge(own.microbiome.melt, master.df)
# 
# 
# master.df.melt <- melt(master.df, id.vars = c("Taxa", "UserName"))
# 
# # get summary stats for plotting error bars
# summary.master <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "mean")
# colnames(summary.master) <- c("variable", "mean")
# summary.sd <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "sd")
# colnames(summary.sd) <- c("variable", "sd")
# 
# summary.master <- merge(summary.master, summary.sd)
# summary.master$se <- summary.master$sd/sqrt(28) # se is sd/sqrt(n)
# 
# # make the plot
# mybarplot <- ggplot(data = master.df.melt, aes(x = variable, y = value)) +
#   geom_bar(stat = "summary", fun.y = mean, aes(fill = as.factor(variable))) +
#   geom_errorbar(data = summary.master, aes(x = variable, ymin = mean - se, ymax = mean + se), width = 0.4, inherit.aes = F) +
#   theme_classic()
# 
# # pring the plot
# mybarplot
# 
# 
# 
# # make paired difference plots by person
# 
# require(gridExtra)
# 
# mb.dt <- master.df.melt[master.df.melt$variable %in% c("own.mb.corr", "own.diet.corr"),]
# 
# mb.dt.wide <- reshape(mb.dt, idvar = c("Taxa", "UserName"), direction = 'wide', timevar = "variable")
# mb.dt.wide$diff <- mb.dt.wide$value.own.diet.corr - mb.dt.wide$value.own.mb.corr
# mb.dt.wide$diffquant <- cut(mb.dt.wide$diff, breaks = 4, labels = F)
# 
# mb.dt.for.merge <- mb.dt.wide[colnames(mb.dt.wide) %in% c("Taxa", "UserName", "diffquant")]
# 
# mb.dt <- merge(mb.dt, mb.dt.for.merge)
# # mb.dt$color <- ifelse(mb.dt$diffquant >=3, "color", "no")
# # mb.dt <- mb.dt[order(mb.dt$color, decreasing = F),]
# 
# 
# # make multiple plots in a loop
# 
# single_plots <- as.list(NULL)
# 
# plot <- mb.dt[order(mb.dt$diffquant, decreasing = T),]
# plot <- plot[order(plot$UserName),]
# 
# 
# for (i in 1:length(unique(plot$UserName))) {
#   j <- as.character(unique(plot$UserName)[[i]])
#   plot_un <- plot[plot$UserName == j,]
#   
#   single_plots[[i]] <- ggplot(data = plot_un, aes(x = variable, y = value)) + 
#     geom_point() +
#     geom_line(aes(group = Taxa, color = factor(diffquant))) +
#     scale_color_manual(values = c("grey", "grey", "red", "red")) +
#     theme_classic() +
#     theme(legend.position = "none",
#           legend.text = element_text(size = 7),
#           legend.title = element_blank()) 
#   
# }
# 
# 
# grid.arrange(grobs = single_plots[1:28], nrow=7)
# 
# 
# 
# 
