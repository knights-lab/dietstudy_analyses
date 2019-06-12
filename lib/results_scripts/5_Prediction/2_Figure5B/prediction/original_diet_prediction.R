
## Soylent participants and participants with less than 5 longitudinal timepoints are excluded ##
## n = 29 ##

setwd(dir = "/Users/abby/Documents/Projects/dietstudy_analyses/")
require(ape)
require(tidyverse)
# 
# # load map
# map <- read.delim(file = "data/maps/SampleID_map.txt")
# map <- map[colnames(map) %in% c("X.SampleID", "UserName","StudyDayNo","Sample.Week.Day")]
# 
# # get soylent ids
# soylent <- map[map$UserName %in% c("MCTs11", "MCTs12"),]
# soylentIDs <- droplevels(soylent$X.SampleID)
# 
# 
# # load diet 3 day beta div for all people
# diet3daydist <- read.delim("data/diet/processed_food/3dayfood_beta/unweighted_unifrac_master3daydiet.txt", row.names = 1)
# #diet1daydist <- read.delim("data/diet/processed_food/dhydrt_beta/unweighted_unifrac_dhydrt.txt", row.names = 1)
# # drop soylents 
# diet3daydist <- diet3daydist[!rownames(diet3daydist) %in% soylentIDs, !colnames(diet3daydist) %in% soylentIDs]
# 
# # get diet pcs
# diet3daypc <- as.data.frame(pcoa(diet3daydist)$vectors)
# 
# # load mb table for everyone
# mb <- read.delim("data/microbiome/processed_sample/taxonomy_norm_s.txt", row = 1)
# # drop soylents
# mb <- mb[,!colnames(mb) %in% soylentIDs]
# 
# # get mb pcs
# mbdist <- dist(t(mb))
# mbpc <- as.data.frame(pcoa(mbdist)$vectors)
# 
# # limit mb and mbpc to just the people with diet data
# mb <- mb[,colnames(mb) %in% rownames(diet3daypc)]
# mbpc <- mbpc[rownames(mbpc) %in% rownames(diet3daypc),]
# 
# # limit diet to just the days with mb data
# diet3daypc <- diet3daypc[rownames(diet3daypc) %in% rownames(mbpc),]
# 
# # limit map to just these people too
# map <- map[map$X.SampleID %in% rownames(mbpc),]
# map <- droplevels(map)
# 
# # limit both diet and mb pcs to the top 5 pcs
# mbpc <- mbpc[,1:5]
# diet3daypc <- diet3daypc[,1:5]
# 
# # fix naming
# colnames(mbpc) <- paste0("Mb.Axis.", 1:5)
# colnames(diet3daypc) <- paste0("Dt.Axis.", 1:5)
# 
# # clean up naming in the mb df
# # fix the naming of each species
# 
# 
# fix.names <- function(x) {
#   x <- gsub("?.*s__", "", x)
#   x <- gsub("?.*g__", "Uncl. Genus ", x)
#   x <- gsub("?.*f__", "Uncl. Family ", x)
#   x <- gsub("?.*o__", "Uncl. Order ", x)
#   x <- gsub("?.*p__", "Uncl. Phylum ", x)
#   x <- gsub("?.*k__", "Uncl. Kingdom ", x)
#   x <- gsub(";NA", "", x)
#   x <- gsub("\\[", "", x)
#   x <- gsub("\\]", "", x)
#   x <- gsub(" ", "_", x)
#   x <- gsub("-", "_", x)
#   x <- gsub("\\/", "_", x)
#   return(x)
# }
# 
# # store taxonomy strings for later use
# taxonomy <- as.data.frame(rownames(mb))
# colnames(taxonomy) <- "taxastring"
# 
# # fix naming in mb
# rownames(mb) <- fix.names(rownames(mb))
# 
# # fix naming in the reference taxonomy strings table
# taxonomy$mbname <- fix.names(taxonomy$taxastring)
# 
# save(taxonomy, file = "data/test_personal_diet_taxastrings.Rdata")
# 
# # sort mb by most variable
# mb.means <- rowMeans(mb)
# mb.sd <- apply(mb, 1, function(x) sd(x))
# mb.var <- mb.sd/abs(mb.means)
# 
# # get top 100 most common
# most.common <- names(head(sort(mb.means, decreasing = T), n = 30))
# 
# 
# #species.to.test <- c(most.common, most.var)
# species.to.test <- c(most.common)
# 
# # subset mb to these and transpose for later analysis
# mb.of.interest <- as.data.frame(t(mb[species.to.test,]))
# 
# save(species.to.test, file = "data/species.to.test.Rdata")
# 
# 
# # listify by person for the predictive modeling #
# 
# 
# list.mb<- NULL
# list.mbpc <- NULL
# list.dietpc <- NULL
# 
# for (i in seq_along(unique(map$UserName))){
#   a <- as.character(unique(map$UserName)[i])
#   ids <- droplevels(map$X.SampleID[map$UserName == a])
#   
#   sub.mb <- mb.of.interest[rownames(mb.of.interest) %in% ids,]
#   list.mb[[i]] <- sub.mb
# 
#   sub.mbpc <- mbpc[rownames(mbpc) %in% ids,]
#   list.mbpc[[i]] <- sub.mbpc
#   
#   sub.dietpc <- diet3daypc[rownames(diet3daypc) %in% ids,]
#   list.dietpc[[i]] <- sub.dietpc
#   
# } 
# 
# 
# 
# 
# # join the microbiome and diet dataframes and only keep days we have data for both
# predictors <- lapply(1:length(list.mbpc), function(i) merge(list.mbpc[[i]], list.dietpc[[i]], by = 0))
# 
# 
# ################################
# # we care about the previous day ra for the model, so add these to the predictors
# ra.predictors <- lapply(list.mb, function(x) as.data.frame(x))
# ra.predictors <- lapply(ra.predictors, function(x) {
#   x$Row.names <- rownames(x)
#   return(x)
# })
# 
# # we also care about these as the response variables
# response <- lapply(list.mb, function(x) as.data.frame(x))
# response <- lapply(response, function(x) {
#   x$Row.names <- rownames(x)
#   return(x)
# })
# 
# 
# # add the day variable to both predictors and response variables
# predictors_day <- lapply(1:length(predictors), function(i) merge(map, predictors[[i]], by.x = "X.SampleID", by.y = "Row.names", all =F))
# response_day <- lapply(1:length(response), function(i) merge(map, response[[i]], by.x = "X.SampleID", by.y = "Row.names", all = F))
# 
# # get the indicies to use for filtering
# possible.days <- c(1:17) # vector to compare to 1:17 becuase there are 17 possible days
# pred_days <- lapply(1:length(predictors_day), function(i) possible.days %in% predictors_day[[i]]$StudyDayNo) # this gives the list of days we values for
# 
# # make a list of the dataframes made from the vector
# possible.days.df <- as.data.frame(possible.days, drop = F)
# 
# possible.days.list <- lapply(1:length(pred_days), function(i) {cbind(possible.days.df, pred_days[[i]])})
# possible.days.list <- lapply(possible.days.list, function(x){
#   colnames(x) <- c("possible.days", "pred.days")
#   return(x)
# })
# possible.days.list <- lapply(possible.days.list, function(x){
#   x$pred_days <- ifelse(x$pred.days == TRUE, x$possible.days, NA)
#   x$resp_days <- x$pred_days + c(diff(x$pred_days), NA)
#   x$pred_days <- x$resp_days - 1
#   keeps <- c("possible.days", "pred_days", "resp_days")
#   x <- x[keeps]
#   return(x)
# })
# 
# # possible.days.list can now act as a framework to add the appropriate prediction and response days
# 
# # merge to get all the information in one table that we will use for testing
# # either boosted regression or randomforests
# dat <- lapply(1:length(possible.days.list), function(i){merge(possible.days.list[[i]], predictors_day[[i]], by.x = "pred_days", by.y = "StudyDayNo")})
# dat <- lapply(1:length(dat), function(i){merge(dat[[i]], ra.predictors[[i]], by.x = "X.SampleID", by.y = "Row.names")})
# dat <- lapply(1:length(dat), function(i){merge(dat[[i]], response_day[[i]], by.x = "resp_days", by.y = "StudyDayNo")})
# 
# # addnames
# mynames <- as.character(unique(map$UserName))
# 
# dat <- setNames(dat, mynames)
# 
# # identify people with < 5 time points
# drops <- which(lapply(dat, function(x) dim(x)[1])<=5)
# 
# # drop the people with too few time points
# dat <- dat[!names(dat) %in% names(drops)]
# 
# save(dat, file = "data/test_personal_diet_dat.Rdata")

################### TEST DIET IMPACT ###########
load(file = "data/test_personal_diet_dat.Rdata")
load(file = "data/species.to.test.Rdata")

species <- keepspecies


"test.diet.impact" <- function(){
  set.seed(42)
  
  Master.own.diet <- NULL
  Master.others.diet <- NULL
  Master.own.microbiome <- NULL
  Master.scramble.diet <- NULL
  
  
  for (n in 1:length(dat)) {
    
    a <- names(dat[n])
    lengtha <- dim(dat[[a]])[1]
    
    #person B index
    possible.b <- names(dat)
    possible.b <- possible.b[possible.b != a]
    
    alllength <- lapply(dat,function(x) dim(x)[1])
    alllength <- alllength[possible.b]
    samplefrom <- possible.b[which(alllength >= lengtha)]
    
    b <- sample(samplefrom,1)
    
    # get the species names
    today.species <- paste0(species, ".x")
    tomorrow.species <- paste0(species, ".y")
    
    # create empty list for storing values
    own.diet.predictions <- NULL
    others.diet.predictions <- NULL
    own.microbiome.predictions <- NULL
    scramble.diet.predictions <- NULL
    
    
    for (i in seq_along(tomorrow.species)) {
      
      tomorrow <- tomorrow.species[i]
      today <- today.species[i] 
      
      mb.features <- c("Mb.Axis.1","Mb.Axis.2", "Mb.Axis.3", "Mb.Axis.4", "Mb.Axis.5")
      dt.features <- c("Dt.Axis.1","Dt.Axis.2", "Dt.Axis.3", "Dt.Axis.4", "Dt.Axis.5")
      
      # predict person A's microbiome from microbiome alone
      # species tomorrow from today's species and microbiome mbA.feats
      # predict tomorrow 
      mbA.feats <- dat[[a]][,mb.features]
      actual.tomorrow <- dat[[a]][,tomorrow]
      mb.only.model <- lm(actual.tomorrow ~., mbA.feats)
      pred.tomorrow <- predict(mb.only.model)
      #pred.tomorrow <- dat[[a]][,today]  ###TODO: test this to see if it's about the same as the regression result
      # compare predicted values to actual values for microbiome only
      own.microbiome <- cor(pred.tomorrow, actual.tomorrow)
      
      
      # get residuals between predicted species tomorrow and actual species tomorrow
      resid.tomorrow <- pred.tomorrow - actual.tomorrow
      
      # regress the residuals on diet
      dtA.feats <- dat[[a]][,dt.features]
      dt.mb.model <- lm(resid.tomorrow ~., dtA.feats)
      pred.resid <- predict(dt.mb.model)
      
      # subtract predicted residuals (dt.mb.model) from yhat
      # and correlate with acutal value for tomorrow
      own.diet <- cor(pred.tomorrow - pred.resid, actual.tomorrow)
      
      
      # model the residuals from scrambled diet features
      # scramble the diet features
      scramble.dtA.feats <- dtA.feats[sample(nrow(dtA.feats)),]
      scrambe.modelA <- lm(resid.tomorrow~.,scramble.dtA.feats)
      
      # now predict using the model for residuals trained on normal diet
      # but pass in scrambled diet features 
      pred.scramble.resid <- predict(scrambe.modelA, dtA.feats)
      #pred.scramble.resid <- predict(scrambe.modelA)
      # determine the correlation
      scramble.diet <- cor(pred.tomorrow - pred.scramble.resid, actual.tomorrow)
      
      ##
      # predict from person B
      # predict person B's species from B microbiome alone
      mbB.feats <- dat[[b]][,mb.features]
      actual.tomorrowB <- dat[[b]][,tomorrow]
      mb.only.modelB <- lm(actual.tomorrowB ~., mbB.feats)
      pred.tomorrowB <- predict(mb.only.modelB)
      #pred.tomorrowB <- dat[[b]][,today]
      
      # get residuals for B
      resid.tomorrowB <- pred.tomorrowB - actual.tomorrowB
      
      #run predict on the model trained with the data from person B
      #but passing in the diet data for person A
      #to get the externally-predicted residuals
      dtB.feats <- dat[[b]][,dt.features]
      dt.mb.modelB <- lm(resid.tomorrowB~.,dtB.feats)
      pred.residB.from.dtA <- predict(dt.mb.modelB, dtA.feats)
      
      #subtract yr.B.hat from yhat.B and correlate the final result with y.B.
      others.diet <- cor(pred.tomorrowB - pred.residB.from.dtA, actual.tomorrowB)
      
      
      # store variables
      own.diet.predictions[[tomorrow]] <- own.diet
      others.diet.predictions[[tomorrow]] <- others.diet
      own.microbiome.predictions[[tomorrow]] <- own.microbiome
      scramble.diet.predictions[[tomorrow]] <- scramble.diet
    }
    Master.own.diet[[n]] <- own.diet.predictions
    Master.others.diet[[n]] <- others.diet.predictions
    Master.own.microbiome[[n]] <- own.microbiome.predictions
    Master.scramble.diet[[n]] <- scramble.diet.predictions
    
    
  }
  
  return(list(Master.own.diet = Master.own.diet, 
              Master.others.diet = Master.others.diet, 
              Master.own.microbiome = Master.own.microbiome,
              Master.scramble.diet = Master.scramble.diet))
  
}


myresults <- test.diet.impact()

save(myresults, file ='data/original.prediction.results.Rdata')

