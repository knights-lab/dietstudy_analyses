
## Soylent participants and participants with less than 5 longitudinal timepoints are excluded ##
## n = 29 ##

setwd(dir = "/Users/abby/Documents/Projects/dietstudy_analyses/")
require(ape)
require(tidyverse)
require(caret)

load(file = "data/test_personal_diet_dat.Rdata")
load(file = "data/species.to.test.Rdata")

################### TEST DIET IMPACT ###########


species <- keepspecies



"leave.one.out" <- function(){
  set.seed(42)
  
  Master.own.microbiome <- NULL
  Master.own.diet <- NULL
  Master.scramble.diet <- NULL
  Master.others.diet <- NULL
  
  for (n in 1:length(dat)) {
    # person A index
    a <- names(dat[n])
    lengtha <- dim(dat[[a]])[1]
    
    #person B index
    possible.b <- names(dat)
    possible.b <- possible.b[possible.b != a]

    alllength <- lapply(dat,function(x) dim(x)[1])
    alllength <- alllength[possible.b]
    samplefrom <- possible.b[which(alllength <= lengtha)]

    b <- sample(samplefrom,1)
    
    # get the species names
    #today.species <- colnames(dat[[a]])[grep("\\.x", colnames(dat[[a]]))]
    #today.species <- today.species[-c(1:3)]
    #tomorrow.species <- colnames(dat[[a]])[grep("\\.y", colnames(dat[[a]]))]
    #tomorrow.species <- tomorrow.species[-c(1:3)]
    today.species <- paste0(species, ".x")
    tomorrow.species <- paste0(species, ".y")
    
    # create empty list for storing values
    own.microbiome.predictions <- NULL
    own.diet.predictions <- NULL
    scramble.diet.predictions <- NULL
    others.diet.predictions <- NULL
    
    for (i in seq_along(tomorrow.species)) {
      
      tomorrow <- tomorrow.species[i]
      today <- today.species[i] 
      
      dt.features <- c("Dt.Axis.1","Dt.Axis.2", "Dt.Axis.3", "Dt.Axis.4", "Dt.Axis.5")
      mb.features <- c("Mb.Axis.1", "Mb.Axis.2", "Mb.Axis.3", "Mb.Axis.4", "Mb.Axis.5")
      
      # obtain variables for training models
      mb.residuals.A <- dat[[a]][,tomorrow] - dat[[a]][,today]
      mbA.feats <- dat[[a]][,mb.features]
      #actual.tomorrowA <- dat[[a]][,tomorrow]
      dtA.feats <- dat[[a]][,dt.features]
      scramble.dtA.feats <- dtA.feats[sample(nrow(dtA.feats)),]
      
      #mbB.feats <- dat[[b]][,mb.features]
      #actual.tomorrowB <- dat[[b]][,tomorrow]
      dtB.feats <- dat[[b]][,dt.features]
      #mb.residuals.B <- dat[[b]][,tomorrow] - dat[[b]][,today]

      # use caret to run the LOOCV models
      #predict mb residual from microbiome features
      mySummary <- function(data, lev = NULL, model = NULL) {
        out <- cor(data$obs, data$pred)
        names(out) <- "COR"
        out
      }

      control <- trainControl(method = "LOOCV", summaryFunction = mySummary)
      #control <- trainControl(method = "LOOCV")
      method <- "glm"
      d.A.mb <- NULL
      d.A.mb <- data.frame(mbA.feats, y = mb.residuals.A)
      #mtry <- sqrt(ncol(mbA.feats))
      #tunegrid <- expand.grid(.mtry = mtry)
      mod.A.mb <- train(y ~.,
                        method = method,
                        data = d.A.mb,
                        #tuneGrid = tunegrid,
                        trControl = control)
      # 
      # predict mb residual from microbiome and dietary features 
      d.A <- NULL
      d.A <- data.frame(dtA.feats, mbA.feats, y = mb.residuals.A)
      #mtry <- sqrt(ncol(data.frame(dtA.feats,mbA.feats)))
      #tunegrid <- expand.grid(.mtry = mtry)
      mod.A <- train(y ~., 
                     method = method, 
                     data = d.A, 
                     #tuneGrid = tunegrid,
                     trControl = control)
      
      
      # predict mb residual from scrambled dietary features
      d.mc.A <- NULL
      d.mc.A <- data.frame(scramble.dtA.feats, mbA.feats, y = mb.residuals.A)
      #mtry <- sqrt(ncol(data.frame(scramble.dtA.feats,mbA.feats)))
      #tunegrid <- expand.grid(.mtry = mtry)
      mod.mc.A <-  train(y ~., 
                         method = method, 
                         data = d.mc.A, 
                         #tuneGrid = tunegrid,
                         trControl = control)
      
      
      # # predict from someone elses diet features
      l <- dim(dtB.feats)[1] # length of samples needed
      d.B <- NULL
      d.B <- data.frame(dtB.feats, mbA.feats[1:l,], y = mb.residuals.A[1:l])
      #mtry <- sqrt(ncol(data.frame(dtB.feats,mbA.feats[1:l,])))
      #tunegrid <- expand.grid(.mtry = mtry)
      mod.B <- train(y ~., 
                     method = method, 
                     data = d.B, 
                     #tuneGrid = tunegrid,
                     trControl = control)
      

      # store these values
      own.microbiome.predictions[[tomorrow]] <- mod.A.mb$results$COR
      own.diet.predictions[[tomorrow]] <- mod.A$results$COR
      scramble.diet.predictions[[tomorrow]] <- mod.mc.A$results$COR
      others.diet.predictions[[tomorrow]] <- mod.B$results$COR

    }
    
    Master.own.microbiome[[n]] <- own.microbiome.predictions
    Master.own.diet[[n]] <- own.diet.predictions
    Master.scramble.diet[[n]] <- scramble.diet.predictions
    Master.others.diet[[n]] <- others.diet.predictions
    
    
  }
  
  return(list(Master.own.microbiome = Master.own.microbiome,
              Master.own.diet = Master.own.diet, 
              Master.scramble.diet = Master.scramble.diet,
              Master.others.diet = Master.others.diet))
  
}
      


test <- leave.one.out()



save(test, file ='data/LOOCV.prediction.results.Rdata')

