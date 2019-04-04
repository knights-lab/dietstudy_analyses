
## Soylent participants and participants with less than 5 longitudinal timepoints are excluded ##
## n = 28 ##

setwd(dir = "/Users/abby/Documents/Projects/dietstudy_analyses/")
require(ape)
require(tidyverse)

load(file = "data/test_personal_diet_dat.Rdata")
load(file = "data/species.to.test.Rdata")
################### TEST DIET IMPACT ###########


species <- species.to.test


"leave.one.out" <- function(){
  set.seed(42)
  
  Master.own.diet <- NULL
  Master.own.microbiome <- NULL
  Master.scramble.diet <- NULL
  Master.others.diet <- NULL
  
  for (n in 1:length(dat)) {
    # person A index
    a <- n
    
    # person B index
    possible.b <- 1:length(dat)
    possible.b <- possible.b[possible.b != a]
    b <- sample(possible.b,1)
    
    # get the species names
    today.species <- paste0(species, ".x")
    tomorrow.species <- paste0(species, ".y")
    
    # create empty list for storing values
    own.diet.predictions <- NULL
    own.microbiome.predictions <- NULL
    scramble.diet.predictions <- NULL
    others.diet.predictions <- NULL
    
    for (i in seq_along(tomorrow.species)) {
      
      tomorrow <- tomorrow.species[i]
      today <- today.species[i] 
      
      mb.features <- c(today, "Mb.Axis.1","Mb.Axis.2", "Mb.Axis.3")
      dt.features <- c("Dt.Axis.1","Dt.Axis.2", "Dt.Axis.3")
      
      # obtain variables from models
      mbA.feats <- dat[[a]][,mb.features]
      actual.tomorrow <- dat[[a]][,tomorrow]
      dtA.feats <- dat[[a]][,dt.features]
      scramble.dtA.feats <- dtA.feats[sample(nrow(dtA.feats)),]
      mbB.feats <- dat[[b]][,mb.features]
      actual.tomorrowB <- dat[[b]][,tomorrow]
      dtB.feats <- dat[[b]][,dt.features]
      
      # create values to store variables
      pred.tomorrow<-NULL
      resid.tomorrow<-NULL 
      
      # hold one out (HOO) prediction of microbiome from mb, diet, scrambled diet, and other's diet
      for (j in 1:dim(mbA.feats)[1]) {
        
        mbA.feats_j <- mbA.feats[-j,] # make reduced features 
        actual.tomorrow_j <- actual.tomorrow[-j] # make reduced actual value
        mb.only.model <- lm(actual.tomorrow_j ~., mbA.feats_j)
        
        # predict person A's microbiome from microbiome alone
        pred.tomorrow_j <- as.numeric(predict(mb.only.model, (mbA.feats_j = mbA.feats[j,])))# predict the jth obs
        pred.tomorrow[[j]] <- pred.tomorrow_j
        
        # calcuate residuals
        resid.tomorrow_j <- pred.tomorrow_j - actual.tomorrow[j]
        resid.tomorrow[[j]] <- resid.tomorrow_j
        
      }
        
      pred.tomorrowB <- NULL
      resid.tomorrowB <- NULL

      
      for (l in 1:dim(mbB.feats)[1]){
        # predict for B (here the issue is that B wont always have the same number of features)
        mbB.feats_l <- mbB.feats[-l,] # make reduced features 
        actual.tomorrowB_l <- actual.tomorrowB[-l] # make reduced actual value
        mb.only.modelB <- lm(actual.tomorrowB_l ~., mbB.feats_l)
        pred.tomorrowB_l <- as.numeric(predict(mb.only.modelB, (mbB.feats_l = mbB.feats[l,])))
        pred.tomorrowB[[l]] <- pred.tomorrowB_l
        
        # get residuals for B
        resid.tomorrowB_l <- pred.tomorrowB_l - actual.tomorrowB[l]
        resid.tomorrowB[[l]] <- resid.tomorrowB_l
        
      }
      
      
      pred.resid <- NULL
      pred.scramble.resid <- NULL
      dt.mb.model <- vector(mode = "list", length = dim(mbA.feats)[1])
      
      for(k in 1:dim(mbA.feats)[1]) {
        
        # regress the residuals on diet (have to have stored residuals first)
        dtA.feats_k <- dtA.feats[-k,]
        resid.tomorrow_k <- resid.tomorrow[-k]
        dt.mb.model <- lm(resid.tomorrow_k ~., dtA.feats_k)
        pred.resid_k <- as.numeric(predict(dt.mb.model, (dtA.feats_k = dtA.feats[k,])))
        pred.resid[[k]] <- pred.resid_k
        
        # # store the models because we use them later
        # dt.mb.model[[k]] <- dt.mb.model
        
        # model the residuals from scrambled diet features (again, have to have stored residuals first)
        scramble.dtA.feats_k <- scramble.dtA.feats[-k,] # make reduced features 
        scramble.model <- lm(resid.tomorrow_k ~., scramble.dtA.feats_k)
        pred.scramble.resid_k <- as.numeric(predict(scramble.model, (scramble.dtA.feats_k = scramble.dtA.feats[k,])))
        pred.scramble.resid[[k]] <- pred.scramble.resid_k
        
      }


      pred.resid.from.dtB <- NULL
      
      for (m in 1:dim(mbB.feats)[1]) {
      
        #run predict on the model trained with the data from person A
        #but passing in the diet data for person B 
        #to get the externally-predicted residuals
        dtA.feats_m <- dtA.feats[-m,]
        resid.tomorrow_m <- resid.tomorrow[-m]
        dt.mb.model <- lm(resid.tomorrow_m ~., dtA.feats_m)
        pred.resid_m <- as.numeric(predict(dt.mb.model, (dtA.feats_m = dtA.feats[m,])))
        
        dtB.feats_m <- dtB.feats[-m,] # make reduced features 
        pred.resid.from.dtB_m <- as.numeric(predict(dt.mb.model, (dtA.feats_m = dtB.feats[m,])))
        pred.resid.from.dtB[[m]] <- pred.resid.from.dtB_m
        
      }
 
      ##Evaluate the HOO models###
      # calculate correlation between hold one out mb predictions and actual values
      own.microbiome <- cor(pred.tomorrow, actual.tomorrow)
      
      # subtract predicted residuals
      # and correlate with acutal value for tomorrow
      own.diet <- cor(pred.tomorrow - pred.resid, actual.tomorrow)

      # determine the correlation with scrambled diet
      scramble.diet <- cor(pred.tomorrow - pred.scramble.resid, actual.tomorrow)
 
      #subtract yr.B.hat from yhat.B and correlate the final result with y.B.
      others.diet <- cor(pred.tomorrowB - pred.resid.from.dtB, actual.tomorrowB)
      
      # store these values
      own.microbiome.predictions[[tomorrow]] <- own.microbiome
      own.diet.predictions[[tomorrow]] <- own.diet
      scramble.diet.predictions[[tomorrow]] <- scramble.diet
      others.diet.predictions[[tomorrow]] <- others.diet
      
    }
    Master.own.diet[[n]] <- own.diet.predictions
    Master.own.microbiome[[n]] <- own.microbiome.predictions
    Master.scramble.diet[[n]] <- scramble.diet.predictions
    Master.others.diet[[n]] <- others.diet.predictions
    
    
  }
  
  return(list(Master.own.diet = Master.own.diet, 
              Master.own.microbiome = Master.own.microbiome,
              Master.scramble.diet = Master.scramble.diet,
              Master.others.diet = Master.others.diet))
  
}
      


test <- leave.one.out()



save(myresults, file ='data/prediction.results.Rdata')

