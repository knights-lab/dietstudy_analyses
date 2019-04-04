

# The input for this analysis is a list of dataframes, 'dat'
# 'dat' contains one person's data per dataframe
# each dataframe has been arranged to pair predictor and response variables within the same row


# key variables needed and included for each person are:
### Predictors ####
# ## microbiome pricipal coordinates for today (Mb.Axis 1 - 5)
# ## diet principal coordinates for the 3-days prior to today (Dt.Axis 1 - 5)
# ## species abundance data for today (species_name.x) - predictor

### Response ###
# ## species abundance data for tomorrow (species_name.y) - response


# the character vector 'species' contains the list of the top 30 most abundant species to loop through and test
# these are the 30 predictor and response species that are included in each person's dataframe

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