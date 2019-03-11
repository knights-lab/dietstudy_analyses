library(dplyr)
############### THIS CODE TAKES THE RAW ASA24 OUTPUT AND CLEANS IT ###################

# Import diet data files
# Description below file import line

 
Items <- read.csv("/Users/abby/Documents/Projects/MCTs/Data/All Diet Data/MCTs_23887_Items.csv", 
                  sep = ",", 
                  stringsAsFactors = FALSE)

# FNDDS food codes, gram weights, food pattern equivalents, for each item reported

INS <- read.csv("~/Documents/Projects/MCTs/Data/All Diet Data/MCTs_23888_INS.csv", 
                sep = ",",
                stringsAsFactors = FALSE)
# Individual supplement analysis file, supplement codes with their nutrients for each supplement reported
# only includes supplements found in FNDDS

# Note before loading responses will need to make sure that there are no ',' that mess up loading
Responses <- read.csv("~/Documents/Projects/MCTs/Data/All Diet Data/MCTs_23889_Responses.csv", 
                      sep = ",", 
                      stringsAsFactors = FALSE)
# Food and supplement names from the quick list, probe questions and answers

Totals <- read.csv("~/Documents/Projects/MCTs/Data/All Diet Data/MCTs_23890_Totals.csv", 
                   sep = ",",
                   stringsAsFactors = FALSE)
# Daily total of FNDDS nutrients and FPED food groups for a consumption day

TNS <- read.csv("~/Documents/Projects/MCTs/Data/All Diet Data/MCTs_23891_TNS.csv", 
                sep = ",", 
                stringsAsFactors = FALSE)
# Daily total nutrients for foods and supplements analysis file 
# FNDDS nutrients from all foods and supplements reported

TS <- read.csv("~/Documents/Projects/MCTs/Data/All Diet Data/MCTs_23892_TS.csv", 
               sep = ",", 
               stringsAsFactors = FALSE)
# Daily total supplements analysis file - total nutrient intake from all supplements consumed per day

######################################################

###### RENAME MAKE-UP DATA FROM SUBJECTS WHO COULDN'T RECORD DATA IN ORGINAL RECORD ###########

# Notes on diet files that need to be matched up correctly with participants/dates
# For each we need to fix UserName, UserID, RecordDayNo, IntakeStartDateTime, IntakeEndDateTime, 
# RecordDay, ReportingDate

Items$UserName <- as.factor(Items$UserName)
users <- levels(Items$UserName)

# limit Items.fixed to just real usernames - will add back the edited values
Items.correct.ids <- Items %>% filter(UserName %in% users[1:37])

# Fix the following in the Items file

# WrongID  Date   CorrectID
# MCTs46 = 2/13 = MCTs06 (RecordDay 15)

MCTs46 <- Items %>% filter(UserName == "MCTs46")

MCTs46$IntakeStartDateTime <- lapply(MCTs46$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/15/2017", 
                                     replacement = "02/13/2017")

MCTs46$IntakeEndDateTime <- lapply(MCTs46$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/16/2017", 
                                   replacement = "02/14/2017")

MCTs46$ReportingDate <- "AfterStudyEnd"

MCTs46$UserID <- "9a4f2d22-663b-40e7-a826-11a4c87dc7f3"

MCTs46$UserName <- "MCTs06"

MCTs46$RecordDayNo <- 15


# MCTs47 = 2/14 = MCTs13 (Day 16)

MCTs47 <- Items %>% filter(UserName == "MCTs47")

MCTs47$IntakeStartDateTime <- lapply(MCTs47$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/15/2017", 
                                     replacement = "02/14/2017")

MCTs47$IntakeEndDateTime <- lapply(MCTs47$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/16/2017", 
                                   replacement = "02/15/2017")

MCTs47$ReportingDate <- "AfterStudyEnd"

MCTs47$UserID <- "7b74b8f3-390f-48d9-a01e-8dae6a8b38a7"

MCTs47$UserName <- "MCTs13"

MCTs47$RecordDayNo <- 16



# MCTs48 = 2/14 = MCTs28 (Day 16)

MCTs48 <- Items %>% filter(UserName == "MCTs48")

MCTs48$IntakeStartDateTime <- lapply(MCTs48$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/16/2017", 
                                     replacement = "02/14/2017")

MCTs48$IntakeEndDateTime <- lapply(MCTs48$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/17/2017",
                                   replacement = "02/15/2017")

MCTs48$ReportingDate <- "AfterStudyEnd"

MCTs48$UserID <- "74e16c6a-8003-471d-ac1b-50c2f91b4e4c"

MCTs48$UserName <- "MCTs28"

MCTs48$RecordDayNo <- 16



# MCTs49 = 2/15 = MCTs28 (Day 17)

MCTs49 <- Items %>% filter(UserName == "MCTs49")

MCTs49$IntakeStartDateTime <- lapply(MCTs49$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/16/2017", 
                                     replacement = "02/15/2017")
MCTs49$IntakeEndDateTime <- lapply(MCTs49$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/17/2017", 
                                   replacement = "02/16/2017")

MCTs49$ReportingDate <- "AfterStudyEnd"

MCTs49$UserID <- "74e16c6a-8003-471d-ac1b-50c2f91b4e4c"

MCTs49$UserName <- "MCTs28"

MCTs49$RecordDayNo <- 17



# MCTs50 = 2/11 = MCTs19 (day 13)

MCTs50 <- Items %>% filter(UserName == "MCTs50")

MCTs50$IntakeStartDateTime <- lapply(MCTs50$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/26/2017", 
                                     replacement = "02/11/2017")

MCTs50$IntakeEndDateTime <- lapply(MCTs50$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/27/2017", 
                                   replacement = "02/12/2017")

MCTs50$ReportingDate <- "AfterStudyEnd"

MCTs50$UserID <- "8df5f5ed-ac75-4f5b-bde2-f6499b5996d9"

MCTs50$UserName <- "MCTs19"

MCTs50$RecordDayNo <- 13



# MCTs51 = 2/12 = MCTs19 (day 14)

MCTs51 <- Items %>% filter(UserName == "MCTs51")

MCTs51$IntakeStartDateTime <- lapply(MCTs51$IntakeEndDateTime, 
                                     gsub, 
                                     pattern = "02/26/2017", 
                                     replacement = "02/12/2017")

MCTs51$IntakeEndDateTime <- lapply(MCTs51$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/27/2017", 
                                   replacement = "02/13/2017")

MCTs51$ReportingDate <- "AfterStudyEnd"

MCTs51$UserID <- "8df5f5ed-ac75-4f5b-bde2-f6499b5996d9"

MCTs51$UserName <- "MCTs19"

MCTs51$RecordDayNo <- 14


# MCTs52 = 2/13 = MCTs27 (Day 15)

MCTs52 <- Items %>% filter(UserName == "MCTs52")

MCTs52$IntakeStartDateTime <- lapply(MCTs52$IntakeEndDateTime, 
                                     gsub, 
                                     pattern = "02/16/2017", 
                                     replacement = "02/13/2017")

MCTs52$IntakeEndDateTime <- lapply(MCTs52$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/17/2017", 
                                   replacement = "02/14/2017")

MCTs52$ReportingDate <- "AfterStudyEnd"

MCTs52$UserID <- "e82b519d-39b5-462b-bf92-ff320cd38830"

MCTs52$UserName <- "MCTs27"

MCTs52$RecordDayNo <- 15


# MCTs53 = 2/14 = MCTs27 (Day 16)

MCTs53 <- Items %>% filter(UserName == "MCTs53")

MCTs53$IntakeStartDateTime <- lapply(MCTs53$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/16/2017", 
                                     replacement = "02/14/2017")

MCTs53$IntakeEndDateTime <- lapply(MCTs53$IntakeEndDateTime,
                                   gsub, 
                                   pattern = "02/17/2017", 
                                   replacement = "02/15/2017")

MCTs53$ReportingDate <- "AfterStudyEnd"

MCTs53$UserID <- "e82b519d-39b5-462b-bf92-ff320cd38830"

MCTs53$UserName <- "MCTs27"

MCTs53$RecordDayNo <- 16


# MCTs54 = 2/13 = MCTs04 (day 15)

MCTs54 <- Items %>% filter(UserName == "MCTs54")

MCTs54$IntakeStartDateTime <- lapply(MCTs54$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/16/2017", 
                                     replacement = "02/13/2017")

MCTs54$IntakeEndDateTime <- lapply(MCTs54$IntakeEndDateTime,
                                   gsub,
                                   pattern = "02/17/2017", 
                                   replacement = "02/14/2017")

MCTs54$ReportingDate <- "AfterStudyEnd"

MCTs54$UserID <- "5f034b96-5bc8-48ad-a1f6-b8f0e5855af7"

MCTs54$UserName <- "MCTs04"

MCTs54$RecordDayNo <- 15


# MCTs55 = 2/14 = MCTs04 (Day 16)

MCTs55 <- Items %>% filter(UserName == "MCTs55")

MCTs55$IntakeStartDateTime <- lapply(MCTs55$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/20/2017", 
                                     replacement = "02/14/2017")

MCTs55$IntakeEndDateTime <- lapply(MCTs55$IntakeEndDateTime, 
                                   gsub,
                                   pattern = "02/21/2017", 
                                   replacement = "02/15/2017")


MCTs55$ReportingDate <- "AfterStudyEnd"

MCTs55$UserID <- "5f034b96-5bc8-48ad-a1f6-b8f0e5855af7"

MCTs55$UserName <- "MCTs04"

MCTs55$RecordDayNo <- 16


# MCTs56 = 2/12 = MCTs36 (Day 14)

MCTs56 <- Items %>% filter(UserName == "MCTs56")

MCTs56$IntakeStartDateTime <- lapply(MCTs56$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/17/2017", 
                                     replacement = "02/12/2017")

MCTs56$IntakeEndDateTime <- lapply(MCTs56$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/18/2017", 
                                   replacement = "02/13/2017")

MCTs56$ReportingDate <- "AfterStudyEnd"

MCTs56$UserID <- "6deb8c57-878e-42f2-991a-b23743be3101"

MCTs56$UserName <- "MCTs36"

MCTs56$RecordDayNo <- 14

# MCTs57 = 2/13 = MCTs36 (Day 15)

MCTs57 <- Items %>% filter(UserName == "MCTs57")

MCTs57$IntakeStartDateTime <- lapply(MCTs57$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/19/2017", 
                                     replacement = "02/13/2017")

MCTs57$IntakeEndDateTime <- lapply(MCTs57$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/20/2017", 
                                   replacement = "02/14/2017")

MCTs57$ReportingDate <- "AfterStudyEnd"

MCTs57$UserID <- "6deb8c57-878e-42f2-991a-b23743be3101"

MCTs57$UserName <- "MCTs36"

MCTs57$RecordDayNo <- 15


# MCTs58 = 2/14 = MCTs36 (Day 16)

MCTs58 <- Items %>% filter(UserName == "MCTs58")

MCTs58$IntakeStartDateTime <- lapply(MCTs58$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/19/2017", 
                                     replacement = "02/14/2017")

MCTs58$IntakeEndDateTime <- lapply(MCTs58$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/20/2017", 
                                   replacement = "02/15/2017")

MCTs58$ReportingDate <- "AfterStudyEnd"

MCTs58$UserID <- "6deb8c57-878e-42f2-991a-b23743be3101"

MCTs58$UserName <- "MCTs36"

MCTs58$RecordDayNo <- 16


# MCTs59 = 2/15 = MCTs36 (Day 17)

MCTs59 <- Items %>% filter(UserName == "MCTs59")

MCTs59$IntakeStartDateTime <- lapply(MCTs59$IntakeStartDateTime, 
                                     gsub, 
                                     pattern = "02/19/2017", 
                                     replacement = "02/15/2017")

MCTs59$IntakeEndDateTime <- lapply(MCTs59$IntakeEndDateTime, 
                                   gsub, 
                                   pattern = "02/20/2017", 
                                   replacement = "02/16/2017")

MCTs59$ReportingDate <- "AfterStudyEnd"

MCTs59$UserID <- "6deb8c57-878e-42f2-991a-b23743be3101"

MCTs59$UserName <- "MCTs36"

MCTs59$RecordDayNo <- 17


################### COMBINE THE EDITED FILES WITH THE ORIGINAL DATA -incorrect names ##########
Items.fixed <- rbind(Items.correct.ids, 
                     MCTs46, MCTs47, MCTs48, MCTs49, MCTs50, MCTs51, MCTs52, 
                     MCTs53, MCTs54, MCTs55, MCTs56, MCTs57, MCTs58, MCTs59)

# remove the partial datasets that we don't need anymore
rm(list = ls(pattern = "^MCTs"))

################## RUN ABDI's CODE TO FIX SOYLENT VALUES#############


#When/if we get more details about the Soylent nutritional composition add it here
Soylent <- read.csv("~/Documents/Projects/MCTs/Data/Soylent.csv",header = TRUE, sep = ",")

# there is a miscode in MCTs12 where 95120010 was entered instead of 95120000
# need to replace just within MCTs12 for this value so it is properly coded as soylent

Items.fixed$FoodCode[which(Items.fixed$FoodCode == "95120010" & Items.fixed$UserName == "MCTs12")] <- "95120000"

items_to_change <- Items.fixed[Items.fixed$FoodCode == "95120000",]


# this gives us the responses that we need to edit

drinklist <- Responses[ Responses$Variable %in% 
                          c("BevSupplType") & Responses$Response %in% 
                          c("Other" , "Regular") & Responses$UserName %in% 
                          c("MCTs11","MCTs12"), ]

# this binds drinklist and test to give all items in responses needed to be changed
n <- cbind(drinklist,items_to_change)

n$FoodCode[n$ResponseOs %in% c("Soylent (Original)", "Soylent (original)", "Soylent (Origin)", "Solylent (Original)", NA, "")] <- 95120051
n$FoodCode[n$ResponseOs %in% c("Soylent (Nectar)", "Soylent (nectar)","nectar")] <- 95120052
n$FoodCode[n$ResponseOs %in% c("Soylent (Coffiest)", "Soylent (coffiest)", "coffiest")] <- 95120054
n$FoodCode[n$ResponseOs %in% c("Soylent (Cacao)", "Soylent (cacao)", "cacao", "Cacao", "Soylent (Cacao)")] <- 95120053

# this is the for loop that we really need
for(k in names(n)){
  for(m in names(Soylent)){
    if(k == m){
      n[k] <- Soylent[m][match(n$FoodCode, Soylent$FoodCode),]
    }
  }
}
# result of loop is n
# n has the correct food codes in columns 34:166 

# limits to just the correct columns
n1 <- n[34:166]

# this puts n1 at the end of the original items list
Items.final <- rbind(Items.fixed,n1)

# to get to the final data that we need we just removed the wrong food code
Items.final <- Items.final %>% filter(FoodCode != 95120000)

Items.final$UserName <- factor(Items.final$UserName)



######################## CREATE NEW TOTALS ####################


Totals.final <- Items.final %>% 
  select(UserName,RecordDayNo, 26:127) %>% 
  group_by(UserName, RecordDayNo) %>% 
  summarise_each(funs(sum(., na.rm = TRUE)))


# colnames(Totals.final)[2] <- "StudyDayNo"

RecordDayNo <- c(1:17)
StudyDayNo <- sprintf("Day.%02d", rep(1:17))

df.for.merge <- data.frame(RecordDayNo, StudyDayNo)

Totals.new <- left_join(Totals.final, df.for.merge)

Totals.new <- Totals.new %>% select(UserName, StudyDayNo, RecordDayNo, everything())

# Totals.new is the complete data for each day's nutritional intake 
# Needs to be matched with the fecal sample identifiers

fecal.mapping <- read.csv("~/Documents/Projects/MCTs/Data/source files for mapping/Fecal_Sample_Mapping.csv")

fecal.mapping <- fecal.mapping %>% filter(Replicate == 1)

Totals.new.with.mapping <- left_join(Totals.new, fecal.mapping)

Totals.new.with.mapping <- Totals.new.with.mapping %>% select(X.SampleID, UserName, StudyDayNo, RecordDayNo, everything())


# Add the mapping data to the Item.final file too

Items.new <- left_join(Items.final, df.for.merge)

Items.new <- Items.new %>% select(UserName, StudyDayNo, RecordDayNo, everything())



Items.new.with.mapping <- left_join(Items.new, fecal.mapping)

Items.new.with.mapping <- Items.new.with.mapping %>% select(X.SampleID, UserName, StudyDayNo, RecordDayNo, everything())

############### EXPORT THE TOTALS ###############################

write.table(Totals.new.with.mapping, 
            file = "/Users/abby/Documents/Projects/MCTs/Data/All Diet Data/Totals_to_use.txt", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


######################## EXPORT THE ITEMS FILE TO USE ######################
Items.new.with.mapping$IntakeStartDateTime <- as.character(Items.final$IntakeStartDateTime)
Items.new.with.mapping$IntakeEndDateTime <- as.character(Items.final$IntakeEndDateTime)

write.table(Items.new.with.mapping, 
            file = "/Users/abby/Documents/Projects/MCTs/Data/All Diet Data/Items_to_use.txt", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



