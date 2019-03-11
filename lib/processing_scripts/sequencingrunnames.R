# List of sample names to combine each person's sequencing for CAZyme analysis

map <- read.delim(file = "~/Documents/Projects/dietstudy/data/maps/SampleID_map.txt")
map$X.SampleID <- gsub("\\.", "_", map$X.SampleID)

run1names <- read.delim(file = "~/Documents/Projects/dietstudy/data/run1names.txt", header = F, col.names = "Run1_name")
run2names <- read.delim(file = "~/Documents/Projects/dietstudy/data/run2names.txt", header = F, col.names = "Run2_name")


# add a column that contains just the first part of each samples name
run1names$X.SampleID <- run1names$Run1_name
run1names$X.SampleID <- gsub("_S.*", "", run1names$X.SampleID)


run2names$X.SampleID <- run2names$Run2_name
run2names$X.SampleID <- gsub("_S.*", "", run2names$X.SampleID)
run2names$X.SampleID <- gsub("-", "_", run2names$X.SampleID)


require(dplyr)

all <- full_join(run1names,run2names)
all <- full_join(all, map)
all$Run1_name <- as.character(all$Run1_name)
all$Run2_name <- as.character(all$Run2_name)


# drop lines that contain the word Blank

all <- all[grep("Blank", all$X.SampleID, invert = T),]

# drop the NA's
all <- all[!is.na(all$Run1_name),]

all$Seq_name_used <- ifelse(is.na(all$Run2_name)==TRUE, all$Run1_name, all$Run2_name)

length(all$Seq_name_used)

write.table(all$Seq_name_used, file = "Documents/Projects/dietstudy/data/maps/seq_name_used.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# TODO: export all the names per person for CAZyme analysis
