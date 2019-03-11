# load raw counts
run1 <- read.delim("~/Documents/Projects/dietstudy/data/seq_counts/human_reads_removed/wcrun1.txt")
run2 <- read.delim("~/Documents/Projects/dietstudy/data/seq_counts/human_reads_removed/wcrun2.txt")

rownames(run1) <- gsub(".fastq", "", run1$Sample)
rownames(run1) <- gsub("_", ".", rownames(run1))
rownames(run2) <- gsub(".fastq", "", run2$Sample)
rownames(run2) <- gsub("_", ".", rownames(run2))
rownames(run2) <- gsub("-", ".", rownames(run2))

bothruns <- rbind(run1,run2)
bothruns <- as.data.frame(t(bothruns), stringsAsFactors = FALSE)

# manicure to just keep the samples that we acutally use for downstream analysis
colnames(bothruns) = gsub(".S[0-9]+.R1.001","",colnames(bothruns));      # Clean old plate IDs
bothruns = bothruns[,order(colnames(bothruns))];              # Sort nicely by sample ID
bothruns = bothruns[,-(grep("L0",colnames(bothruns))-1)];     # Keep new runs only
colnames(bothruns) = gsub(".S[0-9]+.L001.R1.001","",colnames(bothruns)); # Clean new plate IDs

bothruns <- as.data.frame(t(bothruns),stringsAsFactors = FALSE)
bothruns$seq.count <- as.numeric(bothruns$seq.count)

# drop the blanks
bothruns <- bothruns[-(grep("Blank", rownames(bothruns))),]

# calculate the mean + sd for each sample
mean(bothruns$seq.count)
sd(bothruns$seq.count)

# calculate the mean + sd for each participant
map <- read.table("~/Documents/Projects/dietstudy/data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")
map <- map[map$X.SampleID %in% rownames(bothruns),]

perpersonsums <- NULL

# sum reads per person
for (i in unique(map$UserName)){
  subgroup <- map[map$UserName == i,]
  perperson <- bothruns[rownames(bothruns) %in% subgroup$X.SampleID,]
  tmp <- as.data.frame(sum(perperson$seq.count, na.rm = T))
  rownames(tmp) <- i
  perpersonsums <- rbind(perpersonsums,tmp)
}

colnames(perpersonsums) <- "sum"

mean(perpersonsums$sum)
sd(perpersonsums$sum)
