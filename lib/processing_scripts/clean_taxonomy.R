homedir <- "/Users/abby/Documents/Projects/dietstudy/"

require(robCompositions)
require(tibble)

setwd(dir = homedir)

# load map
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# load taxonomy file from October 2, 2017
tax <- read.delim("raw/hoho.tax.txt", row = 1, stringsAsFactors = F)

# colsums tax
tax_check <- as.data.frame(sort(colSums(tax)))

# Drop samples with low read counts (for original table this was >=20000)
# Keep everything after the bad high-read blank
tax <- tax[, colSums(tax) >=23500]

# Drop dropouts 
dropouts =  as.character(map[map$UserName %in% c("MCTs02", "MCTs17", "MCTs30"),"X.SampleID"])
tax <- tax[, !(colnames(tax) %in% dropouts)]

# see how many non-bacterial hits
tax_non_bact <- tax[-grep("k__Bacteria", rownames(tax)),]

# only a few, but a couple have decent counts in some people, so they may survive filtering and I will keep them in.

# DONT DO ANYMORE: remove non-bacteria taxa
# tax <- tax[grep("k__Bacteria", rownames(tax)),]

# make a limited map that matches all samples left
tax_map <- map[map$X.SampleID %in% colnames(tax),]

# Drop rare taxa that don't appear in at least 25% of people in the study
bugs_per_person <- list()

for (i in unique(tax_map$UserName)) {
  submap <- tax_map[tax_map$UserName == i,]
  subset <- tax[, colnames(tax) %in% submap$X.SampleID]
  subset <- subset[rowSums(subset > 0) > ncol(subset)/4, ] # these are the bugs present in this person at least 1/4 days in the study
  mybugs <- rownames(subset)
  bugs_per_person[[i]] <- rownames(subset)
}

n25 <- round(length(bugs_per_person) * 0.25) # 10% of people is about 8
nonuni <- unlist(bugs_per_person)
counts <- table(nonuni)
counts <- counts[counts > n25]
keep <- as.vector(names(counts))

tax <- tax[rownames(tax) %in% keep,]  # limit to taxa we care about
tax <- tax[rowSums(tax) > 0,]         # double check all 0 sum taxa are gone


# write this version to a file for use in other scripts
write.table(tax, file="data/processed_tax/taxonomy_preprocessed.txt", quote=F, sep="\t", row.names = T, col.names = T)


# Remove vary rare taxa from the tax table
  #normalize
tax_norm <- sweep(tax, 2, colSums(tax), "/")
  #drop low abundance
tax_norm_dla <- tax_norm[rowMeans(tax_norm) >= 0.0001,] # keeps 290 species-level anotations

# now limit tax to this
tax <- tax[rownames(tax) %in% rownames(tax_norm_dla),]

# Summarizing at different levels
split <- strsplit(rownames(tax),";")             # Split and rejoin on lv7 to get species level

# Species
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_s <- rowsum(tax,taxaStrings)
rownames(tax_s) = sapply(strsplit(rownames(tax_s),";"),function(x) paste(x[1:7],collapse=";"));

# Genus
taxaStrings <- sapply(split,function(x) paste(x[1:6],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_g <- rowsum(tax,taxaStrings) 
rownames(tax_g) = sapply(strsplit(rownames(tax_g),";"),function(x) paste(x[1:6],collapse=";"));

# Family
taxaStrings <- sapply(split,function(x) paste(x[1:5],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_f <- rowsum(tax,taxaStrings)  
rownames(tax_f) = sapply(strsplit(rownames(tax_f),";"),function(x) paste(x[1:5],collapse=";"));

# Order
taxaStrings <- sapply(split,function(x) paste(x[1:4],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_o <- rowsum(tax,taxaStrings) 
rownames(tax_o) = sapply(strsplit(rownames(tax_o),";"),function(x) paste(x[1:4],collapse=";"));

# Class
taxaStrings <- sapply(split,function(x) paste(x[1:3],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_c <- rowsum(tax,taxaStrings) 
rownames(tax_c) = sapply(strsplit(rownames(tax_c),";"),function(x) paste(x[1:3],collapse=";"));

# Phylum
taxaStrings <- sapply(split,function(x) paste(x[1:2],collapse=";"))
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
tax_p <- rowsum(tax,taxaStrings) 
rownames(tax_p) = sapply(strsplit(rownames(tax_p),";"),function(x) paste(x[1:2],collapse=";"));


# create a normalized taxa tables at each level to use with some downstream analysis
##### Make relative abundance table for plotting bar graphs ####
# Normalize
tax_norm_s <- sweep(tax_s, 2, colSums(tax_s), "/")
tax_norm_g <- sweep(tax_g, 2, colSums(tax_g), "/")
tax_norm_f <- sweep(tax_f, 2, colSums(tax_f), "/")
tax_norm_o <- sweep(tax_o, 2, colSums(tax_o), "/")
tax_norm_c <- sweep(tax_c, 2, colSums(tax_c), "/")
tax_norm_p <- sweep(tax_p, 2, colSums(tax_p), "/")

# Sort by average abundance
tax_norm_s <- tax_norm_s[order(rowMeans(tax_norm_s),decreasing=T),]
tax_norm_g <- tax_norm_g[order(rowMeans(tax_norm_g),decreasing=T),]
tax_norm_f <- tax_norm_f[order(rowMeans(tax_norm_f),decreasing=T),]
tax_norm_o <- tax_norm_o[order(rowMeans(tax_norm_o),decreasing=T),]
tax_norm_c <- tax_norm_c[order(rowMeans(tax_norm_c),decreasing=T),]
tax_norm_p <- tax_norm_p[order(rowMeans(tax_norm_p),decreasing=T),]


### Make clr-adjusted taxa table
# 1. start with table of counts without bad samples or low abundance taxa
# to limit to just the taxa that remain the the tax_norm table
# columns are samples
tax_clr_s <- tax_s[rownames(tax_s) %in% rownames(tax_norm_s),]
tax_clr_g <- tax_g[rownames(tax_g) %in% rownames(tax_norm_g),]
tax_clr_f <- tax_f[rownames(tax_f) %in% rownames(tax_norm_f),]
tax_clr_o <- tax_o[rownames(tax_o) %in% rownames(tax_norm_o),]
tax_clr_c <- tax_c[rownames(tax_c) %in% rownames(tax_norm_c),]
tax_clr_p <- tax_p[rownames(tax_p) %in% rownames(tax_norm_p),]


# 2. interpolate 0s
mynames <- colnames(tax)    # store colnames (sample names) for later

# change int > num
tax_clr_s <- tax_clr_s*1.0      
tax_clr_g <- tax_clr_g*1.0  
tax_clr_f <- tax_clr_f*1.0  
tax_clr_o <- tax_clr_o*1.0  
tax_clr_c <- tax_clr_c*1.0  
tax_clr_p <- tax_clr_p*1.0 

# transpose for imputation
# columns are taxa
transposed_s <- t(tax_clr_s)
transposed_g <- t(tax_clr_g)
transposed_f <- t(tax_clr_f)
transposed_o <- t(tax_clr_o)
transposed_c <- t(tax_clr_c)
transposed_p <- t(tax_clr_p)

#this is the imputation step and it takes ages to finish, so just run when needed
myimpR_s = impRZilr(transposed_s, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_s)), verbose = T)
myimpR_g = impRZilr(transposed_g, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_g)), verbose = T)
myimpR_f = impRZilr(transposed_f, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_f)), verbose = T)
myimpR_o = impRZilr(transposed_o, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_o)), verbose = T)
myimpR_c = impRZilr(transposed_c, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_c)), verbose = T)
myimpR_p = impRZilr(transposed_p, maxit = 3, method = "lm", dl = rep(2,ncol(transposed_p)), verbose = T)

# save as rdata to use again later without re-running imputation
save(myimpR_s, file = "data/processed_tax/myimpR_s.rdata")
save(myimpR_g, file = "data/processed_tax/myimpR_g.rdata")
save(myimpR_f, file = "data/processed_tax/myimpR_f.rdata")
save(myimpR_o, file = "data/processed_tax/myimpR_o.rdata")
save(myimpR_c, file = "data/processed_tax/myimpR_c.rdata")
save(myimpR_p, file = "data/processed_tax/myimpR_p.rdata")

# load files for use
load(file = "data/processed_tax/myimpR_s.rdata")
load(file = "data/processed_tax/myimpR_g.rdata")
load(file = "data/processed_tax/myimpR_f.rdata")
load(file = "data/processed_tax/myimpR_o.rdata")
load(file = "data/processed_tax/myimpR_c.rdata")
load(file = "data/processed_tax/myimpR_p.rdata")

# 3. clr tranformation
tax_clr_s = t(cenLR(myimpR_s$x)$x.clr)  # clr tranformation of imputed table
tax_clr_g = t(cenLR(myimpR_g$x)$x.clr)
tax_clr_f = t(cenLR(myimpR_f$x)$x.clr)
tax_clr_o = t(cenLR(myimpR_o$x)$x.clr)
tax_clr_c = t(cenLR(myimpR_c$x)$x.clr)
tax_clr_p = t(cenLR(myimpR_p$x)$x.clr)

# re-add names
colnames(tax_clr_s) <- mynames
colnames(tax_clr_g) <- mynames
colnames(tax_clr_f) <- mynames
colnames(tax_clr_o) <- mynames
colnames(tax_clr_c) <- mynames
colnames(tax_clr_p) <- mynames

# 4. now good to go for export/use with stats
tax_clr_s <- as.data.frame(tax_clr_s)
tax_clr_g <- as.data.frame(tax_clr_g)
tax_clr_f <- as.data.frame(tax_clr_f)
tax_clr_o <- as.data.frame(tax_clr_o)
tax_clr_c <- as.data.frame(tax_clr_c)
tax_clr_p <- as.data.frame(tax_clr_p)

med_depth <- median(colSums(tax_s))

#### Make counts table for alpha diversity from normalized table #####
# multiply re-normalized reduced table by a factor to remove decimals and round
tax_counts_s <- round(sweep(tax_norm_s, 2, colSums(tax_norm_s),'/')*med_depth)
tax_counts_g <- round(sweep(tax_norm_g, 2, colSums(tax_norm_s),'/')*med_depth)
tax_counts_f <- round(sweep(tax_norm_f, 2, colSums(tax_norm_s),'/')*med_depth)
tax_counts_o <- round(sweep(tax_norm_o, 2, colSums(tax_norm_s),'/')*med_depth)
tax_counts_c <- round(sweep(tax_norm_c, 2, colSums(tax_norm_s),'/')*med_depth)
tax_counts_p <- round(sweep(tax_norm_p, 2, colSums(tax_norm_s),'/')*med_depth)

# prep for export
tax_norm_s <- rownames_to_column(tax_norm_s, var = "#taxonomy")
tax_norm_g <- rownames_to_column(tax_norm_g, var = "#taxonomy")
tax_norm_f <- rownames_to_column(tax_norm_f, var = "#taxonomy")
tax_norm_o <- rownames_to_column(tax_norm_o, var = "#taxonomy")
tax_norm_c <- rownames_to_column(tax_norm_c, var = "#taxonomy")
tax_norm_p <- rownames_to_column(tax_norm_p, var = "#taxonomy")

tax_counts_s <- rownames_to_column(tax_counts_s, var = "#taxonomy")
tax_counts_g <- rownames_to_column(tax_counts_g, var = "#taxonomy")
tax_counts_f <- rownames_to_column(tax_counts_f, var = "#taxonomy")
tax_counts_o <- rownames_to_column(tax_counts_o, var = "#taxonomy")
tax_counts_c <- rownames_to_column(tax_counts_c, var = "#taxonomy")
tax_counts_p <- rownames_to_column(tax_counts_p, var = "#taxonomy")

# fix names on taxonomy map
colnames(tax_map)[1] <- "#SampleID"

# prep clr tables for export
tax_clr_s <- rownames_to_column(tax_clr_s, var = "#taxonomy")
tax_clr_g <- rownames_to_column(tax_clr_g, var = "#taxonomy")
tax_clr_f <- rownames_to_column(tax_clr_f, var = "#taxonomy")
tax_clr_o <- rownames_to_column(tax_clr_o, var = "#taxonomy")
tax_clr_c <- rownames_to_column(tax_clr_c, var = "#taxonomy")
tax_clr_p <- rownames_to_column(tax_clr_p, var = "#taxonomy")

# export cleaned map
write.table(tax_map, file="data/maps/taxonomy_norm_map.txt", quote = F, sep="\t", row.names = F, col.names = T)

# export normalized tables
write.table(tax_norm_s, file="data/processed_tax/taxonomy_norm_s.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(tax_norm_g, file="data/processed_tax/taxonomy_norm_g.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(tax_norm_f, file="data/processed_tax/taxonomy_norm_f.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(tax_norm_o, file="data/processed_tax/taxonomy_norm_o.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(tax_norm_c, file="data/processed_tax/taxonomy_norm_c.txt", quote=F, sep="\t", row.names = F, col.names = T)
write.table(tax_norm_p, file="data/processed_tax/taxonomy_norm_p.txt", quote=F, sep="\t", row.names = F, col.names = T)

# export counts tables (really probably only need counts at the species level)
write.table(tax_counts_s, file = "data/processed_tax/taxonomy_counts_s.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(tax_counts_g, file = "data/processed_tax/taxonomy_counts_g.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(tax_counts_f, file = "data/processed_tax/taxonomy_counts_f.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(tax_counts_o, file = "data/processed_tax/taxonomy_counts_o.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(tax_counts_c, file = "data/processed_tax/taxonomy_counts_c.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(tax_counts_p, file = "data/processed_tax/taxonomy_counts_p.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# export clr transformed taxonomy tables for each level
write.table(tax_clr_s, file = "data/processed_tax/taxonomy_clr_s.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(tax_clr_g, file = "data/processed_tax/taxonomy_clr_g.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(tax_clr_f, file = "data/processed_tax/taxonomy_clr_f.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(tax_clr_o, file = "data/processed_tax/taxonomy_clr_o.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(tax_clr_c, file = "data/processed_tax/taxonomy_clr_c.txt", quote = F, sep="\t", row.names = F, col.names = T)
write.table(tax_clr_p, file = "data/processed_tax/taxonomy_clr_p.txt", quote = F, sep="\t", row.names = F, col.names = T)
