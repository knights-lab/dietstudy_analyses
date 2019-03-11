require(tibble)

# set wd
setwd(dir = "Documents/Projects/dietstudy/")

# map
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")


# Read in burst tables
#taxa <- read.delim(file = "raw/mcTax.tsv", row=1, stringsAsFactors = F, sep ="\t")
taxa <- read.delim(file = "raw/hoho.tax", row=1, stringsAsFactors = F, sep ="\t")

# clean up the names
#### Manicure the samplenames, grab latest (mct study only!)
colnames(taxa) = gsub(".S[0-9]+.R1.001","",colnames(taxa))      # Clean old plate IDs
taxa = taxa[,order(colnames(taxa))]              # Sort nicely by sample ID
taxa = taxa[,-(grep("L0",colnames(taxa))-1)]     # Keep new runs only
colnames(taxa) = gsub(".S[0-9]+.L001.R1.001","",colnames(taxa)) # Clean new plate IDs

taxa = taxa[order(rowMeans(taxa),decreasing=T),] # Sort by avg. abundance

# get names for subsetting
keep <- names(which(colSums(taxa) >= 23000))

# add #taxonomy label
taxa <- rownames_to_column(taxa, var = "#taxonomy")

# Write table to use in other analyses
# write.table(taxa, file = "raw/taxonomy_burst_omegaGWG.txt", sep = "\t", quote = F, row.names = F)
write.table(taxa, file = "raw/hoho.tax.txt", sep = "\t", quote = F, row.names = F)

# do the same for capitalist strain table
#strain <- read.delim(file = "raw/mcStrain.tsv", row=1, stringsAsFactors = F, sep = "\t")
strain <- read.delim(file = "raw/hoho.otu.str", row = 1, stringsAsFactors = F, sep = "\t")

med_samp_depth <- median(colSums(strain))

# clean up the names
#### Manicure the samplenames, grab latest (mct study only!)
colnames(strain) = gsub(".S[0-9]+.R1.001","",colnames(strain))      # Clean old plate IDs
strain = strain[,order(colnames(strain))]              # Sort nicely by sample ID
strain = strain[,-(grep("L0",colnames(strain))-1)]     # Keep new runs only
colnames(strain) = gsub(".S[0-9]+.L001.R1.001","",colnames(strain)) # Clean new plate IDs

# drop blanks from strain table
strain <- strain[,-(grep("Blank", colnames(strain)))]

# drop samples that have low counts
strain <- strain[,colnames(strain) %in% keep]

# and drop dropouts
dropouts =  as.character(map[map$UserName %in% c("MCTs02", "MCTs17", "MCTs30"),"X.SampleID"])
strain <- strain[, !colnames(strain) %in% dropouts]


#### Massage for use #####
strain.n <- sweep(strain, 2, colSums(strain), "/")
#strain = strain[order(rowMeans(strain),decreasing=T),] # Sort by avg. abundance
strain <- strain[rowMeans(strain.n)>= 0.001,]
strain.n <- strain.n[rowMeans(strain.n)>=0.001,]
strain <- strain[rowSums(strain.n > 0) > 20,]
strain.n <- strain.n[rowSums(strain.n > 0) > 20,]
strain.n <- round(strain.n * med_samp_depth)


# fix names for export (and so the match the tree for qiime)
strain <- rownames_to_column(strain, var = "#strain")
strain$`#strain` <- gsub("_", " ", strain$`#strain`)

strain.n <- rownames_to_column(strain.n, var = "#strain")
strain.n$`#strain` <- gsub("_", " ", strain.n$`#strain`)



# Write table to use in other analyses
write.table(strain, file = "raw/hoho.otu.txt", sep = "\t", quote = F, row.names = F)
write.table(strain.n, file = "raw/hoho.otu_norm.txt", sep = "\t", quote = F, row.names = F)
