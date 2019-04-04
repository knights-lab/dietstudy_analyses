require(tibble)


# map
map <- read.table("~/Documents/Projects/dietstudy/data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")


# Read in burst output table for greengenes
taxa <- read.delim(file = "~/Documents/Projects/dietstudy/data/greengenes/strain.txt", row=1, stringsAsFactors = F, sep ="\t")

# clean up the names
#### Manicure the samplenames, grab latest (mct study only!)
colnames(taxa) = gsub(".S[0-9]+.R1.001","",colnames(taxa))      # Clean old plate IDs
taxa = taxa[,order(colnames(taxa))]              # Sort nicely by sample ID
taxa = taxa[,-(grep("L0",colnames(taxa))-1)]     # Keep new runs only
colnames(taxa) = gsub(".S[0-9]+.L001.R1.001","",colnames(taxa)) # Clean new plate IDs

taxa = taxa[order(rowMeans(taxa),decreasing=T),] # Sort by avg. abundance

median(colSums(taxa))

# get names for subsetting from other data
keep <- colnames(read.delim(file = "~/Documents/Projects/dietstudy/data/processed_tax/taxonomy_clr_s.txt", row.names = 1))

# subset table
taxa <- taxa[,colnames(taxa) %in% keep]

minimal.processed <- taxa
minimal.processed <- rownames_to_column(minimal.processed, var = "#OTU")
# output minimally processed table
write.table(minimal.processed, file = "~/Documents/Projects/dietstudy/data/greengenes/strain.minimal.processed.txt", sep = "\t", quote = F, row.names = F)


# calculate relative abundance
taxa.ra <- sweep(taxa, 2, colSums(taxa), "/")

# drop rare taxa
taxa.sub <- taxa.ra[rowMeans(taxa.ra) >= 0.000001,]
median.depth <- median(colSums(taxa[rowMeans(taxa.ra) >= 0.000001,]))

taxa.norm.counts <- round(taxa.sub * median.depth)

taxa <- rownames_to_column(taxa.norm.counts, var = "#OTU")

# Write table to use with QIIME
write.table(taxa, file = "~/Documents/Projects/dietstudy/data/greengenes/strain.processed.txt", sep = "\t", quote = F, row.names = F)

# Now make a version that is collapsed by person
map <- read.delim("~/Documents/Projects/dietstudy/data/maps/SampleID_map.txt")
map <- map[map$X.SampleID %in% colnames(minimal.processed),]

# Now collapse by UserName (UN)
UN_taxa <- data.frame(matrix(nrow=nrow(minimal.processed), ncol=0))
rownames(UN_taxa) <- rownames(minimal.processed)

# collapse by username (sum raw counts across the rows)
for (subject in unique(map$UserName)) {
  submap <- map[map$UserName == subject,]
  subtax <- minimal.processed[,colnames(minimal.processed) %in% submap$X.SampleID]
  mytax <- as.data.frame(rowSums(subtax))
  colnames(mytax) <- subject
  UN_taxa <- cbind(UN_taxa, mytax)
}

UN_taxa <- round(UN_taxa)


# add back OTU names
rownames(UN_taxa) <- minimal.processed$`#OTU`
UN_taxa <- rownames_to_column(UN_taxa, var = "#OTU")

write.table(UN_taxa, file = "~/Documents/Projects/dietstudy/data/greengenes/UN_strain.processed.txt", sep = "\t", quote = F, row.names = F)
