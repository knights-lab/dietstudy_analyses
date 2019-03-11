homedir <- "/Users/abby/Documents/Projects/dietstudy/"

require(robCompositions)
require(tibble)

setwd(dir = homedir)


# load map
tax_map <- read.table("data/maps/taxonomy_norm_map.txt", sep = "\t", header = T, comment = "")

# Load cleaned pre-processed normalized (relative abundance) values at the species level
tax_norm_s <- read.delim("data/processed_tax/taxonomy_norm_s.txt", row = 1, stringsAsFactors = F)

tax <- tax_norm_s
# Now collapse by UserName (UN)
UN_tax <- data.frame(matrix(nrow=nrow(tax), ncol=0))
rownames(UN_tax) <- rownames(tax)

# collapse by username (average CLR distances)
for (i in unique(tax_map$UserName)) {
  submap <- tax_map[tax_map$UserName == i,]
  subtax <- tax[,colnames(tax) %in% submap$X.SampleID]
  mytax <- as.data.frame(rowMeans(as.matrix(subtax)))
  colnames(mytax) <- i
  UN_tax <- cbind(UN_tax, mytax)
}

# transpose for imptation
transposed_UN_tax <- t(UN_tax)

myimpR = impRZilr(transposed_UN_tax, maxit = 3, method = "lm", dl = rep(0.000001,ncol(transposed_UN_tax)), verbose = T) # pick a good detection limit...

UN_tax_clr_s = t(cenLR(myimpR$x)$x.clr)

colnames(UN_tax_clr_s) <- colnames(UN_tax)


# prep for export
UN_tax_clr_s <- rownames_to_column(UN_tax, var = "#taxonomy")

# export median CLR table
write.table(UN_tax_clr_s, file="data/processed_UN_tax/UN_tax_CLR_mean_norm_s.txt", quote=F, sep="\t", row.names = F, col.names = T)
# 
# 
# comptax <- read.delim("data/processed_UN_tax/UN_taxonomy_clr_s.txt", row = 1, header = T)
# 
# require(ape)
# test1 <- pcoa(dist(t(UN_tax_clr_s)))$vectors[,1:2]
# test2 <- pcoa(dist(t(comptax)))$vectors[,1:2]
# 
# plot(test1)
# plot(test2)
