homedir <- "/Users/abby/Documents/Projects/dietstudy/"

require(tibble)
require(robCompositions)


setwd(dir = homedir)

# load map
map <- read.table("data/maps/SampleID_map.txt", sep = "\t", header = T, comment = "")

# Load KEGG module file from APRIL 2018
func<- read.delim("data/KEGG/Both_runs/kegg_module_table.txt", row = 1, stringsAsFactors = F)

# limit to correct run samples
# Manicure the func samplenames, grab latest
colnames(func) = gsub(".S[0-9]+.R1.001","",colnames(func));      # Clean old plate IDs
func = func[,order(colnames(func))];              # Sort nicely by sample ID
func = func[,-(grep("L0",colnames(func))-1)];     # Keep new runs only
colnames(func) = gsub(".S[0-9]+.L001.R1.001","",colnames(func)); # Clean new plate IDs


# Drop dropouts from func
dropouts =  as.character(map[map$UserName %in% c("MCTs02", "MCTs17", "MCTs30"),"X.SampleID"])
func <- func[, !(colnames(func) %in% dropouts)]

# export for plotting
write.table(func, file = "data/processed_KEGG/clean_both_KEGG_modules_table.txt", quote = F, sep = "\t", row.names = T, col.names = NA)

# make a limited map/func df that matches all samples left (517 total once blanks are removed)
func_map <- map[map$X.SampleID %in% colnames(func),]
func <- func[,colnames(func) %in% map$X.SampleID]

# find minium value in each sample (trick to first get rid of all zeros)
func_fix <- func
func_fix[func_fix == 0] <- NA

func_fix_mins <- apply(func_fix, 2, function(x) min(x, na.rm=T)) 

# divide all counts in that sample by that minimum value
func_fix_scale <- func_fix/func_fix_mins
func_fix_scale[is.na(func_fix_scale)] <- 0

func_fix_clr <- func_fix_scale

# Do CLR (multiplicitive count replacement code from Gabe)
func_fix_clr = t(func_fix_clr); eps = 0.5
func_fix_clr = func_fix_clr*(1 - rowSums(func_fix_clr==0)*eps/rowSums(func_fix_clr))
func_fix_clr[func_fix_clr==0]=eps
func_fix_clr = sweep(func_fix_clr,1,rowSums(func_fix_clr),'/');
ls = log(func_fix_clr)
func_fix_clr = t(ls - rowMeans(ls))
func_fix_clr = func_fix_clr[,!is.nan(colSums(func_fix_clr))]


# filter out low coeficient of variation (SD/mean) 
# want standard deviation pretty high relative to the mean 
func_filt <- func_fix_clr # columns are samples
func_means <- rowMeans(func_filt)
func_sd <- apply(func_filt, 1, function(x) sd(x))
func_var <- func_sd/abs(func_means)  # range 0.04 to 211.11, mean = 1.2, median = 0.32
quantile(func_var)

# find the highest variance functions
tail(sort(func_var),n = 20)

# drop rows with low variance
func_filt <- func_filt[func_var >= 3.6,] # keeps 20 highest variance functions after CLR
func <- func[rownames(func) %in% rownames(func_filt),]


# write this version to a file for use in other scripts
write.table(func, file="data/processed_KEGG/preprocessed_KEGG_modules.txt_filtered", quote=F, sep="\t", row.names = T, col.names = NA)
write.table(func_filt, file="data/processed_KEGG/preprocessed_KEGG_modules_filtered_clr.txt", quote=F, sep="\t", row.names = T, col.names = NA)
write.table(func_fix_clr, file = "data/processed_KEGG/preprocessed_KEGG_modules_all_clr.txt", quote = F, sep = "\t", row.names = T, col.names = NA)

### create average table per person ###
# set up the summary dataframes
func_ave <- data.frame(matrix(nrow=nrow(func_fix_clr), ncol=0))

# add rownames to each (rows are KEGG modules)
rownames(func_ave) <- rownames(func_fix_clr)

# average per person
for (i in unique(func_map$UserName)){
  subgroup <- func_map[func_map$UserName == i,]
  funcsub <- func_fix_clr[,colnames(func_fix_clr) %in% subgroup$X.SampleID]
  tmp <- as.data.frame(rowMeans(funcsub, na.rm = T))                    # need to choose if means or sums is better for analysis
  colnames(tmp) <- i
  func_ave <- cbind(func_ave,tmp)
}

# make no soylent version
func_ave_no_soy <- func_ave[,!colnames(func_ave) %in% c("MCTs11", "MCTs12")]

# export files for use
write.table(func_ave, file = "data/processed_KEGG/ave_KEGG_modules.txt", quote = F, sep = "\t", row.names = T, col.names = NA)
write.table(func_ave_no_soy, file = "data/processed_KEGG/ave_KEGG_modules_no_soy.txt", quote = F, sep = "\t", row.names = T, col.names = NA)
