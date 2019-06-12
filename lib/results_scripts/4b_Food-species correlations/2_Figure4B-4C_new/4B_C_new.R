# make some plots to show the unique individualized relationships
# Highlight these two relationships, for starters


# Darkgreen_vegetables Bacteroides sp. 1_1_6 0.81 0.002 0.1800 20
# Darkgreen_vegetables Bacteroides sp. 1_1_6 -0.73 < 0.001 0.1100 18
# Darkgreen_vegetables Bacteroides sp. 1_1_6 0.72 0.009 0.1800 35

# load data for these three subjects
# read in individual summarized food tables for each person and combine them
setwd("/Users/abby/Documents/Projects/dietstudy_analyses/data/procrustes/data_username_decay/")
source("../../../lib/colors/UserNameColors.R")
names(UserNameColors) <- gsub("MCTs", "", names(UserNameColors))


# read all food files into a list
temp <- list.files(pattern = "*_food.txt")
foods_3d <- lapply(temp, function(x) read.delim(x, row.names = "taxonomy"))
foods_3d <- lapply(foods_3d, function(x) x[!names(x) == "X.FoodID"])
names(foods_3d) <- gsub("_food.txt", "", temp)
split = strsplit(rownames(foods_3d[[1]]),";")
foodStrings = sapply(split,function(x) paste(x[1:2],collapse=";"))
foods_3d <- lapply(foods_3d, function(x) rowsum(x,foodStrings))

# read in tax files into a list
temp <- list.files(pattern = "*_tax.txt")
tax_3d <- lapply(temp, function(x) read.delim(x, row.names = 1))
names(tax_3d) <- gsub("_tax.txt", "", temp)


# limit to the people we are interested in and their foods
dkveg <- foods_3d[c("MCTs20", "MCTs18", "MCTs35")]
dkveg <- lapply(dkveg, function(x) x[rownames(x) == "L1_Vegetables;L2_Darkgreen_vegetables",])

bact <- tax_3d[c("MCTs20", "MCTs18", "MCTs35")]
bact <- lapply(bact, function(x) x[(grep("1_1_6", rownames(x))),])

dkveg <- unlist(dkveg)
bact <- unlist(bact)

plot = data.frame(darkveg = dkveg, bacteroides = bact)
plot$Subject <- substr(rownames(plot), 5,6)
plot$ID <- substr(rownames(plot), 8, length(rownames(plot)))

require(ggplot2)

myplot <- ggplot(data = plot, aes(y = bacteroides, x = darkveg)) + 
  geom_point(aes(color = Subject), size = 2) +
  scale_color_manual(values = UserNameColors) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  facet_wrap(.~Subject, scales = "free") +
  theme_classic() +
  ylab("Bacteroides sp. 1 1 6\n(CLR adjusted abundance)") +
  xlab("Vegetables: Darkgreen vegetables\n(dry weight, g)") +
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 3), 
        legend.position = "none",
        strip.text = element_text(size = 6),
        panel.grid.major = element_line(colour = "lightgrey"))

myplot

ggsave(filename = "../../../output/Figure4/Figure4B_new.pdf", height = 1.7, width = 2.5)

# Pork Bacteroides oleiciplenus 0.87 < 0.001 0.0220 33 
# Pork Bacteroides oleiciplenus -0.73 < 0.001 0.1300 23
# Pork Bacteroides oleiciplenus -0.77 0.001 0.2000 22


# limit to the people we are interested in and their foods
pork <- foods_3d[c("MCTs33", "MCTs23", "MCTs22")]
pork <- lapply(pork, function(x) x[(grep("Pork", rownames(x))),])

ole <- tax_3d[c("MCTs33", "MCTs23", "MCTs22")]
ole <- lapply(ole, function(x) x[(grep("oleiciplenus", rownames(x))),])

pork <- unlist(pork)
ole <- unlist(ole)

plot = data.frame(Pork = pork, "Bacteroides oleiciplenus" = ole)
plot$Subject <- substr(rownames(plot), 5, 6)
plot$ID <- substr(rownames(plot), 8, length(rownames(plot)))

require(ggplot2)

myplot <- ggplot(data = plot, aes(y = `Bacteroides.oleiciplenus`, x = Pork)) + 
  geom_point(aes(color = Subject), size = 2) +
  scale_color_manual(values = UserNameColors) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  facet_wrap(.~Subject, scales = "free") +
  theme_classic() +
  ylab("Bacteroides oleiciplenus\n(CLR adjusted abundance)") +
  xlab("Meats: Pork\n(dry weight, g)") +
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 3), 
        legend.position = "none",
        strip.text = element_text(size = 6),
        panel.grid.major = element_line(colour = "lightgrey"))
myplot

ggsave(filename = "../../../output/Figure4/Figure4C_new.pdf", height = 1.7, width = 2.5)



#Cakes_cookies_pies_pastries_bars	Roseburia inulinivorans	0.74	0.006	0.16	35
#Cakes_cookies_pies_pastries_bars	Roseburia inulinivorans	0.73	0.001	0.12	07
#Cakes_cookies_pies_pastries_bars	Roseburia inulinivorans	0.71	0.001	0.13	18


# limit to the people we are interested in and their foods
cakes <- foods_3d[c("MCTs35", "MCTs07", "MCTs18")]
cakes <- lapply(cakes, function(x) x[(grep("Cakes_cookies", rownames(x))),])

rose <- tax_3d[c("MCTs35", "MCTs07", "MCTs18")]
rose <- lapply(rose, function(x) x[(grep("Roseburia inulinivorans", rownames(x))),])

cakes <- unlist(cakes)
rose <- unlist(rose)

plot = data.frame(Cakes = cakes, "Roseburia inulinicorans" = rose)
plot$Subject <- substr(rownames(plot), 5, 6)
plot$ID <- substr(rownames(plot), 8, length(rownames(plot)))

require(ggplot2)

myplot <- ggplot(data = plot, aes(y = `Roseburia.inulinicorans`, x = Cakes)) + 
  geom_point(aes(color = Subject), size = 2) +
  scale_color_manual(values = UserNameColors) +
  geom_smooth(method = "lm", se = F, color = "black", size = 0.5) +
  facet_wrap(.~Subject, scales = "free") +
  theme_classic() +
  ylab("Roseburia inulinivorans\n(CLR adjusted abundance)") +
  xlab("Grains: Cakes, cookies, bars, etc.\n(dry weight, g)") +
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 3), 
        legend.position = "none",
        strip.text = element_text(size = 6),
        panel.grid.major = element_line(colour = "lightgrey"))
myplot

ggsave(filename = "../../../output/Figure4/Figure4D_new.pdf", height = 1.7, width = 2.5)

