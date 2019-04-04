setwd("~/Documents/Projects/dietstudy_analyses/")

# plot predictions results
load(file = "data/LOOCV.prediction.results.Rdata")
load(file = "data/species.to.test.Rdata")
load(file = "data/test_personal_diet_dat.Rdata")
myresults <-test
species <- keepspecies

# format myresults as dataframes for plotting
# Note: should really be a function - but this will do.

mynames <- gsub("MCTs", "", names(dat))

myrown <- length(species)


own.diet <- data.frame(matrix(unlist(myresults$Master.own.diet), nrow = myrown, byrow = F))
rownames(own.diet) <- species
colnames(own.diet) <- mynames

others.diet <- data.frame(matrix(unlist(myresults$Master.others.diet), nrow = myrown, byrow = F))
rownames(others.diet) <- species
colnames(others.diet) <- mynames

own.microbiome <- data.frame(matrix(unlist(myresults$Master.own.microbiome), nrow = myrown, byrow = F))
rownames(own.microbiome) <- species
colnames(own.microbiome) <- mynames

scramble.diet <- data.frame(matrix(unlist(myresults$Master.scramble.diet), nrow = myrown, byrow = F))
rownames(scramble.diet) <- species
colnames(scramble.diet) <- mynames


# format for ggplot to make plots
require(ggplot2)
require(reshape2)

own.diet.melt <- rownames_to_column(own.diet, var = "Taxa")
own.diet.melt <- melt(own.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.corr")

others.diet.melt <- rownames_to_column(others.diet, var = "Taxa")
others.diet.melt <- melt(others.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "others.diet.corr")

own.microbiome.melt <-rownames_to_column(own.microbiome, var = "Taxa")
own.microbiome.melt <- melt(own.microbiome.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.mb.corr")

scramble.diet.melt <- rownames_to_column(scramble.diet, var = "Taxa")
scramble.diet.melt <- melt(scramble.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "scramble.diet.corr")


# make master df

master.df <- merge(own.diet.melt,own.microbiome.melt)
master.df <- merge(others.diet.melt, master.df)
master.df <- merge(scramble.diet.melt, master.df)

# master.df <- merge(own.diet.melt, scramble.diet.melt)
# master.df <- merge(others.diet.melt, master.df)

master.df.melt <- melt(master.df, id.vars = c("Taxa", "UserName"))

require(agricolae)
#Tukey
hsd=HSD.test(aov(value~variable,data=master.df.melt), "variable", group=T)

# get summary stats for plotting error bars
summary.master <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "mean")
colnames(summary.master) <- c("variable", "mean")
summary.sd <- aggregate(master.df.melt$value, by = list(master.df.melt$variable), FUN = "sd")
colnames(summary.sd) <- c("variable", "sd")

summary.master <- merge(summary.master, summary.sd)
summary.master$se <- summary.master$sd/sqrt(length(mynames)) # se is sd/sqrt(n)
summary.master <- merge(summary.master, as.data.frame(hsd$groups), by.x = "mean", by.y = "value")

summary.master

master.df.melt$variable <- factor(master.df.melt$variable, levels = c("own.diet.corr","own.mb.corr", "scramble.diet.corr", "others.diet.corr"))

# make the plot
mybarplot <- ggplot(data = master.df.melt, aes(x = variable, y = value)) +
  geom_bar(stat = "summary", fun.y = mean, aes(fill = as.factor(variable)), show.legend = F) +
  scale_x_discrete(labels = c("Diet +\nmicrobiome", "Microbiome\nonly", "Shuffled diet\n+ microbiome", "Microbiome\n+ another\nsubject's diet")) +
  scale_fill_manual(values = c("#5a2071", "#5f86b7","#bb7ad2","#64baaa")) +
  geom_errorbar(data = summary.master, aes(x = variable, ymin = mean - se, ymax = mean + se), width = 0.4, inherit.aes = F) +
  theme_classic() +
  ylab("Prediction accuracy\n(Pearson correlation of\nLOOCV species abundance)") +
  xlab("Model training set") +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text.x = element_text(size = 5)) 
  #geom_text(data = summary.master, aes(x = variable, y = mean + se, label=groups), vjust=-1, size = 3)

# pring the plot
mybarplot

ggsave("output/Figure4/Figure4C_LOOCV.pdf", mybarplot, height = 2.5, width = 2.5, device = "pdf")

