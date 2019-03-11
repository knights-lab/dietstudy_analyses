
setwd(dir = "/Users/abby/Documents/Projects/dietstudy_analyses/")
require(ape)
require(tidyverse)

load(file = "data/test_personal_diet_dat.Rdata")
load(file = "data/species.to.test.Rdata")
################### TEST DIET IMPACT ###########

# 
# 
# # want to look at only same speices across people
# # so find the species in everyonen - can probably do this more eligantly later if needed 
# all_species <- unlist(most.abund)
# common_species <- names(which(table(all_species)==28)) # uses common species found in everyone
# 

species <- species.to.test



"test.diet.impact" <- function(){
  set.seed(42)
  
  Master.own.diet <- NULL
  Master.others.diet <- NULL
  Master.own.microbiome <- NULL
  Master.adj.own.diet.pvalues <- NULL
  
  
  for (k in 1:length(dat)) {
    # person A index
    a <- k
    
    # person B index
    possible.b <- 1:length(dat)
    possible.b <- possible.b[possible.b != a]
    b <- sample(possible.b,1)
    
    # get the species names
    today.species <- paste0(species, ".x")
    tomorrow.species <- paste0(species, ".y")
    
    # create empty list for storing values
    own.diet.predictions <- NULL
    others.diet.predictions <- NULL
    own.microbiome.predictions <- NULL
    own.diet.pvalues <- NULL
    
    for (i in seq_along(tomorrow.species)) {
      
      tomorrow <- tomorrow.species[i]
      today <- today.species[i]
      
      mb.features <- c(today, "Mb.Axis.1","Mb.Axis.2", "Mb.Axis.3")
      dt.features <- c("Dt.Axis.1","Dt.Axis.2", "Dt.Axis.3")
      
      # predict person A's microbiome from microbiome alone
      # regress tomorrow (yA) from today's species and microbiome
      # predict tomorrow (yAhat)
      yA.feats <- dat[[a]][,mb.features]
      yA <- dat[[a]][,tomorrow]
      yA.model <- lm(yA ~., yA.feats)
      yAhat <- predict(yA.model)
      
      # get residuals
      yAr <- yAhat - yA
      
      # regress the residuals on diet
      modelA.feats <- dat[[a]][,dt.features]
      modelA <- lm(yAr ~., modelA.feats)
      yrAhat <- predict(modelA)

      # subtract predicted residuals (modelA) from yhat
      # and correlate with y
      own.diet <- cor(yAhat - yrAhat, yA)
     
      # compare yA and yAhat to 
      own.microbiome <- cor(yAhat, yA)
      
      # delta between own diet and own mb
      own.diet.delta <- own.diet - own.microbiome
      
      ##### permutation part #####
      # predict for A from scrambled diet features
      # scramble the diet features
      scramble.cor <- NULL
      scramble.cor.delta <- NULL
      
      for (j in 1:1000) { # TODO dont forget to change back
        scramble.Afeats <- modelA.feats[sample(nrow(modelA.feats)),]
        yrAhat.scramble <- predict(modelA, scramble.Afeats)
        scramble.cor[[j]] <- cor(yAhat - yrAhat.scramble, yA)
        scramble.cor.delta[[j]] <- scramble.cor[[j]] - own.microbiome
      } 

      own.diet.p <- mean(scramble.cor.delta > own.diet.delta)
      
      ##
      # predict person B
      # predict person B's species from B microbiome alone
      yB.feats <- dat[[b]][,mb.features]
      yB <- dat[[b]][,tomorrow]
      yB.model <- lm(yB ~., yB.feats)
      yBhat <- predict(yB.model)
      # get residuals
      yBr <- yBhat - yB
      
      #run predict on the model with person A, but passing in the diet data for person B 
      #(you can pass in a new data matrix to predict) 
      #to get their externally-predicted residuals yr.B.hat.
      modelB.feats <- dat[[b]][,dt.features]
      yrBhat <- predict(modelA, modelB.feats)

      #subtract yr.B.hat from yhat.B and correlate the final result with y.B.
      others.diet <- cor(yBhat - yrBhat, yB)
      
     
      
      # store variables
      
      own.diet.predictions[[tomorrow]] <- own.diet
      own.diet.pvalues[[tomorrow]] <- own.diet.p
      others.diet.predictions[[tomorrow]] <- others.diet
      own.microbiome.predictions[[tomorrow]] <- own.microbiome

    }
    Master.own.diet[[k]] <- own.diet.predictions
    Master.others.diet[[k]] <- others.diet.predictions
    Master.own.microbiome[[k]] <- own.microbiome.predictions
    Master.adj.own.diet.pvalues[[k]] <- p.adjust(own.diet.pvalues, "fdr")  

  }
  
  return(list(Master.own.diet = Master.own.diet, 
              Master.others.diet = Master.others.diet, 
              Master.own.microbiome = Master.own.microbiome, 
              Master.adj.own.diet.pvalues = Master.adj.own.diet.pvalues))
  
}


#myresults <- test.diet.impact()
save(myresults, file ='data/diet.responsive.bugs.Rdata')


# format myresults as dataframes for plotting
own.diet <- data.frame(matrix(unlist(myresults$Master.own.diet), nrow = length(species), byrow = F))
rownames(own.diet) <- species

mynames <- gsub("MCTs", "", names(dat))
colnames(own.diet) <- mynames

others.diet <- data.frame(matrix(unlist(myresults$Master.others.diet), nrow = length(species), byrow = F))
rownames(others.diet) <- species
colnames(others.diet) <- mynames

own.microbiome <- data.frame(matrix(unlist(myresults$Master.own.microbiome), nrow = length(species), byrow = F))
rownames(own.microbiome) <- species
colnames(own.microbiome) <- mynames

own.diet.pvalues <- data.frame(matrix(unlist(myresults$Master.adj.own.diet.pvalues), nrow = length(species), byrow = F))
rownames(own.diet.pvalues) <- species
colnames(own.diet.pvalues) <- mynames





# format for ggplot to make plots
require(ggplot2)
require(reshape2)

own.diet.melt <- rownames_to_column(own.diet, var = "Taxa")
own.diet.melt <- melt(own.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.corr")

others.diet.melt <- rownames_to_column(others.diet, var = "Taxa")
others.diet.melt <- melt(others.diet.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "others.diet.corr")

own.microbiome.melt <-rownames_to_column(own.microbiome, var = "Taxa")
own.microbiome.melt <- melt(own.microbiome.melt, id.vars = "Taxa", variable.name = "UserName", value.name = "own.mb.corr")

own.diet.pvalues <- rownames_to_column(own.diet.pvalues, var = "Taxa")
own.diet.pvalues.melt <- melt(own.diet.pvalues, id.vars = "Taxa", variable.name = "UserName", value.name = "own.diet.p")



# make master df

master.df <- merge(own.diet.melt,others.diet.melt)
master.df <- merge(own.microbiome.melt, master.df)

master.df.melt <- melt(master.df, id.vars = c("Taxa", "UserName"))


# make paired difference plots by person

require(gridExtra)
require(directlabels)

mb.dt <- master.df.melt[master.df.melt$variable %in% c("own.mb.corr", "own.diet.corr"),]

mb.dt.wide <- merge(master.df, own.diet.pvalues.melt)
mb.dt.wide <- mb.dt.wide[,!colnames(mb.dt.wide) %in% c("others.diet.corr")]
mb.dt.wide$sig<- ifelse(mb.dt.wide$own.diet.p < 0.05, "sig", "ns")
mb.dt.wide$groups <- ifelse(mb.dt.wide$sig == "sig", mb.dt.wide$Taxa, NA)

mb.dt.for.merge <- mb.dt.wide[colnames(mb.dt.wide) %in% c("Taxa", "UserName", "groups")]

mb.dt <- merge(mb.dt, mb.dt.for.merge)


# make multiple plots in a loop

single_plots <- as.list(NULL)

plot <- mb.dt[order(mb.dt$UserName, decreasing = F),]
plot <- plot[order(plot$groups, na.last = F),]

group.names <- unique(mb.dt$groups)
group.colors <- c("#858c69", "#3d7600", "#69ef7b", 
                  "#02531d", "#12982d", "#ffede7",
                  "#c25945", "#ffaae8", "#ff0000",
                  "#ffa787", "#b6601e", "#e242a8",
                  "#8f0000", "#dd6c43", "#b04fbc",
                  "#b04fbc", "#4b03a9", "#ec9bfa",
                  "#a717d9", "#5bc07f", "#0024ff",
                  "#8e4fbc", "#aaa4e1", "#afd35a",
                  "#6f1253", "#3c3e00", "yellow",
                  "#1f72d4", "#01ceff", "#71619c",
                  "#aaa4e1")
names(group.colors) <- group.names



for (i in 1:length(unique(plot$UserName))) {
  j <- as.character(unique(plot$UserName)[[i]])
  plot_un <- plot[plot$UserName == j,]
  
  plot_un <- plot_un[order(plot_un$groups, na.last = T),]
  
  single_plots[[i]] <- ggplot(data = plot_un, aes(x = variable, y = value)) + 
    facet_free(.~UserName, scales = "free") +
    geom_line(aes(group = Taxa), color = "grey", alpha = 0.5) + 
    geom_line(aes(group = Taxa, color = factor(groups)), size = 1) +
    scale_color_manual(values = group.colors) +
    #geom_dl(aes(label = groups), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.5, hjust = 1)) +
    geom_dl(aes(label = groups, color = groups), method = "last.polygons") +
    geom_point() +
    theme_classic() +
    ylim(0,1) +
    scale_x_discrete(expand = c(0,0.1), labels = c("Microbiome\nonly", "Diet and\nmicrobiome")) +
    expand_limits(x = 2.5) +
    theme(legend.position = "none",
          legend.text = element_text(size = 7),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = 4)) 
}

# get find out what people have significant diet-related species
keep_plots <- unique(subset(plot, !is.na(plot$groups))$UserName)
#keep_plots <- keep_plots[!keep_plots %in% c("11", "12")] # don't plot soylent people
loc <- which(mynames %in% keep_plots)


multiplot <- arrangeGrob(grobs = single_plots[loc], ncol = 4)

ggsave(plot = multiplot, filename = "output/multiplot.pdf", width = 7.5, height = 10 )

