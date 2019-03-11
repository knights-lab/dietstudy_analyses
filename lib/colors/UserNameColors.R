# named colors for UserName 
# use these colors to be consistent throughout study plots

# map <- read.delim("~/Documents/Projects/dietstudy/data/maps/UserName_map.txt")
# map <- map[map$Study.Status == "Complete",]
# map <- droplevels(map)
# 
# map$UserName

UserNameColors <- c("#e6194b",
               "#3cb44b",
               "#ffe119",
               "#4363d8",
               "#f58231",
               "#911eb4",
               "#42d4f4",
               "#f032e6",
               "#bfef45",
               "#fabebe",
               "#469990",
               "#e6beff",
               "#9a6324",
               "#fffac8",
               "#800000",
               "#aaffc3",
               "#808000",
               "#ffd8b1",
               "#000075",
               "#a9a9a9",
               "#000000",
               "#ffffff",
               "#7460d2",
               "#ff0000",
               "#ec943c",
               "#faee69",
               "#0055fd",
               "#ff007e",
               "#4d4949",
               "#a6f1ff",
               "#2dff00",
               "#ff6300",
               "#965353",
               "#4600f8")

names(UserNameColors) <- c("MCTs03", "MCTs07", "MCTs22", "MCTs10", "MCTs13", "MCTs14", 
                           "MCTs15", "MCTs16", "MCTs18", "MCTs20", "MCTs21", "MCTs09", "MCTs23",
                           "MCTs24", "MCTs25", "MCTs26", "MCTs32", "MCTs33", "MCTs35",
                           "MCTs01", "MCTs04", "MCTs05", "MCTs06", "MCTs08", "MCTs11", 
                           "MCTs12", "MCTs19", "MCTs27", "MCTs28", "MCTs29", "MCTs31", "MCTs34",
                           "MCTs36", "MCTs37")
                           

UserNameColors


#Visualize colors
#pie(rep(1,length(UserNameColors)), labels = names(UserNameColors), col=UserNameColors)
