require(zoo)

set.seed(42)

#make the fake dataset
mytestdat <- data.frame(SampleID = rep(1,17), 
                        MbSampDayNo = seq(1:17), 
                        DietFromDayNo = c(0:16), 
                        Fakeyfake = sample(1:17))

#function for use with rollapply
my_fun = function(x) {d = (length(x):1)-1; f = 2^-d; sum(f*x)}

#make the weighted average
mytestdat$weightedave <- rollapply(mytestdat$Fakeyfake, 
                                   width = 17, 
                                   FUN = my_fun, 
                                   by.column=T, 
                                   align="right", 
                                   partial = T)
