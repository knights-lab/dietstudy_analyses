require(zoo)

set.seed(42)

#make the fake dataset that is arranged how my real data look
mytestdat <- data.frame(SampleID = rep(1,17), 
                        MbSampDayNo = seq(1:17), 
                        DietFromDayNo = c(0:16), 
                        Fakeyfake = sample(1:17))

mytestdat


#function for use with rollapply
my_fun = function(x) {d = (length(x):1)-1; f = 2^-d; sum(f*x)}

# lets break down the function in parts:
# input 3 numbers that make this make sense
# want 100% of 0, 50% of 1, and 25% of 2, should sum to 1
testfakes <- c(2, 1, 0)


# 1. calculate the day offset index d
testd <- (length(testfakes):1)-1
testd

# 2. calculate the decay weight f
testf <- 2^-testd
testf

# 3. sum the numbers in the correct proportion
sum(testf*testfakes)

# confirm the function does the same thing when it's all done together
my_fun(testfakes)

# cool that works, so how to apply it repeatedly with a sliding window across the entire vector?

# Enter rollapply from the package zoo
#make the weighted average
mytestdat$weightedave <- rollapply(mytestdat$Fakeyfake, # pass your vector or dataframe
                                   width = 17,# set the width of the sliding window
                                   FUN = my_fun, # your function to apply
                                   by.column=T, # if T applied to each column seperately, redundant here, but good for scaling up
                                   align="right", # I think this controls how the output results are named
                                   partial = T) # allows for partial lengths of the rolling window - perfect for this example

mytestdat[,c(4,5)]
