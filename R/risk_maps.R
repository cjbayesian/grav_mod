########## Generate Risk maps ########

## This script builds 3 maps.
##    1. A map of environmental suitability
#         a) mean
#         b) variance
##    2. A map of invasion risk




# Read in mcmc posterior
#trace<-read.table("output/lib.mcmc")
trace<-array(rnorm(1000*18),dim=c(1000,18)) ##standin sim - remove when lib.mcmc complete.

# Read in lake data
lakes<-read.table("../2010_bytho_data/lakes_processed.csv")


n_lakes<-nrow(lakes)
post_means<-apply(trace,2,mean)

z<-sapply(1:n_lakes,function(i){ 
      post_means[5]+ sum( post_means[6:18] * lakes[i,7:19] )}) 

alpha = -log(1-(1/(1+exp(-z))))

# Normalize
nalpha<- floor( ( alpha - min(alpha) ) / (max(alpha)-min(alpha)) * 200)


colors<-c('blue','yellow','red')
cus_col<-colorRampPalette(colors=colors, 
   bias = 5,  #spread out values more at the high end - for right skewed distributions.
   space = c("rgb", "Lab"),
   interpolate = c("linear", "spline")) 

plot(lakes[,2],lakes[,3],
   pch=20,
   cex=log(lakes[,1]),
   col=cus_col(200)[nalpha])



### MAP 1: Env suitability ###




### MAP 2: Invasion Risk ###
## Pull in sim data from c++ --TODO code sim from posterior

