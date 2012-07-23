########## Generate Risk maps ########

## This script builds 3 maps.
##    1. A map of environmental suitability
#         a) mean
#         b) variance
##    2. A map of invasion risk

datetime<-system("date +%y%m%d%H%M",intern=TRUE)


# Read in mcmc posterior
trace<-read.table("output/lib.mcmc")
#trace<-array(rnorm(1000*18),dim=c(1000,18)) ##standin sim - remove when lib.mcmc complete.


trace<-trace[2000:nrow(trace),] ##burn-in

# Read in lake data
lakes<-read.table("../2010_bytho_data/lakes_processed.csv")


n_lakes<-nrow(lakes)
post_means<-apply(trace,2,mean)

z<-sapply(1:n_lakes,function(i){ 
      post_means[5] + sum( post_means[6:18] * lakes[i,7:19] )}) 

mean_alpha = -log(1-(1/(1+exp(-z))))

mean_alpha<-log(mean_alpha)

hist(mean_alpha,breaks=50)

# Normalize
nalpha<- floor( ( mean_alpha - min(mean_alpha) ) / (max(mean_alpha)-min(mean_alpha)) * 200)

### MAP 1: Env suitability ###

colors<-c('blue','green','red')
cus_col<-colorRampPalette(colors=colors, 
   bias = 1,  #spread out values more at the high end - for right skewed distributions.
   space = c("rgb", "Lab"),
   interpolate = c("linear", "spline")) 

plot(lakes[,2],lakes[,3],
   pch=20,
   #cex=log(lakes[,1]),
   cex=1.5,
   col=cus_col(200)[nalpha])

source('~/R_bits/cust_colorlegend.R')
colorlegend(posx=c(0.8,0.85),posy=c(0.02,0.35),col=cus_col(100)[25:100],zlim=c(0,1),digit=2) #zval=seq(0,1,0.25) )

   ### VARIANCE ###
N<-1000
library(parallel)
sub_trace<-sample(1:nrow(trace),N)
post_z<-mcmapply(sub_trace,mc.cores = 2,function(j){   
   z<-sapply(1:n_lakes,function(i){ 
      trace[j,5] + sum( trace[j,6:18] * lakes[i,7:19] )}) 
   print(j)
   return(z)
   })

dim(post_z)

post_alpha = -log(1-(1/(1+exp(-z))))

post_var<-apply(post_z,1,var)

plot(lakes[,2],lakes[,3],
   pch=20,
   cex=post_var
   #cex=1.5,
   #col=cus_col(200)[nalpha]
   )



### MAP 2: Invasion Risk ###
## Pull in sim data from c++ --TODO code sim from posterior

