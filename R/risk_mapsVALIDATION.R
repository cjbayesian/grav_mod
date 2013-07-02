### Risk and uncertainty maps VALIDATION data ###

## Read in lake data
lakes<-read.table("../2010_bytho_data/lakes_processed.csv")

## Read in 2EB ws boundary
library(maptools)
ws<-readShapePoly('../lakepolygons/2EB Boundary/Dissolved_Boundary')

## eXample ##
plot(ws)
points(lakes[,2],lakes[,3],col='grey')
index_val_lakes <- which(lakes[,7]!=0)
points(lakes[index_val_lakes_abs,2],lakes[index_val_lakes_abs,3],pch=20,col='red')
index_val_lakes_abs <- which(lakes[,7]==1)
points(lakes[index_val_lakes_abs,2],lakes[index_val_lakes_abs,3],pch=20,col='black')
#################



## Predictions ###
p_hat<-as.matrix(read.table('output/val_sim_props.tab'))

## Observed (validation)
o<-read.table('output/val_lakes.dat')
o<-o[,1]


## For plotting with colorRampPalettes
normalizeIntRange<-function(x,size=200)
{
 floor( ( x - min(x) ) / (max(x)-min(x)) * size)
}
colors<-c('blue','yellow','red')
cus_col<-colorRampPalette(colors=colors, 
   bias = 1,  #spread out values more at the high end - for right skewed distributions.
   space = c("rgb", "Lab"),
   interpolate = c("linear", "spline")) 



par(mfrow=c(2,2),mar=c(0,0,0,0))
### MAP 1: Expected p_hat ### 
plot(ws,lwd=0.5)
points(lakes[,2],lakes[,3],col='grey',cex=0.5,pch=20)
expected_prob<-apply(p_hat,2,mean)

# Normalize
n_ep <- normalizeIntRange(expected_prob)

points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
   col=cus_col(200)[n_ep],
    pch=20,
    cex=2)



### MAP 2: variance p_hat ###
plot(ws,lwd=0.5)
points(lakes[,2],lakes[,3],col='grey',cex=0.5,pch=20)
var_prob<-apply(p_hat,2,var)

# Normalize
n_varp <- normalizeIntRange(var_prob)

points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
   col=cus_col(200)[n_varp],
    pch=20,
    cex=1.5)


