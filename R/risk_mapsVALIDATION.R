### Risk and uncertainty maps VALIDATION data ###

## Read in lake data
lakes<-read.table("../2010_bytho_data/lakes_processed.csv")

## Read in 2EB ws boundary
library(maptools)
ws<-readShapePoly('../lakepolygons/2EB Boundary/Dissolved_Boundary')

## Custom legend
source('R/cust_colorlegend.R')


## eXample ##
if(F)
{
    plot(ws)
    points(lakes[,2],lakes[,3],col='grey')
    index_val_lakes <- which(lakes[,7]!=0)
    points(lakes[index_val_lakes_abs,2],lakes[index_val_lakes_abs,3],pch=20,col='red')
    index_val_lakes_abs <- which(lakes[,7]==1)
    points(lakes[index_val_lakes_abs,2],lakes[index_val_lakes_abs,3],pch=20,col='black')
}
#################



## Predictions ###
p_hat<-as.matrix(read.table('output/val_sim_props.tab'))

## Observed (validation)
o<-read.table('output/val_lakes.dat')
o<-o[,1]

## For plotting with colorRampPalettes
normalizeIntRange<-function(x,size=200)
{
 floor( ( x - min(x) ) / (max(x)-min(x)) * (size-1)) + 1
}
colors<-c('blue','yellow','red')
cus_col<-colorRampPalette(colors=colors, 
   bias = 1,  #spread out values more at the high end - for right skewed distributions.
   space = c("rgb", "Lab"),
   interpolate = c("linear", "spline")) 



par(mfrow=c(2,2),mar=c(0,0,0,0))
### MAP 1: Expected p_hat ### 
plot(ws,lwd=0.5)
points(lakes[,2],lakes[,3],col='grey',cex=0.4,pch=20)
expected_prob <- as.vector(apply(p_hat,2,mean))

# Normalize
n_ep <- normalizeIntRange(expected_prob)

points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
   col=cus_col(200)[n_ep],
    pch=20,
    cex=20*expected_prob)
points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
    cex=15*expected_prob)

colorlegend(posx=c(0.1,0.12),
    posy=c(0.6,0.9),
    col=cus_col(200),
    zlim=range(expected_prob),
    digit=1,
    main='Probability',
    zval=seq(0,max(expected_prob),0.1) )






### MAP 2: variance p_hat ###
plot(ws,lwd=0.5)
points(lakes[,2],lakes[,3],col='grey',cex=0.4,pch=20)
var_prob <- as.vector(apply(p_hat,2,var))

# Normalize
n_varp <- normalizeIntRange(var_prob)

points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
    col=cus_col(200)[n_varp],
    pch=20,
    cex=20*expected_prob)
points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
    cex=15*expected_prob)

colorlegend(posx=c(0.1,0.12),
    posy=c(0.6,0.9),
    col=cus_col(200),
    zlim=range(var_prob),
    digit=3,
    main="Variance",
    zval=seq(0,max(var_prob),0.001) )



### MAP 3: variance p_hat ###
if(F){
v <- vmapp(pred=p_hat,d=o)
delta <- as.vector(apply(v$delta,2,mean))
plot(expected_prob,delta)
}
plot(ws,lwd=0.5)
points(lakes[,2],lakes[,3],col='grey',cex=0.4,pch=20)

# Normalize
n_delta <- normalizeIntRange(delta)

points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
    col=cus_col(200)[n_delta],
    pch=20,
    cex=20*expected_prob)
points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
    cex=15*expected_prob)

colorlegend(posx=c(0.1,0.12),
    posy=c(0.6,0.9),
    col=cus_col(200),
    zlim=range(delta),
    digit=3,
    main=expression(delta),
    zval=seq(min(delta),max(delta),0.005) )


################################################################################
################################################################################


if(F){
## Which of the validation lakes were sampled previously:
last_obs_abs <- lakes[index_val_lakes,6]

## Plot expected prob vs variance 
## Color code according to which lakes had been previously sampled.
x11()
plot(expected_prob,
    var_prob,
    col=(last_obs_abs==0)+1,
    pch=(last_obs_abs==0)+1)

legend("topleft",
    col=1:2,
    pch=1:2,
    legend=c('Sampled previously','Not sampled previously'))


## Plot histograms of the bootstrapped predictions
## Color code by year sampled.
pdf('plots/boot_pred_hists.pdf')
    par(mfrow=c(3,3))
    for(i in 1:length(index_val_lakes))
    {
        hist(p_hat[,i],
            col=(last_obs_abs[i]==0)+1,
            xlim=c(0,0.5),
            probability=TRUE,
            main=paste('Lake',index_val_lakes[i],'Sampled in',last_obs_abs[i]))
    }
dev.off()
}
