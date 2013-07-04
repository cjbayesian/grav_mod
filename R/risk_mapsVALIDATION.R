### Risk and uncertainty maps VALIDATION data ###
library(rvmapp)
library(maptools)


## Read in lake data
lakes<-read.table("../2010_bytho_data/lakes_processed.csv")

## Read in 2EB ws boundary
ws<-readShapePoly('../lakepolygons/2EB Boundary/Dissolved_Boundary')

## Custom legend
source('R/cust_colorlegend.R')


    index_val_lakes <- which(lakes[,7]!=0)
    index_val_lakes_abs <- which(lakes[,7]==1)

## eXample ##
if(F)
{
    plot(ws)
    points(lakes[,2],lakes[,3],col='grey')
    points(lakes[index_val_lakes_abs,2],lakes[index_val_lakes_abs,3],pch=20,col='red')
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

## Location calculation for legends
loc_prop <- function(pos=c(0.5,0.5))
{
    pars <- par("usr")
    dx <- pars[2] - pars[1]
    posx <- pars[1] + pos[1] * dx

    dy <- pars[4] - pars[3]
    posy <- pars[3] + pos[2] * dy

    return(c(posx,posy))
}


pdf('plots/risk_maps.pdf',width=16,height=8)
    par(mfrow=c(1,2),mar=c(0,0,0,0))

    ### MAP 1: variance p_hat ###
    plot(ws,lwd=0.5)
    points(lakes[,2],lakes[,3],col='grey',cex=0.4,pch=20)
    var_prob <- as.vector(apply(p_hat,2,var))

    colors<-c('blue','yellow','red')
    cus_col<-colorRampPalette(colors=colors, 
       bias = 1,  #spread out values more at the high end - for right skewed distributions.
       space = c("rgb", "Lab"),
       interpolate = c("linear", "spline")) 

    # Normalize
    n_varp <- normalizeIntRange(var_prob)

    points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
        col=cus_col(200)[n_varp],
        pch=20,
        cex=20*expected_prob)
    ## Add circles around each point
    points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
        cex=15*expected_prob)
    ## Add crosses at invaded lakes
    points(lakes[index_val_lakes[o==1],2],lakes[index_val_lakes[o==1],3],
        cex=10*expected_prob[o==1],pch=4)

    ## Legends ##
    p_levels <- seq(0.05,0.5,0.15)
    y_locs <- seq(0.6,0.85,length=4)
    loc <- loc_prop(c(0.25,0.91))
    text(loc[1],
        loc[2],
        labels='Expected\nProbability')
    for(i in 1:4)
    {
        loc <- loc_prop(c(0.25,y_locs[i]))
        points(loc[1],
            loc[2],
            cex=15*p_levels[i])
        loc <- loc_prop(c(0.3,y_locs[i]))
        text(loc[1],
            loc[2],
            labels=p_levels[i])
    }

    colorlegend(posx=c(0.1,0.12),
        posy=c(0.6,0.9),
        col=cus_col(200),
        zlim=range(var_prob),
        digit=3,
        main="Variance",
        zval=seq(0,max(var_prob),0.001) )
    ## ---- ##





    ### MAP 2: Delta ###
    if(F){
        v <- vmapp(pred=p_hat,d=o)
        delta <- as.vector(apply(v$delta,2,mean))
        plot(expected_prob,delta,ylim=c(-0.1,0))
        abline(0,-1,lty=2)
    }
    plot(ws,lwd=0.5)
    points(lakes[,2],lakes[,3],col='grey',cex=0.4,pch=20)

    colors<-c('Purple','blue','white')
        cus_col<-colorRampPalette(colors=colors, 
        bias = 1,  #spread out values more at the high end - for right skewed distributions.
        space = c("rgb", "Lab"),
       interpolate = c("linear", "spline")) 

    # Normalize
    n_delta <- normalizeIntRange(delta)

    points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
        col=cus_col(200)[n_delta],
        pch=20,
        cex=20*expected_prob)
    points(lakes[index_val_lakes,2],lakes[index_val_lakes,3],
        cex=15*expected_prob)
    ## Add crosses at invaded lakes
    points(lakes[index_val_lakes[o==1],2],lakes[index_val_lakes[o==1],3],
        cex=10*expected_prob[o==1],pch=4)

    ## Legends ##
    p_levels <- seq(0.05,0.5,0.15)
    y_locs <- seq(0.6,0.85,length=4)
    loc <- loc_prop(c(0.25,0.91))
    text(loc[1],
        loc[2],
        labels='Expected\nProbability')
    for(i in 1:4)
    {
        loc <- loc_prop(c(0.25,y_locs[i]))
        points(loc[1],
            loc[2],
            cex=15*p_levels[i])
        loc <- loc_prop(c(0.3,y_locs[i]))
        text(loc[1],
            loc[2],
            labels=p_levels[i])
    }

    colorlegend(posx=c(0.1,0.12),
        posy=c(0.6,0.9),
        col=cus_col(200),
        zlim=range(delta),
        digit=3,
        main=expression(delta),
        zval=seq(min(delta),max(delta),0.005) )
    ## ---- ##

dev.off()
################################################################################
################################################################################



## Which of the validation lakes were sampled previously:
last_obs_abs <- lakes[index_val_lakes,6]

## Plot expected prob vs variance 
## Color code according to which lakes had been previously sampled.
x11()
point_symbols <- numeric(length(expected_prob))

point_symbols[!o & last_obs_abs==0] <- 1 # Not inv. Not sampled
point_symbols[!o & last_obs_abs!=0] <- 2 # Not inv. Prev sampled
point_symbols[o & last_obs_abs==0] <- 16 # Inv2010, Not sampled
point_symbols[o & last_obs_abs!=0] <- 17 # Inv2010, Prev sampled


plot(expected_prob,
    var_prob,
    col=(last_obs_abs!=0)+1,
    pch=point_symbols,
    xlab='E[Predicted Probability]',
    ylab='Var[Predicted Probability]' )

legend("topleft",
    col=c(1,2,'white'),
    pch=1:3,
    legend=c('Not Sampled previously',
        'Sampled previously',
        '\n (Filled symbols indicate Invaded Sites)'),
    bty='n')

if(F){
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
