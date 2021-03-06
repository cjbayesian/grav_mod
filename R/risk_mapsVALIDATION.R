### Risk and uncertainty maps VALIDATION data ###
library(rvmapp)
library(maptools)
library(raster) ## For image import (brick)
library(png) ## simpler image import.
#library(maps) ## For scale bar in map (doesn't work since map is in UTM)

## Read in lake data
lakes<-read.table("../2010_bytho_data/lakes_processed.csv")
#lakes<-read.table("sims/simmed_lakes.csv")

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
p_hat<-as.matrix(read.table('output/val_sim_propsA.tab'))
#p_hat<-as.matrix(read.table('sims/gb_output/val_sim_props.tab'))

## Observed (validation)
o<-read.table('output/val_lakes.dat')
#o<-read.table('sims/gb_output/val_lakes.dat')
o<-o[,1]


### Ordered preds ###
var_prob <- as.vector(apply(p_hat,2,var))
expected_prob <- as.vector(apply(p_hat,2,mean))
ci <- apply(p_hat,2,quantile, probs=c(0.025,0.975))

ord <- order(expected_prob)
par(cex=1.5)
plot(expected_prob[ord],col=(o==1)[ord]+1,pch=20,ylim=c(0,0.8),xlab='Rank',ylab='Predicted Probability')
for(i in 1:length(ord))
{
    arrows(i,expected_prob[ord[i]],i,ci[1,ord[i]],col=(o==1)[ord[i]]+1,angle=90, length = 0.05)
    arrows(i,expected_prob[ord[i]],i,ci[2,ord[i]],col=(o==1)[ord[i]]+1,angle=90, length = 0.05)
    if((o==1)[ord[i]])
        arrows(i,ci[2,ord[i]] + 0.2 ,i,ci[2,ord[i]] + 0.1,col=(o==1)[ord[i]]+1,angle=45, length = 0.05,lwd=3)
}
legend('topleft',legend=c('Uninvaded','Invaded'),pch=20,lty=1,col=1:2,cex=1.5)
###


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


pdf('plots/risk_maps.pdf',width=16,height=16)
    par(mfrow=c(2,2),mar=c(0,0,0,0))

    ### MAP 1: variance p_hat ###
    plot(ws,lwd=0.5)
    points(lakes[,2],lakes[,3],col='grey',cex=0.4,pch=20)
    var_prob <- as.vector(apply(p_hat,2,var))
    expected_prob <- as.vector(apply(p_hat,2,mean))
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

    loc_panel_letter <- loc_prop(c(0.05,0.95))
    text(loc_panel_letter[1],loc_panel_letter[2],labels="A)",cex=2)
    ## ---- ##





    ### MAP 2: Delta ###
    if(T){
        v <- vmapp(pred=p_hat,d=o)
        delta <- as.vector(apply(v$delta,2,mean))
        #plot(expected_prob,delta,ylim=c(-0.1,0))
        #abline(0,-1,lty=2)
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

    text(loc_panel_letter[1],loc_panel_letter[2],labels="B)",cex=2)
    ## ---- ##

    
    ### Map Inset for context of 2eb in Ontario ###
    #ont<-read.table('../survey/analysis/ont_outline.dat')
    #ontimg<-brick('../ontstub.jpg')  ## Stub -- Will replace with GIS map
    par(mar=c(3,3,0,2))
    plot.new()
    lim <- par()
    ima <- readPNG('../ontario_map2eb_noalpha_sm.png')
    rasterImage(ima,lim$usr[1], lim$usr[3], lim$usr[2], lim$usr[4])
    #par(mar=c(4,3,3,1))
    #plot(ont,type='l',xlab='Longitude',ylab='Latitude')
    #plot(ws,lwd=0.5,add=T)
    ## ---- ##

    
    ### VMAPP Delta plot ###
    par(mar=c(5,4,4,2),cex=1.5,lwd=3)
    plot(v,main="VMAPP estimation\nof prediction bias")
    title(main="C)",adj=0)
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
