### Plot 3 parameter space with llhood as heat using rgl ###

## ll_smoothed.csv is a 5 column tab separated file
## par1, par2, par3, par4, llhood
## can plot 3 of the four pars at once
require(rgl)

ll<-read.csv('ll_smoothed.csv',sep='\t',header=FALSE)

##Trim the dataframe to exclude the wildly poor points
ll<-subset(ll,ll[,5]>-150)


colors<-c('blue','lightblue','yellow','red')
cus_col<-colorRampPalette(colors=colors, bias = 1, space = c("rgb", "Lab"),interpolate = c("linear", "spline"))



plot3d(ll[,1],ll[,2],ll[,3],col=cus_col(100)[floor(100*(ll[,5]-min(ll[,5]))/(max(ll[,5])-min(ll[,5])))],pch=20,type='l',lwd=5)

plot3d(ll[,1],ll[,2],ll[,4],col=cus_col(100)[floor(100*(ll[,5]-min(ll[,5]))/(max(ll[,5])-min(ll[,5])))],pch=20,type='l',lwd=5)

plot3d(ll[,1],ll[,3],ll[,4],col=cus_col(100)[floor(100*(ll[,5]-min(ll[,5]))/(max(ll[,5])-min(ll[,5])))],pch=20,type='l',lwd=5)

plot3d(ll[,2],ll[,3],ll[,4],col=cus_col(100)[floor(100*(ll[,5]-min(ll[,5]))/(max(ll[,5])-min(ll[,5])))],pch=20,type='l',lwd=5)


