

W<-c(5,10)
i=1
curve(x^2/sum(W^2),xlim=c(0,20),col=i)
for(e in seq(2,-2,length.out=10))
{
   curve(x^e/sum(W^e),xlim=c(0,20),col=i,add=TRUE)
   i<-i+1
}





############################################

pp<-read.table('output/pp.dat')

pp_mat<-array(dim=c(21,1646))
for(lake in 1:1646)
{
   for(year in 1990:2010)
   {
      index<-pp[,1]==year & pp[,2]==lake
      if(sum(index)==1)
         pp_mat[year-1989,lake]<-pp[index,3]
      else
          pp_mat[year-1989,lake]<-NA
   }
}

x11()
to<-max(pp_mat,na.rm=TRUE)
plot(pp_mat[,1],type='l',ylim=c(0,to),xlim=c(1,21))
for(i in 2:1646)
{
   lines(pp_mat[,i],lty=i)
}


max_pp<-apply(pp_mat,2,max,na.rm=TRUE)
max_pp 

x11()
par(mfrow=c(2,1))
plot(max_pp,lakes[,4]) ##known invs

inv<-read.csv('output/inv_stat.dat',header=FALSE)
plot(max_pp,inv[,1]) ##after simulation







##############################################
library(rgl)

ll<-read.table('output/traf_ll.dat')

# static plot
#library(scatterplot3d)
#scatterplot3d(ll[, 1], ll[, 2], ll[, 3])

 ll<-ll[ll[,3]>(-1000),]

plot3d(x=ll[, 1], y=ll[, 2], z=ll[, 3],
   col=heat.colors(1000)[ floor( (ll[,3]-min(ll[,3])) / ( max(ll[,3])-min(ll[,3]) ) * (500-1))+1],
   xlab='e',
   ylab='c',
   zlab='Log L')

plot(x=ll[, 1], y=ll[, 2],
   col=heat.colors(1000)[ floor( (ll[,3]-min(ll[,3])) / ( max(ll[,3])-min(ll[,3]) ) * (500-1))+1] )






