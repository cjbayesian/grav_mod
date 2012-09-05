

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

 ll<-ll[ll[,4]>(-1000),]
 ll<-ll[ll[,2]>0,]

ht_col<-heat.colors(1000)[ floor( (ll[,4]-min(ll[,4])) / ( max(ll[,4])-min(ll[,4]) ) * (800-1))+1]

plot3d(x=ll[, 1], y=ll[, 2], z=ll[, 4],
   col=ht_col,
   xlab='e',
   ylab='c',
   zlab='Log L',
   cex=2)
 
writeWebGL(width=500,height=500)

#### Experimental: x=e,y=c,time=alpha ####
alphas<-unique(ll[,3])
mll_alpha<-numeric(length(alphas))
mle_e<-numeric(length(alphas))
mle_c<-numeric(length(alphas))
i<-1
for(alpha in alphas)
{
   tmpll<-ll[ll[,3]==alpha,]

   ht_col<-heat.colors(1000)[ floor( (tmpll[,4]-min(tmpll[,4])) / ( max(tmpll[,4])-min(tmpll[,4]) ) * (800-1))+1]

   png(paste('plots/',1+alpha,'.png',sep=''))
   plot(x=tmpll[, 1], y=tmpll[, 2],
      col=ht_col,
      xlab='e',
      ylab='c',
      pch=15,
      cex=3,
      main=paste('alpha = ',alpha,'\nMax LL = ',max(tmpll[,4])))

   max_e<-tmpll[which.max(tmpll[,4]),1]
   max_c<-tmpll[which.max(tmpll[,4]),2]
   points(max_e,max_c)
   abline(v=max_e,lty=2)
   abline(h=max_c,lty=2)
   dev.off()

   plot(x=tmpll[, 1], y=tmpll[, 2],
      col=ht_col,
      xlab='e',
      ylab='c',
      pch=15,
      cex=3,
      main=paste('alpha = ',alpha,'\nMax LL = ',max(tmpll[,4])))

   
   mle_e[i]<-max_e
   mle_c[i]<-max_c
   mll_alpha[i]<-max(tmpll[,4])
   i<-i+1
}

## Is there a ridge? (expect flat if there is.)
par(mfrow=c(3,1))
plot(alphas,mll_alpha)
plot(alphas,mle_e,type='b')
plot(alphas,mle_c,type='b')








########################### Examine the stochastic spread component #########################

## Rows are across a dimension of parameter space.
n_inv<-as.matrix(read.table('output/n_inv_file.dat'))
ll<-as.matrix(read.table("output/ll.dat"))
  
  colors<-c('blue','lightblue','yellow','red')
  cus_col<-colorRampPalette(colors=colors, bias = 1, space = c("rgb", "Lab"),interpolate = c("linear", "spline"))
  ht_col<-cus_col(1000)[ floor( (ll[,2]-min(ll[,2])) / ( max(ll[,2])-min(ll[,2]) ) * (1000-1))+1]

for(i in 1:ncol(n_inv))
{
   plot(ll[,1],n_inv[,i],
      col=ht_col,
      ylim=c(0,max(n_inv)),
      xlab='e_par',
      ylab='# of sites invaded',
      main=i+1989)
   colorlegend(posx=c(0.8,0.84),
      posy=c(0.5,0.9),
      col=cus_col(1000),
      zlim=c(min(ll[,2]),max(ll[,2])),
      digit=2,
      zval=seq(floor(min(ll[,2])),floor(max(ll[,2])),length.out=5),
      main='Log L')
}

for(i in 1:ncol(n_inv))
{
   plot(n_inv[,i],ll[,2],main=i+1989)

}


e_levels<-unique(ll[,1])
for(j in 1:length(e_levels))
{
   plot(n_inv[es ,],ylim=c(0,max(n_inv)))
   for(es in which(ll[,1]==e_levels[j]))
      lines(n_inv[es ,],col=j)
}


for(i in 1:ncol(n_inv))
   plot(ll[,1],,main=i+1989)



plot(ll[,1],ll[,2],
   xlab='e',
   ylab='Log L')



########################################################################################
################### Colorfull spread vis ###############################################

d<-as.matrix(read.table('output/t_mcmc.dat'))
lakes<-as.matrix(read.csv('../2010_bytho_data/lakes_processed.csv',sep="\t",header=FALSE))

invaded_mat<-array(0,dim=c(ncol(d),max(d)-min(d)+1)) ##n_lakes x years


mm<-mclapply(1:ncol(d),function(i){
   tmp_vec<-NULL
   for(year in min(d):max(d))
      tmp_vec[(year-min(d))+1]<-sum(d[,i]==year)
   print(i)
   return(tmp_vec)
})

mm<-do.call(cbind, mm)



  colors<-c('white','blue','lightblue','yellow','red')
  cus_col<-colorRampPalette(colors=colors, bias = 1, space = c("rgb", "Lab"),interpolate = c("linear", "spline"))
  

plot(rank(log(lakes[,1])),tti,col='white')
for(i in 1:ncol(d))
{
   ht_col<-cus_col(1000)[ floor( (mm[1:(nrow(mm)-1),i]) / ( max(mm[1:(nrow(mm)-1),i]) ) * (1000-1) + 1 ) ]
   points(rep(rank( log(lakes[,1]) )[i] ,max(d)-min(d)) ,min(d):(max(d)-1),col=ht_col,pch=15,cex=0.2)
}


plot(mm[1:(nrow(mm)-1),i],type='l',ylim=c(0,1000))
for(i in 1:ncol(mm))
lines(mm[1:(nrow(mm)-1),i],col=i)


