
d<-read.table('output/lib.mcmc')
burn_in<-5000
thinned<- seq(burn_in,nrow(d)-1,10)
d <- d[thinned,]
n_iter<-nrow(d)-1
length( unique(d[,2]))/length(d[,2] )


vars<-c("i",
    "d",
    "e",
    "c",
    "B_o",
    "NAUT",
    "KKUT",
    "MGUT",
    "CAUT",
    "PPUT1",
    "SIO3UR",
    "DOC",
    "COLTR",
    "ALKTI",
    "ALKT",
    "PH",
    "COND25",
    "SECCHI.DEPTH")

names(d) <- vars

par(mfrow=c(5,4))
for(i in 2:19)
{    plot(d[,i],type='l',main=names(d)[i])
    abline(h=0,lty=2)
}

#for(i in 2:19)
#    plot(d[sample(1:nrow(d)),i],type='l',main=names(d)[i])

par(mfrow=c(5,4))
for(i in 2:19)
    hist(d[,i],probability=TRUE,main=names(d)[i])

x11()
plot(d$c,d$B_o)# ,ylim=c(-14,-6),xlim=c(0.3,1.5))
 

plot(d)

### Traces ###
#cur_time<-format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
pdf(paste('plots/traces',cur_time,'.pdf',sep=''))
par(mfrow=c(2,1))
for(i in 2:ncol(d))
{
   plot(d[burn_in:n_iter,i],type='l',main=vars[i])
   hist(d[burn_in:n_iter,i],main='')
}

dev.off()

  
acceptance_rate<-length(unique(d[,2]))/length(d[,2] )
acceptance_rate

plot(d[seq(burn_in,nrow(d),100),]) ## Thinned correlogram


### Metaanalysis style (post mean and 95% BCI) ####
e_x<-apply(d[burn_in:n_iter,],2,mean)
q_x<-apply(d[burn_in:n_iter,],2,quantile,c(0.025,0.975))

par(mar=c(5, 7.7, 2, 2))
plot(e_x[6:18],6:18,
   xlim=c(min(q_x[1,6:18]),max(q_x[2,6:18]) ),
   ylim=c(5,19),
   ylab='',
   xlab='Posterior',
   pch=20,
   yaxt="n")

y<-6:18
axis(2, at=y,labels=vars[y], las=2)


for(i in 6:18)
{
   arrows(e_x[i],i,q_x[1,i],i,angle=90,length=0.07)
   arrows(e_x[i],i,q_x[2,i],i,angle=90,length=0.07)
}

abline(v=0,lty=2)


### correlogram ###
sub_index<-sample(2000:n_iter,1000)

d_sub<-d[sub_index,2:18]
plot(d_sub)


