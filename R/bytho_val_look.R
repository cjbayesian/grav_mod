
p_pars<-read.table('output/pred_pars.tab')

par(mfrow=c(2,2))
for(i in 1:4)
    hist(p_pars[,i],breaks=30)
x11()
plot(p_pars)


library(rvmapp)

p_hat<-read.table('output/val_sim_props.tab')
o<-read.table('output/val_lakes.dat')
o<-o[,1]
v<-vmapp(pred=p_hat,d=o)

plot(t(p_hat[1,]),t(v$delta[1,]),ylim=c(-0.5,0.5),xlim=c(0,1))
for(i in 2:195)
    points(t(p_hat[i,]),t(v$delta[i,]))
abline(h=0)

mean_p_hat<-as.vector(apply(p_hat,2,mean))
mean(mean_p_hat)
mean(o)


## AUC ##
source('~/SchoolBackUp/AUC/AUC.R')
#source('~/AUC/AUC.R')
aucs<-numeric(nrow(p_hat))
for(i in 1:nrow(p_hat))
    aucs[i]<-AUC(o,p_hat[i,])$auc

hist(aucs,breaks=20,xlim=c(0.5,1))
 

