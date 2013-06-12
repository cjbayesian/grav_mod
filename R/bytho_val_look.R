
p_pars<-read.table('output/pred_pars.tab')

par(mfrow=c(2,2))
for(i in 1:4)
    hist(p_pars[,i],breaks=30)
x11()
plot(p_pars)


library(rvmapp)

p_hat<-as.matrix(read.table('output/val_sim_props.tab'))
o<-read.table('output/val_lakes.dat')
o<-o[,1]
v<-vmapp(pred=p_hat,d=o)

## if delta is outside of what is possible, replace by cutoffs
 v$delta[v$delta > v$pred]<-v$pred[v$delta > v$pred]
 v$delta[v$delta < -(1-v$pred)]<- -(1-v$pred[v$delta < -(1-v$pred)])


x11()
plot(p_hat[1,],v$delta[1,],ylim=c(-0.5,0.5),xlim=c(0,1))
for(i in 2:195)
    points(p_hat[i,],v$delta[i,])
abline(h=0)

## p-hat vs p (estimated)
plot(t(p_hat[1,]),t(v$delta[1,])+p_hat[1,],ylim=c(0,1),xlim=c(0,1))
for(i in 2:195)
    points(t(p_hat[i,]),t(v$delta[i,])+p_hat[i,])
abline(0,1)
##

mean_p_hat<-as.vector(apply(p_hat,2,mean))
mean(mean_p_hat)
mean(o)


## On a null model (p = p_hat)
p_hat<-runif(1000)
o<-bs(p_hat)
p_hat<-t(array(rep(p_hat,1000),dim=c(1000,1000)))

v<-vmapp(pred=p_hat,d=o)


par(mfrow=c(1,2))
mean_delta<-apply(v$delta,2,mean)
plot(p_hat[1,],mean_delta+p_hat[1,],ylim=c(0,1),xlim=c(0,1))
abline(0,1,lty=2)

## if delta is outside of what is possible, replace by cutoffs
 v$delta[v$delta > v$pred]<-v$pred[v$delta > v$pred]
 v$delta[v$delta < -(1-v$pred)]<- -(1-v$pred[v$delta < -(1-v$pred)])

mean_delta<-apply(v$delta,2,mean)
mean_p_hat<-apply(p_hat,2,mean)
plot(mean_p_hat,mean_delta+mean_p_hat,ylim=c(0,1),xlim=c(0,1))
abline(0,1,lty=2)

for(i in 1:1000)
    points(p_hat[i,],v$delta[i,]+p_hat[i,],col=rgb(0,0,1,0.01),pch=20)
#########################################



pdf('plots/bytho_val_look.pdf')
m<-nrow(v$discrepencies)
par(mfrow=c(3,1))
for(i in 1:m)
{
    plot(v$pred[i,],v$discrepencies[i,],main='f2',xlim=c(0,1))
    curve(F2(x,pars=v$f2_pars[i,]),add=T)

    plot(v$pred[i,],v$direction_discrepencies[i,],main='f1',xlim=c(0,1))
    curve(F1(x,pars=v$f1_pars[i,]),add=T)
    abline(v=min(v$pred[i,]),lty=2)
    abline(v=max(v$pred[i,]),lty=2)

    curve(2*(F1(x,pars=v$f1_pars[i,])-0.5 ) * F2(x,pars=v$f2_pars[i,]),
        xlim=c(0,1),
        ylim=c(-0.5,0.5))
    abline(h=0,lty=2)

    abline(v=min(v$pred[i,]),lty=2)
    abline(v=max(v$pred[i,]),lty=2)
}
dev.off()


## AUC ##
if(F)
{
    source('~/SchoolBackUp/AUC/AUC.R')
    #source('~/AUC/AUC.R')
    aucs<-numeric(nrow(p_hat))
    for(i in 1:nrow(p_hat))
        aucs[i]<-AUC(o,p_hat[i,])$auc

    hist(aucs,breaks=20,xlim=c(0.5,1))
} 

