
p_pars<-read.table('output/pred_pars.tab')

par(mfrow=c(2,2))
for(i in 1:4)
    hist(p_pars[,i],breaks=30)
x11()
plot(p_pars)


library(rvmapp)

p_hat<-read.table('output/val_sim_propsA.tab')
o<-read.table('output/val_lakes.dat')
o<-o[,1]
v<-vmapp(pred=p_hat,d=o)

x11()
plot(t(p_hat[1,]),t(v$delta[1,]),ylim=c(-0.5,0.5),xlim=c(0,1))
for(i in 2:195)
    points(t(p_hat[i,]),t(v$delta[i,]))
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
for(i in 1:10)
{
p_hat<-runif(100)
o<-bs(p_hat)
v<-vmapp(pred=p_hat,d=o)

mean_delta<-apply(v$delta,2,mean)
plot(p_hat,mean_delta+p_hat,ylim=c(0,1),xlim=c(0,1),main=i)
abline(0,1,lty=2)

}


pdf('plots/bytho_val_look.pdf')
m<-nrow(v$discrepencies)
par(mfrow=c(3,1))
for(i in 1:m)
{
    plot(t(v$pred[i,]),v$discrepencies[i,],main='f2',xlim=c(0,1))
    curve(F2(x,pars=v$f2_pars[i,]),add=T)

    plot(t(v$pred[i,]),v$direction_discrepencies[i,],main='f1',xlim=c(0,1))
    curve(F1(x,pars=v$f1_pars[i,]),add=T)
    abline(v=min(t(v$pred[i,])),lty=2)
    abline(v=max(t(v$pred[i,])),lty=2)

    curve(2*(F1(x,pars=v$f1_pars[i,])-0.5 ) * F2(x,pars=v$f2_pars[i,]),
        xlim=c(0,1),
        ylim=c(-0.5,0.5))
    abline(h=0,lty=2)

    abline(v=min(t(v$pred[i,])),lty=2)
    abline(v=max(t(v$pred[i,])),lty=2)
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

