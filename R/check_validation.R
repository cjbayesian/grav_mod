
lakes<-read.csv('../2010_bytho_data/lakes_processed.csv',sep='\t',header=FALSE)
pr<-read.csv('output/pred_p.tab',sep='\t')

pr_sim<-read.csv('output/val_sim_props.tab',sep='\t',head=FALSE)

pr<-as.matrix(pr)
pr<-pr[,1:(ncol(pr)-1)]

pr_sim<-as.matrix(pr_sim)
pr_sim<-pr_sim[,1:(ncol(pr_sim)-1)]

mean_pr<-apply(pr,2,mean)
ord<-order(mean_pr)

plot(pr[1,ord],pch=20,ylim=c(0,1))
for(i in 1:nrow(pr))
   points(pr[i,ord],pch=20,col=i)


for(i in 1:nrow(pr_sim))
   points(pr_sim[i,ord],pch=15,col=i+2)



val_index<-as.matrix(read.csv('output/pred_p.tab',sep='\t',header=FALSE))
val_index<-as.integer(val_index[1,1:(ncol(val_index)-1)])
points(lakes[val_index[ord],7]-1,pch=3)


## Overal rate of predicted invasions
#total_predicted<-apply(pr,1,sum)
#hist(total_predicted)
#abline(v=sum(lakes[val_index,4])) ## Add actual

########### AUC #################
source('~/AUC/AUC.R')
#source('~/SchoolBackUp/AUC/AUC.R')
x11()
##DIRECT CALCULATION##
d= lakes[val_index,7]==2
AUC(d=d,pred=pr[1,],plot=TRUE)
cal_auc<-numeric(nrow(pr))
for(i in 1:nrow(pr))
   cal_auc[i]<-AUC(d=d,pred=pr[i,],plot=TRUE,add=TRUE)$auc


x11()
##SIMS##
d= lakes[val_index,7]==2
AUC(d=d,pred=pr_sim[1,],plot=TRUE)
sim_auc<-numeric(nrow(pr_sim))
for(i in 1:nrow(pr_sim))
   sim_auc[i]<-AUC(d=d,pred=pr_sim[i,],plot=TRUE,add=TRUE)$auc


cal_rates<-apply(pr,1,sum)
sim_rates<-apply(pr_sim,1,sum)
actual_rate<-sum(d)

hist(cal_rates,breaks=30,xlim=c(0,25))
hist(sim_rates,breaks=30,add=TRUE,col='lightgrey')
abline(v=actual_rate,lwd=3,lty=2,col='blue')


########## New Metric ##############
source('~/AUC/bs.R')
source('~/AUC/ripley_style2.R')

rip_style2(d,pr,plot_met=TRUE)
x11()
rip_style2(d,pr_sim,plot_met=TRUE)

####################################


