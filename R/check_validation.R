
lakes<-read.csv('../2010_bytho_data/lakes_processed.csv',sep='\t',header=FALSE)
pr<-read.csv('output/pred_p.tab',sep='\t')

pr<-as.matrix(pr)
pr<-pr[,1:(ncol(pr)-1)]


mean_pr<-apply(pr,2,mean)
ord<-order(mean_pr)

plot(pr[1,ord],pch=20,ylim=c(0,1))
for(i in 1:nrow(pr))
   points(pr[i,ord],pch=20,col=i)

val_index<-as.matrix(read.csv('output/pred_p.tab',sep='\t',header=FALSE))
val_index<-as.integer(val_index[1,1:(ncol(val_index)-1)])
points(lakes[val_index[ord],4],pch=3)

x11()
## -- ##
## Inv prob as a function of lake size.
plot(lakes[val_index[ord],1],pr[1,ord],col=lakes[val_index[ord],4]+1,pch=20,ylim=c(0,1))
for(i in 1:nrow(pr))
{
   points(lakes[val_index[ord],1],pr[i,ord],col=lakes[val_index[ord],4]+1,pch=20)
   lines(lakes[val_index[ord],1],pr[i,ord])
}

## Overal rate of predicted invasions
#total_predicted<-apply(pr,1,sum)
#hist(total_predicted)
#abline(v=sum(lakes[val_index,4])) ## Add actual

########### AUC #################
#source('~/AUC/AUC.R')
source('~/SchoolBackUp/AUC/AUC.R')
x11()
AUC(d=lakes[val_index,4],pred=pr[1,],plot=TRUE)

