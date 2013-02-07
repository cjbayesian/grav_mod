

d<-read.csv("output/holdout2006_data_status.csv",header=FALSE)
pred<-read.csv("output/holdout_sim_props.csv",header=FALSE)

source("~/AUC/AUC.R")


tmp_pred<-as.vector(pred[1,])
tmp_d<-as.vector(d[1,])

##plot(tmp_pred,tmp_d)
print(AUC(d=tmp_d,pred=tmp_pred,plot=TRUE))

auc_vals<-numeric(length(pred[,1]))
pred_tot<-numeric(length(pred[,1]))
obs_tot<-numeric(length(pred[,1]))

for(i in 1:length(pred[,1]))
{
   tmp_pred<-as.vector(pred[i,])
   tmp_d<-as.vector(d[i,])

   ##plot(tmp_pred,tmp_d)
   auc<-AUC(d=tmp_d,pred=tmp_pred,plot=TRUE,add=TRUE,col='grey')
   print(auc$auc)
   auc_vals[i]<-auc$auc

   pred_tot[i]<-sum(tmp_pred)
   obs_tot[i]<-sum(tmp_d)
}


hist(auc_vals,breaks=20,xlim=c(0,1))

x11()

plot(pred_tot,obs_tot)
abline(0,1)



