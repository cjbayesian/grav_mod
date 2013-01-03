### Calculate and plot val metric on simulated data ###

source('~/SchoolBackUp/AUC/bs.R')
source('~/SchoolBackUp/AUC/ripley_style2.R')
#source('~/AUC/ripley_style2.R')

## Read in model predictions
pred<-as.matrix(read.table('../output/pred_p.tab',header=TRUE))

## Read in validation outcomes
#val_data<-read.table("../sims/val_lakes.dat")
val_data<-read.table("../output/val_lakes.dat")
val_data<-as.vector(val_data[,1])


### Alt-metric ###
rip_style2(val_data,pred,plot_met=TRUE)


### AUC ###
x11()
source('~/SchoolBackUp/AUC/AUC.R')
AUC(val_data,pred[1,],plot=TRUE)

for(i in 2:nrow(pred))
   AUC(val_data,pred[i,],plot=TRUE,add=TRUE)




