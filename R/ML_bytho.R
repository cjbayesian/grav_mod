################################################################
## Machine learning approaches to the Byho data
## Corey Chivers 
## <corey.chivers@mail.mcgill.ca>
################################################################

ultrabook<-TRUE

### Load AUC ###
sbu<-''
if(ultrabook)
   sbu<-'SchoolBackUp/'
AUC_source<-paste('~/',sbu,'AUC/AUC.R',sep='')
source(AUC_source)
################

lakes<-read.csv('../2010_bytho_data/lakes_processed.csv',sep='\t',header=FALSE)
env_vars<-lakes[,c(1:3,8:20)]


## Which lakes have been sampled (pre 2010)
sampled_index<-lakes[,5] != 0 | lakes[,6] !=0
sampled_label<-lakes[sampled_index,4]
inputs<-env_vars[sampled_index,]

d<-data.frame(sampled_label,inputs)




#####################################################################
library(randomForest)

rf<-randomForest(as.factor(sampled_label) ~ ., data = d,
   ntree = 1000, importance = TRUE)

rf

pred<-predict(rf,type='prob')
AUC(d=d[,1],pred=pred[,2],plot=T)


val_index<-lakes[,7] != 0
val_actual<-lakes[val_index,7] - 1
val_input<-env_vars[val_index,]

pred<-predict(rf,newdata=val_input,type='prob')
AUC(d=val_actual,pred=pred[,2],plot=TRUE,main='2010')

## Overall predicted rate ##
print(sum(pred[,2]))

cbvm(d=val_actual,pred=pred[,2])

#####################################################################
