### Test an ANN on the 2010 validation set ###
source("~/AUC/AUC.R")
library(ANN)

lakes<-read.csv('../2010_bytho_data/lakes_processed.csv',sep='\t',header=FALSE)

## Define which env variables to include ##
var_cols<-c(1,10,15,16,17)
n_vars<-length(var_cols)

## Normalize predictors [0-1] ##
env_vars_normalized<-lakes[,var_cols]
for(i in 1:n_vars)
   env_vars_normalized[,i]<-(lakes[,var_cols[i]]-min(lakes[,var_cols[i]]) )/(max(lakes[,var_cols[i]])-min(lakes[,var_cols[i]]))
## ---- ##


## Which lakes have been sampled (pre 2010)
sampled_index_all<-lakes[,5] != 0 | lakes[,6] !=0
tmp<-sample(which(sampled_index_all),sum(sampled_index_all))


## Split into training/test
sampled_index<-tmp[1:200]
test_index<-tmp[201:sum(sampled_index_all)]

sampled_label<-lakes[sampled_index,4]
inputs<-env_vars_normalized[sampled_index,]
## ---- ##


fit<-ANNGA ( inputs,
      as.numeric(sampled_label),
      design = c(n_vars, 3+n_vars, 1),
      minW=-10,
      maxW=10)

## Internal fit -- performance on training data ##
preds<-predict.ANN(fit, inputs)
AUC(d=sampled_label,pred=preds$predict,plot=TRUE,main='Training')


## Performance on holdout set ##
val_input<-env_vars_normalized[test_index,]
preds<-predict.ANN(fit, val_input)
AUC(d=lakes[test_index,4],pred=preds$predict,plot=TRUE,main='Holdout')



##############################################################
## 2010 Validation lakes ##

val_index<-lakes[,7] != 0
val_actual<-lakes[val_index,7] - 1
val_input<-env_vars_normalized[val_index,]

preds<-predict.ANN(fit, val_input)
AUC(d=val_actual,pred=preds$predict,plot=TRUE,main='2010')

## Overall predicted rate ##
print(sum(preds$predict))

