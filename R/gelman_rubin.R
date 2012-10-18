####Convergence Diagnosis #####

library(coda)

d1<-read.table('output/lib.mcmcA')
d2<-read.table('output/lib.mcmcB')

#plot(d1[,ncol(d1)],d2[,ncol(d2)])


burn_in<-10000

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

names(d1)<-vars
names(d2)<-vars

d1<-mcmc(d1[burn_in:(nrow(d1)-1),2:(ncol(d1)-1)],thin=1)
d2<-mcmc(d2[burn_in:(nrow(d2)-1),2:(ncol(d2)-1)],thin=1)

chains<-as.mcmc.list(list(d1,d2))



gelman.diag(chains, 
   confidence = 0.95,
   transform=FALSE,
   autoburnin=TRUE,
   multivariate=TRUE)



