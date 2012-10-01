####Convergence Diagnosis #####

library(coda)

d1<-read.table('output/lib.mcmc')
d2<-read.table('output/lib.mcmcB')


burn_in<-10000
n_iter<-nrow(d)-1

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

