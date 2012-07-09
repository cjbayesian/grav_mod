

## Get true params
gen_pars<-as.matrix(read.table("sims/ch_params.csv"))

## Get Trace
burn_in=500
par_mcmc<-as.matrix(read.table('output/par_mcmc.dat',header=F))

index<-burn_in:length(par_mcmc[,1])

rep_i<-nrow(gen_pars)
n_par<-length(par_mcmc[1,])

quantile<-numeric(n_par)
for(q in 1:n_par)
    quantile[q]<-length(which(par_mcmc[index,q]<gen_pars[rep_i,q]))/length(par_mcmc[index,q])


write.table(t(quantile),"output/quantiles.out",append=TRUE,col.names=FALSE,row.names=FALSE)


# Posterior Means
mu<-numeric(n_par)
for(i in 1:n_par)
	mu[i]<-mean(par_mcmc[index,i])

#mu<- -log( 1-(1/(1+exp(-mu))) ) # Translate into alpha

cat("Est: ",mu,"\n")
cat("Gen: ",gen_pars[length(gen_pars[,1]),],"\n")
write.table(t(mu),"output/means.out",append=TRUE,col.names=FALSE,row.names=FALSE)


