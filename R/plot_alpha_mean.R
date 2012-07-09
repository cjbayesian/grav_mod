
#gen_alpha<-as.matrix(read.table("sims/alpha.dat"))

gen_pars<-as.matrix(read.table("sims/ch_params.csv"))

post_pars<-as.matrix(read.table("output/means.out"))
 
## Since means.out is in alpha and we 
## want beta_0, need to transform
## post_alpha<--log(1-(1/(1+exp(-post_pars))))

n_pars<-length(post_pars[1,])
pdf('plots/pred_obs_ch.pdf')
    for(i in 1:n_pars)
	    plot(gen_pars[,i],post_pars[,i],pch=20);abline(0,1,lty=2)
dev.off()


## Q-Q plot ##

quantiles<-as.matrix(read.table("output/quantiles.out"))
pdf('plots/qq.pdf')
    for(n in 1:length(quantiles[1,]))
    {
        exp_qs<-seq(0,1,0.01)
        actual_qs<-exp_qs
        i<-1
        for(p in exp_qs)
        {
            actual_qs[i]<-length(which(quantiles[,n]<=p))/length(quantiles[,n])
            i<-i+1
        }

        plot(exp_qs,actual_qs,type='b',xlim=c(0,1),ylim=c(0,1),pch=20,main=n)
        abline(0,1)
    }
dev.off()


