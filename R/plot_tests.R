
sim=TRUE

## Plot traces ##
plot_traces<-function(sim=TRUE)
{
    if(sim)
    {
        sim_inv<-as.matrix(read.csv('sims/inv_year.csv',header=FALSE))
        sim_inv[sim_inv==0]<-2011
        last_abs<-as.matrix(read.csv('sims/last_obs.csv',header=FALSE)) ##last abs.
        ch_params<-as.matrix(read.csv('sims/ch_params.csv',header=FALSE)) ##generating ENV params.
    }
    ### STATESPACE ###
    if(T)
    {

        tl<-as.matrix(read.table('t_mcmc.dat',header=F))
        n_lakes<-length(tl[1,])
        pdf('trace.pdf')
        par(mfrow=c(5,5),mar=c(0,0,0,0))
        n_lakes<-length(tl[1,])
        predict_t<-numeric(n_lakes)
        ucl<-numeric(n_lakes)
        lcl<-numeric(n_lakes)
         for(i in 1:n_lakes )
         {
            plot(tl[,i],ylim=c(1989,2025),type='l')
            abline(h=2010,lty=2)
            

            if(sim)
            {
                abline(h=sim_inv[i],col='blue',lwd=2)
                abline(h=last_abs[i],col='green',lwd=2)
            }

            abline(h=mean(tl[,i]),col='red',lwd=2,lty=3)    
            predict_t[i]<-mean(tl[,i])
            ucl[i]<-quantile(tl[,i],probs=0.975)
            lcl[i]<-quantile(tl[,i],probs=0.025)
            text(0.5*length(tl[,i]),2015,i)
         }
        dev.off()

        if(sim)
        {
            jit<-runif(n_lakes,-1,1) #jitter
            plot(sim_inv+jit,predict_t,xlab='actual',ylab='predicted') 
            for(i in 1:n_lakes )
            {
                arrows(sim_inv[i]+jit[i],ucl[i],sim_inv[i]+jit[i],lcl[i],length=0.05,angle=90,code=3)
            }
            abline(0,1)
        }

    }

}

if(F)
{
    pdf('t_hist.pdf')
    par(mfrow=c(5,5),mar=c(0,0,0,0))
    for(i in 1:n_lakes )
    {
        hist(tl[,i],breaks=seq(1988.5,2011.5,1),probability=TRUE)
    }
    dev.off()
}

### IN GEOSPACE ###
geo_ani<-function(sim=FALSE)
{
    if(sim)
        {lakes<-as.matrix(read.csv('sims/simmed_lakes.csv',sep='\t',header=FALSE))}else{
    lakes<-as.matrix(read.csv('../2010_bytho_data/lakes_processed.csv',sep="\t",header=FALSE))}
    tl<-as.matrix(read.table('output/t_mcmc.dat',header=F))
    tl<-tl[1:(nrow(tl)-1),]
    n_lakes<-length(tl[1,])
    for(t in min(tl):(max(tl)-1))
    {
        png(paste('plots/images/posterior_',t,'.png',sep=''))
        prop_inv<-sapply(1:n_lakes,function(x){sum(tl[,x]<=t)})
        prop_inv<-prop_inv/length(tl[,1])
        risk_order<-order(prop_inv)
        plot(lakes[risk_order,2],lakes[risk_order,3],pch=1,main=t,cex=1,xlab='UTM X',ylab='UTM Y')
        points(lakes[risk_order,2],lakes[risk_order,3],col=grey(1-prop_inv[risk_order]),pch=20,main=t,cex=1)
        dev.off()
    }
    system('convert -delay 60 plots/images/posterior_* plots/images/post.gif')
}


cummulative_total<-function()
{
   tl<-as.matrix(read.table('output/t_mcmc.dat',header=F))
   tl<-tl[1:(nrow(tl)-1),]

   #Discovery Rate
   lakes<-as.matrix(read.csv('../2010_bytho_data/lakes_processed.csv',sep="\t",header=FALSE))
   disc_year<-lakes[,5]

   n_inv<-numeric(length(1989:2010))
   n_disc<-numeric(length(1989:2010))
   for(time in 1989:2010)
   {
      n_inv[time-1988]<-sum(tl<=time)/nrow(tl)
      n_disc[time-1988]<-sum(disc_year<=time & disc_year != 0)
   }
   plot(1989:2010,n_inv,type='l',lty=2,ylim=c(0,100),xlab='Year',ylab='# invaded')    
   lines(1989:2010,n_disc)

   if(FALSE) ## Plot the trace expected time to invasion across all the sites.
   {
      x<-apply(tl,1,sum)
      plot(x/ncol(tl))
   }
}

top_n_risk<-function(top_n=500)
{
   tl<-as.matrix(read.table('output/t_mcmc.dat',header=F))
   tl<-tl[1:(nrow(tl)-1),]

   #Discovery Rate
   lakes<-as.matrix(read.csv('../2010_bytho_data/lakes_processed.csv',sep="\t",header=FALSE))
   disc_year<-lakes[,5]
   last_abs<-lakes[,6]
   avg_inv_year<-apply(tl,2,mean)
   lower_inv_year<-apply(tl,2,quantile, probs = 0.05)
   upper_inv_year<-apply(tl,2,quantile, probs = 0.95)
   ever_sampled<- (last_abs == 0 & disc_year == 0)
   ever_sampled<-ever_sampled+1

   risk_order<-order(avg_inv_year)

   pdf('plots/top_n_risk.pdf',width=6,height=12)
   
   plot(avg_inv_year[risk_order[1:top_n]],1:top_n,
      xlim=c(1988,2012),
      xlab='Year',
      pch=20,
      col=ever_sampled[risk_order[1:top_n]],
      ylab=paste('Top',top_n,'lakes at risk'))
   abline(v=1989:2011,lty=4,col='grey')
   for(i in 1:top_n)
   {
      arrows(avg_inv_year[risk_order[i]],i,lower_inv_year[risk_order[i]],angle=90,length=0.01,lty=3,col=ever_sampled[risk_order[i]])
      arrows(avg_inv_year[risk_order[i]],i,upper_inv_year[risk_order[i]],angle=90,length=0.01,lty=3,col=ever_sampled[risk_order[i]])
   }
   points(disc_year[risk_order[1:top_n]],1:top_n,pch=4,col='black')
   points(last_abs[risk_order[1:top_n]],1:top_n,pch=2,col='black')

   ##
   legend('topleft',bg='white',legend=c('Predicted timing of esablisment','Discovered','Last observed absence'),
      pch=c(20,4,2),col=c('black','black','black'))
   dev.off()
}


## Correlations ##
if(F)
{
    start_at<-300
    par(mar=c(0,0,0,0),mfrow=c(6,6))
    for(i in start_at:(start_at+5))
    {
        for(j in start_at:(start_at+5))
        {
            plot(tl[,i],tl[,j])
        }
    }
}
### Likelihood of MCMC trace ###
plot_lhood<-function(...)
{
    ll<-as.matrix(read.table('l_mcmc.dat'))
    plot(ll,type='l',...)
}
# Call to watch the lhood trace:    
watch_lhood<-function()
{
     while(T){plot_lhood(xlim=c(0,10000));Sys.sleep(2)}
}


## PARS ##
plot_pars<-function()
{

    sim_inv<-as.matrix(read.csv('sims/inv_year.csv',header=FALSE))
    sim_inv[sim_inv==0]<-2011
    last_abs<-as.matrix(read.csv('sims/last_obs.csv',header=FALSE)) ##last abs.
    ch_params<-as.matrix(read.csv('sims/ch_params.csv',header=FALSE)) ##generating ENV params.

    tl<-as.matrix(read.table('t_mcmc.dat',header=F))
    n_lakes<-length(tl[1,])

    par_mcmc<-as.matrix(read.table('par_mcmc.dat',header=F))
    pdf('par_traces.pdf')
    for(i in 1:15)
    {
        plot(par_mcmc[,i],type='l', main=i)  ## alpha (chem_0 intercept)
        if(i <= length(ch_params))
            abline(h=ch_params[i],lty=2) ## Generating values
    }
    dev.off()

    #plot(par_mcmc[,14],type='l',main='d_par') ## d_par
    #plot(par_mcmc[,15],type='l',main='e_par') ## e_par

    alphas_mcmc<-as.matrix(read.table('alpha_mcmc.dat',header=F))
    pdf('alphas_trace.pdf')
    for(i in 1:n_lakes)
        plot(alphas_mcmc[,i],type='l', main=i)  ## alpha (chem_0 intercept)
    dev.off()
}

## Watch MCMC trace as it is running.
watch_trace<-function(dim,delay=4,burn_in=100)
{
    n_pars<-length(dim)
    par(mfrow=c(n_pars,1))
	while(T)
	{
		try(par_mcmc<-as.matrix(read.table('par_mcmc.dat',header=F)))
        for(i in 1:n_pars)
	        try(plot(par_mcmc[burn_in:length(par_mcmc[,1]),dim[i]],type='l'))
        Sys.sleep(delay)
    }
}






if(F)
{

    ## Likelihoods profile testing ##
    l<-read.table('ll_test.dat',header=F)

    ## Does distance matter?
    ## What is the amount of information contained in each lake? Do bigger lakes contain more information (contrib to likelihood)
    ## or is it simply the propagule pressure contribution?
    #library(scatterplot3d)
    #scatterplot3d(l[, 1], l[, 2], l[, 3])

    gridtosurface<-function(d)
    {
        z<-outer(unique(d[,2]),unique(d[,1]),function(x,y){d[,3]})
        return(list(x=unique(d[,2]),y=unique(d[,1]),z=z))
    }

    z<-gridtosurface(l)
    contour(z)





    plot(l[,1],l[,2],pch=20)

    plot(l[,2],l[,3],col=l[,1],pch=20)

    for(lake in unique(l[,1]))
    {
        index<-which(l[,1]==lake)
        lines(l[index,2],l[index,3],col=lake)
    }


}


## 1D Lhood slice plot
simple_ll<-function()
{
    ll<-read.table('ll_test.dat')
    plot(ll[,1],ll[,2])
}

## For Visualizing Posterior parameters.
effect_size<-function()
{
     par_mcmc<-as.matrix(read.table('par_mcmc.dat',header=F))
     ndims<-length(par_mcmc[1,])
     min_par<-min(par_mcmc)
     max_par<-max(par_mcmc)
     plot(0,1,xlim=c(min_par,max_par),ylim=c(0,(ndims+2)),col='white' )
     for(i in 1:ndims)
     {
          mu<-mean(par_mcmc[,i])
          points(mu,i)
          ci<-quantile(par_mcmc[,i],c(0.05,0.95))
          arrows(mu,i,ci[1],i,angle=90,length=0.05)
          arrows(mu,i,ci[2],i,angle=90,length=0.05)
     }
     abline(v=0,lty=3)
}

