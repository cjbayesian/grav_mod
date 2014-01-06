########### Sim data for testing Bayesian Gravity model ############
rm(list=ls())
n_sources<-10
n_lakes<-1000
pres_only<-TRUE
alpha<-runif(1,0,0.25)
error_sd<-0.75

#### sim sources ####

    sources<-data.frame(x=runif(n_sources,574653.0,693691.5),y=runif(n_sources,4972658,5056544)) 

    # sum(oi) = 40937 #actual
    sources$oi<-floor(runif(n_sources,0,(sources$x+sources$y)/10000))
    write.table(sources$oi,file='Oi.csv',row.names=FALSE,col.names=FALSE,sep='\t')

#### sim lakes ####

    lakes<-data.frame(x=runif(n_lakes,574653.0,693691.5),y=runif(n_lakes,4972658,5056544))

    # Fit size distribution to data
    lakes_data<-read.csv("../../2010_bytho_data/lakes_processed.csv",header=FALSE,sep='\t')

    lakes$area<-sample(lakes_data[,1],n_lakes)

    # Chem data
    # for now just random
    chem<-array(rnorm(n_lakes*13,0,2),dim=c(n_lakes,13))

    # Sample from actual Env data to mimic skewed distribution of some variables
    #for(ch in 1:13)
    #    chem[,ch]<-sample(lakes_data[,ch+6],n_lakes)

    ch_params<-runif(13+1,-0.4,0.4)
    #ch_params[1]<- runif(1,-4,-1) #Intercept -> alpha btw 0.01814993 and 0.3132617
    ch_params[1]<- -0.2

    ch_params[2:14]<-0 ## to test on intercept only (equiv to testing common alpha)
    #ch_params[2]<- 1

    ## Logit function (alpha=f(chem))
    logit_fn<-function(z)
    {
        return(-log(1-(1/(1+exp(-z) ) ) ) )
    }
    
    alphas<-numeric(n_lakes)
    z_lake<-numeric(n_lakes)
    epsilon<-rnorm(n_lakes,0,error_sd)
    for(l in 1:n_lakes)
    {
        z_lake[l]<-ch_params[1]+sum(ch_params[2:14]*chem[l,]) #+ epsilon[l]
        alphas[l]<-logit_fn(z_lake[l])
    }

    write.table(t(c(ch_params,2,0.8)),file='ch_params.csv',row.names=FALSE,col.names=FALSE,append=TRUE)
    #write.table(c(alphas,0.8,2),file='ch_params.csv',row.names=FALSE,col.names=FALSE)

    print(1-sd(epsilon)/(sd(z_lake)) )   ## R^squared
    write.table(1-sd(epsilon)/(sd(z_lake)), file='rsq.csv',row.name=FALSE,col.names=FALSE,append=TRUE) 
    
#### gen distmat ###

    euclidian<-function(from_x,from_y,to_x,to_y)
    {
            return( sqrt( (from_x-to_x)^2 + (from_y-to_y)^2 ) )
    }
    dist_mat<-array(dim=c(n_sources,n_lakes))
    for(b in 1:n_sources)
    {
            dist_mat[b,]<-euclidian( sources$x[b],sources$y[b],lakes$x,lakes$y )
    }
    write.table(dist_mat,file='distance_matrix.csv',row.names=FALSE,col.names=FALSE,sep='\t')


##### Gravity function #######
    grav_func<-function(e=0.3,d=2.5)
    {
        pred_p<-array(dim=c(n_sources,n_lakes))
        for(i in 1:n_sources)
        {
            for(l in 1:n_lakes)
            {
                pred_p[i,l]<-lakes$area[l]^e*dist_mat[i,l]^(-d)
            }
            pred_p[i,]<-pred_p[i,] / sum(pred_p[i,])
        }
        return(pred_p)
    }

    pred_p<-grav_func()


##### sim spread ######

    to_year<-2010
    from_year<-1989

    non_spread_pp<-array(dim=c(n_lakes,(to_year-from_year)+1))
    non_spread_pp[,1]<-rgamma(n_lakes,0.5,1) /15
    for(year in 2:(to_year-from_year+1))
        non_spread_pp[,year]<-non_spread_pp[,year-1] #+rgamma(n_lakes,0.5,25) /25
    
    write.table(non_spread_pp,file='pp.dat',col.names=FALSE,row.names=FALSE,sep='\t')

    seed_lake<-which.max(lakes$area)
    c_par<-2
    all_inv<-TRUE
    none_inv<-TRUE
    while(all_inv || none_inv)
    {
        invaded<-rep(FALSE,n_lakes)
        invaded[seed_lake]<-TRUE #Seed invaded lake

        index_invaded<-which(invaded)
        index_uninvaded<-which(!invaded)
        ind_inv_all<-NULL
        n_inv<-from_year:to_year
        n_inv[1]<-1
        inv_year<-rep(0,n_lakes)
        inv_year[seed_lake]<-from_year
        invaded_prior<-rep(FALSE,n_lakes)
        invaded_prior[seed_lake]<-TRUE
        for(year in (from_year+1):to_year)
        {
                    ## Calculate propagule pressure and invaded stochastically ##
                    X_i_t<-sapply(1:n_sources,function(i){sum(pred_p[i,index_invaded])})
                    invaded[index_uninvaded]<-sapply(index_uninvaded,function(x){
                        pp<-sum( sources$oi * (pred_p[,x])* X_i_t ) ##log joint prob of visiting an invaded lake then each uninvaded lake.
                        #pp<-non_spread_pp[x,year-(from_year-1)]

                        return( 1-exp(-((alphas[x]*pp)^c_par) ) > runif(1,0,1) )
                    })

                    index_invaded<-which(invaded)
                    if(length(index_invaded)==n_lakes)
                    {
                        all_inv<-TRUE
                        break
                    }
                    new_inv<-as.logical(invaded-invaded_prior)
                    inv_year[new_inv]<-year
                    invaded_prior<-invaded
                    index_uninvaded<-(1:n_lakes)[-index_invaded]

                    n_inv[year-from_year+1]<-length(index_invaded)
                    ind_inv_all<-c(ind_inv_all,index_invaded)
                    all_inv<-FALSE
        }
        print(inv_year)
        if(n_inv[to_year-from_year+1]>1){none_inv=FALSE}
    }


#### Sim samples ####
    diff_sample=TRUE ## if false, the same set of lakes is sampled at each time interval.
    #N_obs<-c(5,10,15,20,25)
    #N_obs<-c(25,20,15,10,5)
 #   N_obs<-c(75,75,75,75,75)
    sample_years<-c(1995,2005,2010)
    #sample_years<-1991:2010
    N_obs<-rep(100,length(sample_years))

    # Corresponds to columns 4,5,6 in data table
    invaded<-rep(0,n_lakes)
    first_obs_inv<-rep(0,n_lakes)
    last_obs_uninv<-rep(0,n_lakes)

    # Include seed lake as observed at t=0
    invaded[seed_lake]<-1
    first_obs_inv[seed_lake]<-from_year


    sampled_lakes<-sample(1:n_lakes,N_obs[1],prob=log(lakes$area))
    for(i in 1:length(sample_years)) ## 5 sampling times   
    {
        if(diff_sample)
            sampled_lakes<-sample(1:n_lakes,N_obs[i]) #,prob=log(lakes$area))

        for(l in sampled_lakes)
        {
            if(inv_year[l] <= sample_years[i] && inv_year[l] != 0 && first_obs_inv[l] == 0)
            {
                first_obs_inv[l]<-sample_years[i]
                invaded[l]<-1
            }
            if(inv_year[l] > sample_years[i] || inv_year[l] ==0)
            {
                last_obs_uninv[l]<-sample_years[i]
            }
        }
    }

### Identify lakes that were sampled in 2010 (validation set)
is_validation <- last_obs_uninv == 2010 | first_obs_inv == 2010
validation_index <- which(is_validation)

val_col <- rep(0,n_lakes)
val_col[is_validation] <- 1
val_col[is_validation & first_obs_inv] <- 2


#### Mimic the presence only data generation process with some
#### detection probability
if(pres_only)
{
    p2 <- 0.2 ## Detection Probability
    source('presence_only.R')
    first_obs_inv <- found_year
    first_obs_inv[seed_lake] <- from_year
}
####


#### Wrap up lakes info and write to file ####
    write_lakes<-data.frame(signif(lakes$area,4),lakes$x,lakes$y,invaded,first_obs_inv,last_obs_uninv,val_col)
    write_lakes<-cbind(write_lakes,signif(chem,3))

    write_lakes[write_lakes[,5]==2010,4] <- 0
    write_lakes[write_lakes[,6]==2010,6] <- 0
    write_lakes[write_lakes[,5]==2010,5] <- 0

    write.table(write_lakes,file='simmed_lakes.csv',col.names=FALSE,row.names=FALSE,sep='\t')

    # Validation data
    write.table(invaded[validation_index],file='val_lakes.dat',col.names=FALSE,row.names=FALSE,sep='\t')

    # for ploting diagnostics
    write.table(inv_year,file='inv_year.csv',row.names=FALSE,col.names=FALSE)
    write.table(last_obs_uninv,file='last_obs.csv',row.names=FALSE,col.names=FALSE)

    write.table(alphas[1],file='alpha.dat',append=TRUE,col.names=FALSE,row.names=FALSE)
