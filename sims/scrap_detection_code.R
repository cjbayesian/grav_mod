POE <- function(a_par,c_par,x)
{
    1-exp(-a_par*x^c_par)
}

N <- 100
TIME <- 30

pp <- runif(N,0,100)
poe <- POE(0.0001,1.5,pp)

inv <- array(dim=c(TIME,N))

for(time in 1:TIME)
    inv[time,] <- runif(N) < poe

year_inv <- apply(inv,2,function(x){ifelse(sum(x)==0,0,min(which(x))) })

p_detect <- 0.1
year_discovered <- 1:N
for(i in 1:N)
{
    if(year_inv[i] == 0){year_discovered[i]<-0} else{    
        detect <- runif(TIME-year_inv[i]) < p_detect
        year_discovered[i] <- ifelse(sum(detect)==0,0,min(which(detect))) + year_inv[i]
    }
}

plot(year_inv[year_discovered!=0],year_discovered[year_discovered!=0])


lhood <- function(pars)
{
    a_par = pars[1]
    c_par = pars[2]
    p_detect = pars[3]

    

}
