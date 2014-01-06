n_sampled <- floor(0.2*n_lakes) ## sample 20% of lakes
look_lakes <- sample(1:n_lakes,n_sampled)

found_year <- rep(0,n_lakes)
for(l in look_lakes)
{
    for(time in (from_year+1):to_year)
    {
        if(inv_year[l]<=time & inv_year[l] != 0)
        {
            found <- sample(c(TRUE,FALSE),1,prob=c(p2,1-p2))
            if(found)
            {
                found_year[l] <- time
                break
            }
        }
    }
}

plot(inv_year[found_year!=0],found_year[found_year!=0])
abline(0,1)



