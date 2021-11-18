library(dplyr)
library(rstan)
print('here')
#library(geex)
library(parallel)
#library(pbapply)

#load('stanMod.RData')

#rstan_options(auto_write = TRUE)


ncore <-  8#detectCores()/2
#cl <- makeCluster(ncore)
#     clusterEvalQ(cl,library(tidyverse))
#     clusterEvalQ(cl,library(rstan))
#     clusterEvalQ(cl,rstan_options(auto_write=TRUE))
#     clusterEvalQ(cl, stan_model('code/ps.stan'))
#     clusterEvalQ(cl,source('code/regression.r'))
#     clusterEvalQ(cl,source('code/simulation/simFuncs.r'))



source('code/simulation/simFuncs.r')

source('code/regression.r')

### see what simulations have run already
#fl <- list.files('./simResults',pattern="sim[0-9]+\\.RData")
#nums <- readr::parse_number(fl)

#print(paste("starting at case",max(nums,na.rm=TRUE)+1))

fullsim(500,ns=c(500,1000),mu01=c(0,.3),mu10=0,b1s=c(0,0.2,0.5),
	ext='',cl=NULL,start=1)

system('drive push --no-clobber -quiet')
