library(dplyr)
library(rstan)
print('here')
library(geex)
library(parallel)
#library(pbapply)

#load('stanMod.RData')

source('code/simFuncs.r')

if(file.exists('stanMod.RData')) load('stanMod.RData') else{ 
  stanMod <-  stan_model('code/ps.stan')
  save(stanMod,file='stanMod.RData')
  }
print('here')


#rstan_options("auto_write" = TRUE)
#pboptions(type='none')

## if(!exists('cl')){

ncore <-  detectCores()/2
cl <- makeCluster(ncore)
     clusterEvalQ(cl,library(tidyverse))
     clusterEvalQ(cl,library(rstan))
     clusterEvalQ(cl,library(geex))

     clusterExport(cl,c('simStan2step','makeDat','mle','twoStep','oneStep','stanMod'))
## }

### see what simulations have run already
fl <- list.files('./simResults',pattern="sim[0-9]+\\.RData")
nums <- readr::parse_number(fl)

print(paste("starting at case",max(nums,na.rm=TRUE)+1))

fullsim(1000,cl=cl,start=max(nums,na.rm=TRUE)+1)

fullsim(1000,mu11=0,mu01=0.3,mu10=0.3,ext='reverse',cl=cl)

stopCluster(cl)

system('drive push --no-clobber -quiet')