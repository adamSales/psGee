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
     cl <- makeCluster(50)
     clusterEvalQ(cl,library(tidyverse))
     clusterEvalQ(cl,library(rstan))
     clusterEvalQ(cl,library(geex))

     clusterExport(cl,c('simStan2step','makeDat','mle','twoStep','oneStep','stanMod'))
## }

fullsim(1000,cl=cl)

fullsim(1000,mu11=0,mu01=0.3,mu10=0.3,ext='reverse',cl=cl)

stopCluster(cl)

system('drive push --no-clobber -quiet')