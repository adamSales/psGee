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
     cl <- makeCluster(10)
     clusterEvalQ(cl,library(tidyverse))
     clusterEvalQ(cl,library(rstan))
     clusterEvalQ(cl,library(geex))

     clusterExport(cl,c('simStan2step','makeDat','mle','twoStep','oneStep','stanMod'))
## }

fullsim(1000,cl=cl,ns=500,mu01=0,mu10=0,gumbs=FALSE,b1s=0.5,ext='smaller')

stopCluster(cl)

system('drive push --no-clobber -quiet')