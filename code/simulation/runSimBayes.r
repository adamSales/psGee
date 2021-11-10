library(dplyr)
library(rstan)
print('here')
library(geex)
library(parallel)
#library(pbapply)

#load('stanMod.RData')

rstan_options(auto_write = TRUE)

source('code/simFuncs.r')



fullsim(500,ns=c(500,1000),mu01=c(0,.3),mu10=0,b1s=c(0,0.2,0.5))

system('drive push --no-clobber -quiet')
