library(dplyr)
library(rstan)
print('here')
library(geex)
library(parallel)
#library(pbapply)

#load('stanMod.RData')

rstan_options(auto_write = TRUE)

source('code/simFuncs.r')


ncore <-  detectCores()/2
cl <- makeCluster(ncore)
     clusterEvalQ(cl,library(tidyverse))
     clusterEvalQ(cl,library(rstan))
     clusterEvalQ(cl,library(geex))

     clusterExport(cl,c('simOneBayes','makeDat','oneStep','bayes'))
## }

### see what simulations have run already
#fl <- list.files('./simResults',pattern="sim[0-9]+\\.RData")
#nums <- readr::parse_number(fl)

#print(paste("starting at case",max(nums,na.rm=TRUE)+1))

fullsim(500,cl=cl,ns=500,mu01=c(0,.3),mu10=0,gumbs=TRUE,b1s=c(0.2,0.5),ext='Bayes',Bayes=TRUE)


stopCluster(cl)

system('drive push --no-clobber -quiet')