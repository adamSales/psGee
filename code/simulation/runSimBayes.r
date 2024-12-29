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


#system('drive push --no-clobber -quiet')


pswResults=sapply(seq(max(as.numeric(gsub('dat|\\.RData','',list.files('simData/'))),na.rm=TRUE)),
                  function(i){
                    print(i)
                    load(paste0('simData/dat',i,'.RData'))
                    facs=attributes(datasets[[1]])$facs
                    out=as.data.frame(t(vapply(datasets,psw,psw(datasets[[1]]))))
                    attr(out,'facs')=facs
                    out
                  },simplify=FALSE)
save(pswResults,file='simResults/pswResults.RData')


################################
### some further investigations
fullsim(500,ns=100,mu01=.3,mu10=0,errDist=c('norm','unif'),b1s=c(0,.2,.5),ext='n100')

fullsim(500,ns=c(200,300,400),errDist='norm',b1s=.5,intS=FALSE,intZ=FALSE,mu01=.3,mu10=0,ext='ns')

fullsim(500,ns=c(200,300,400,600,700,800),errDist='norm',b1s=.5,intS=FALSE,intZ=FALSE,mu01=.3,mu10=0,ext='ns',start=4)


fullsim(500,ns=500,errDist='norm',b1s=c(.3,.4,seq(.6,1,.1)),intS=FALSE,intZ=FALSE,mu01=.3,mu10=0,ext='b1s')
