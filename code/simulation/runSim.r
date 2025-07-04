##########################################################################################
########### run the simulation
##########################################################################################


library(dplyr)
library(rstan)
library(parallel)


stanPSmod=stan_model("ps.stan",model_name="stanPSmod")

#rstan_options(auto_write = TRUE)

### number of repetitions
nrep <- 500

### For running with multiple cores, change the following to the number of cores you want to use:
ncore <-  1#


if(.Platform$OS.type=='windows' & ncore>1){
  print("making cluster")
  cl <- makeCluster(ncore)
  clusterEvalQ(cl,library(dplyr))
  clusterEvalQ(cl,library(rstan))
  clusterEvalQ(cl,rstan_options(auto_write=TRUE))
  clusterEvalQ(cl, stan_model('ps.stan'))
  clusterEvalQ(cl,source('regression.r'))
  clusterEvalQ(cl,source('simFuncs.r'))
} else cl <- NULL


source('simFuncs.r')

source('regression.r')


### this produces the results in Section 4.2.2
### the simulated datasets are saved in ../simData
### the results are saved in ../simResults
fullsim(nrep,ns=c(500,1000),mu01=c(0,.3),mu10=0,b1s=c(0,0.2,0.5),
	ext='',cl=cl,ncores=ncore,start=1)


### this adds PSW results:
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
### look at wider ranges of n and alpha
### this produces the results in Section 4.2.1

fullsim(nrep,ns=c(100,200,300,400,600,700,800,900),errDist='norm',b1s=.5,intS=FALSE,intZ=FALSE,mu01=.3,mu10=0,ext='ns',cl=cl,
        ncores=ncore)

fullsim(nrep,ns=c(200,300,400,600,700,800),errDist='norm',b1s=.5,intS=FALSE,intZ=FALSE,mu01=0,mu10=0,ext='ns_mu01is0',cl=cl,
            ncores=ncore)

fullsim(nrep,ns=500,errDist='norm',b1s=c(.1,.3,.4,seq(.6,1,.1)),intS=FALSE,intZ=FALSE,mu01=0,mu10=0,ext='b1ss_mu01is0',cl=cl,
            ncores=ncore)


if(!is.null(cl)) stopCluster(cl)


##########################################################################################
########### read in the results, pre-process
##########################################################################################

source("code/simulation/readSim.r")
### this produces files:
## simResults/fullResults.RData
## simResults/resultsNs_mu01is0.RData
## simResults/resultsB1s_mu01is0



##########################################################################################
########### make tables and figures
##########################################################################################
source("code/simulation/displaySimPublish.r")
