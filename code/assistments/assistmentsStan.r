library(rstan)
options(mc.cores = parallel::detectCores())

load('assistmentsData/assistmentsAnalysisData.RData')

sdat <-
    with(newDat,
         list(
             nctl=sum(1-Z),
             ntrt=sum(Z),
             p=ncol(X),
             xt=X[Z==1,],
             xc=X[Z==0,],
             Ytrt=Y[Z==1],
             Yctl=Y[Z==0],
             St=S[Z==1],
             lowBound=0.001,
             upBound=0.999
         )
         )


mod <- stan('code/assistments.stan',data=sdat,chains=8,thin=2)

save(mod,file='assistmentsStan.RData')
