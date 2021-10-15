library(rstan)


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


mod <- stan_model('code/assistments.stan')

fit <- optimizing(mod,data=sdat,hessian=TRUE)

save(mod,fit,file='assistmentsStan.RData')
