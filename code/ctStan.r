options(mc.cores=6)
bayes <- function(data){
    require(rstan)

    form <- Y~state+grade+race+sex+frl+xirt+esl
    data <- droplevels(data)
    X <- model.matrix(form,data=data)[,-1]
    p <- ncol(X)

    X <- scale(X)

    scl <- as.numeric(as.factor(data$schoolid2))

    sdat <- with(data,
                 list(
                     nctl=sum(treatment==0),
                     ntrt=sum(treatment),
                     p=p,
                     nscl=max(scl),
                     xt=X[treatment==1,],
                     xc=X[treatment==0,],
                     sclt=scl[treatment==1],
                     sclc=scl[treatment==0],

                     Ytrt=Y[treatment==1],
                     Yctl=Y[treatment==0],

                     St=ifelse(everCP[treatment==1],1,0)
                 ))

    stan('code/ct.stan',data=sdat)
}
