options(mc.cores=6)
bayes <- function(data,...){
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

    ## get good initial values
    mm <- with(sdat,lmer(Ytrt~St+xt+(1|sclt)))
    mm2 <- with(sdat,glm(St~xt,family=binomial))

    init1 <- list(a2=fixef(mm)[-c(1,2)],
                  a01=fixef(mm)[1],
                  a11=fixef(mm)[2],
                  sigt=summary(mm)$sigma,
                  sigSclY=sqrt(VarCorr(mm)[[1]][1,1]),
                  sigc=summary(mm)$sigma,
                  b0=coef(mm2)[1],
                  b1=coef(mm2)[-1]
                  )

    stan('code/ct.stan',data=sdat,chains=6,init=list(init1,init1,init1,init1,init1,init1),...)
}
