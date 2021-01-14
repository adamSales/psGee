

rstan_options("auto_write" = TRUE)


mle <- function(data){
    sdat <- with(data, list(
                           nctl=sum(1-Z),
                           ntrt=sum(Z),
                           x1t=x1[Z==1],
                           x1c=x1[Z==0],
                           x2t=x2[Z==1],
                           x2c=x2[Z==0],
                           Ytrt=Y[Z==1],
                           Yctl=Y[Z==0],
                           St=S[Z==1]
                       )
                 )
                                        #fit1 <-
    fit <- optimizing(stanMod,data=sdat)
#    fit2 <- optimizing(stanMod,data=sdat)

 #   fit <- if(fit1$value>fit2$value) fit1 else fit2

    c(eff0=fit$par['mu10']-fit$par['mu00'],
      eff1=fit$par['mu11']-fit$par['mu01'])
}

stanMod <- stan_model('ps.stan',auto_write=TRUE)

#fit <- optimizing(stanMod,data=sdat)
