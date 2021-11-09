
mle <- function(data,se=TRUE){
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
    fit <- rstan::optimizing(stanMod,data=sdat,hessian=se)
#    fit2 <- optimizing(stanMod,data=sdat)

 #   fit <- if(fit1$value>fit2$value) fit1 else fit2

    if(se){
        invHess <- -solve(fit$hessian)

        return(
            c(eff0=fit$par['mu10']-fit$par['mu00'],
              se0=sqrt(invHess['mu10','mu10']+invHess['mu00','mu00']-2*invHess['mu10','mu00']),
              eff1=fit$par['mu11']-fit$par['mu01'],
              se1=sqrt(invHess['mu11','mu11']+invHess['mu01','mu01']-2*invHess['mu11','mu01'])
              )
        )
    }

    c(eff0=fit$par['mu10']-fit$par['mu00'],
      eff1=fit$par['mu11']-fit$par['mu01']
      )
}

#stanMod <- rstan::stan_model('ps.stan')#,auto_write=TRUE)

bayes <- function(data,...){
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
    fit <- stan('code/ps.stan',data=sdat)

    summary(fit, par=c('eff0','eff1','effDiff'))$summary
}



twoStep <- function(data,justEst=TRUE){
    mod1 <- lm(Y~x1+x2+S,data=data,subset=Z==1)
    lmod <- glm(S~x1+x2,data=data,subset=Z==1,family=binomial)

    dat1 <- data.frame(
        r=data$Y[data$Z==0]-#coef(mod1)['x1']*data$x1[data$Z==0]-coef(mod1)['x2']*data$x2[data$Z==0],
                                        predict(mod1,subset(data,Z==0)),
        ps=predict(lmod,subset(data,Z==0),type='response')
    )

    estFun2 <- function(data){
        function(theta){
            c(
                theta[1]-data$ps,
                theta[2]-data$ps^2,
                theta[4]*theta[1]+theta[3]*(1-theta[1])-data$r,
                theta[4]*theta[2]+theta[3]*(theta[1]-theta[2])-data$r*data$ps
            )
        }
    }

    if(justEst){
        bbb <- create_basis(estFun2,dat1)
        mycontrol <- new('geex_control', .root = setup_root_control(start = rep(.5,4)))
        bbb@.control <- mycontrol
        theta <- estimate_GFUN_roots(bbb)$root

        theta <- c(theta,
                   -theta[3],#+mean(data$Y[data$Z==1&data$S==0]),
                   coef(mod1)['S']-theta[4]
                   )
        names(theta) <- c('Epi','Epi2','mu0','mu1','eff0','eff1')

        theta <- c(theta,attr(data,"trueEff"))

        return(theta)
    }

    m_estimate(estFun2,data=dat1,root_control=setup_root_control(start = rep(.5,4)))
}


oneStep <- function(dat){
    estFun <- function(data){

        function(theta){
            ## ols z=1
            xb1 <- with(data,(cbind(1,S,x1,x2)%*%theta[1:4])[,1])

            ## xb for z=0
            xb2 <- with(data,(cbind(x1,x2)%*%theta[3:4])[,1]) ## (b is same in trt groups)

            ## logit z=1
            xb3 <- with(data,(cbind(1,x1,x2)%*%theta[5:7])[,1])
            ps <- plogis(xb3)

            out <- c(
### regression for treatment group
                ifelse(data$Z==1,data$Y-xb1,0),
                ifelse(data$Z==1,data$S*(data$Y-xb1),0),
                ifelse(data$Z==1,data$x1*(data$Y-xb1),data$x1*(data$Y-xb2-theta[9]*ps-theta[8]*(1-ps))),
                ifelse(data$Z==1,data$x2*(data$Y-xb1),data$x2*(data$Y-xb2-theta[9]*ps-theta[8]*(1-ps))),
### logistic regression
                ifelse(data$Z==1,data$S-ps,0),
                ifelse(data$Z==1,data$x1*(data$S-ps),0),
                ifelse(data$Z==1,data$x2*(data$S-ps),0),
### mixture model in control group
                ifelse(data$Z==0,theta[9]*ps+theta[8]*(1-ps)-data$Y+xb2,0),
                ifelse(data$Z==0,theta[9]*ps^2+theta[8]*(ps-ps^2)-(data$Y-xb2)*ps,0),
                theta[10]-(theta[1]-theta[8]),
                theta[11]-(theta[1]+theta[2]-theta[9])
            )


            out
        }
    }



    res <- try(m_estimate(estFun,dat,root_control = setup_root_control(start = rep(0.1,11))))
    if(inherits(res,'try-error'))
        res <- try(m_estimate(estFun,dat,root_control = setup_root_control(start = rep(1,11))))
    if(inherits(res,'try-error'))
        return(rep(NA,4))

    est <- coef(res)
    se <- sqrt(diag(vcov(res)))

    c(eff0=est[10],se0=se[10],eff1=est[11],se1=se[11])
}


oneStepInt <- function(dat){
    estFun <- function(data){

        function(theta){
            ## ols z=1
            xb1 <- with(data,(cbind(1,S,x1,x2)%*%theta[1:4])[,1])

            ## xb for z=0
            xb2 <- with(data,(cbind(x1,x2)%*%theta[12:13])[,1]) ##  (b differs btw trt grps)

            ## logit z=1
            xb3 <- with(data,(cbind(1,x1,x2)%*%theta[5:7])[,1])
            ps <- plogis(xb3)

            c(
### regression for treatment group
                ifelse(data$Z==1,data$Y-xb1,0),
                ifelse(data$Z==1,data$S*(data$Y-xb1),0),
                ifelse(data$Z==1,data$x1*(data$Y-xb1),0),
                ifelse(data$Z==1,data$x2*(data$Y-xb1),0),
### logistic regression
                ifelse(data$Z==1,data$S-ps,0),
                ifelse(data$Z==1,data$x1*(data$S-ps),0),
                ifelse(data$Z==1,data$x2*(data$S-ps),0),
### mixture model in control group
                ifelse(data$Z==0,theta[9]*ps+theta[8]*(1-ps)-data$Y+xb2,0),
                ifelse(data$Z==0,theta[9]*ps^2+theta[8]*(ps-ps^2)-(data$Y-xb2)*ps,0),
                theta[10]-(theta[1]-theta[8]),
                theta[11]-(theta[1]+theta[2]-theta[9]),
                ifelse(data$Z==0,data$x1*(data$Y-xb2-theta[9]*ps-theta[8]*(1-ps)),0),
                ifelse(data$Z==0,data$x2*(data$Y-xb2-theta[9]*ps-theta[8]*(1-ps)),0)
            )
        }
    }

    res <- try(m_estimate(estFun,dat,root_control = setup_root_control(start = rep(0.1,13))))
    if(inherits(res,'try-error'))
        res <- try(m_estimate(estFun,dat,root_control = setup_root_control(start = rep(1,13))))
    if(inherits(res,'try-error'))
        return(rep(NA,4))

    est <- coef(res)
    se <- sqrt(diag(vcov(res)))

    c(eff0=est[10],se0=se[10],eff1=est[11],se1=se[11])
}



### mu00=0
makeDat <- function(n,mu01=0.2,mu10=0,mu11=0.5,b1=1,gumb=FALSE,...){

    x1 <- rnorm(2*n)
    x2 <- rnorm(2*n)
    x3 <- rnorm(2*n)

    x1 <- x1-mean(x1)
    x2 <- x2-mean(x2)
    x3 <- x3-mean(x3)

    psTrue <- plogis(b1*(x1+x2+x3))

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)

    error <- if(gumb) evd::rgumbel(2*n,loc=0,scale=0.16) else rnorm(2*n,0,0.2)
    error <- error-mean(error)

    Y <- 0.5*(x1+x2+x3)+mu01*S+mu10*Z+(mu11-mu01-mu10)*Z*S+error

    dat <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2,Strue=S)
    attr(dat,'trueEffs') <- c(
        S0=mean(Y[Z==1&S==0])-mean(Y[Z==0&S==0]),
        S1=mean(Y[Z==1&S==1])-mean(Y[Z==0&S==1])
    )
    dat

}

sim2Step <- function(n,...){
    dat <- makeDat(n,...)
    twoStep(dat)
}

simStan2step <- function(n,...,se=TRUE){
    dat <- makeDat(n,...)
    two <- if(se) try(oneStep(dat)) else twoStep(dat)
    MLE <- try(mle(dat,se=se))

    if(inherits(two,'try-error') | inherits(MLE, 'try-error')){
      out <- rep(NA, ifelse(se,10,6))
    } else {

        out <- setNames(c(attr(dat,'trueEffs'),two,MLE),
                    if(se){
                        c('true0','true1','mom0','mom0se','mom1','mom1se','mle0','mle0se','mle1','mle1se')
                        } else c('true0','true1','mom0','mom1','mle0','mle1')
             )
    }

    out <- c(out,n,unlist(list(...)))
    #if(any(abs(out)>2)) return(list(out,dat))
    out
}


summ <- function(sss){
    rbind(rowMeans(sss),
          apply(sss,1,sd),
          sqrt(rowMeans((sss-sss[rep(1:2,3),])^2))
          )
    }

simOneBayes <- function(n,...){
    dat <- makeDat(n,...)
    mest <- twoStep(dat)
    BAYES <- bayes(dat)

    list(
        true=attr(dat,'trueEffs'),
        mest=mest,
        bayes=BAYES,
        facs=c(n=n,unlist(list(...)))
    )
}


oneCase <- function(nsim,cl, facs,Bayes=FALSE,se=TRUE){ #n,mu00,mu01,mu10,mu11,gumb,b1,cl){

    print(Sys.time())

    facs$se <- se
  #cat(n,mu00,mu01,mu10,mu11,gumb,b1,'\n',sep=' ')
  clusterExport(cl,'facs',envir=environment())

  time <- system.time(
    res <-
      parLapply(cl, #mclapply( #pbreplicate(
        1:nsim,
        function(i) if(Bayes) try(do.call('simOneBayes',facs)) else try(do.call('simStan2step',facs))
      )
   )
  print(time)
    if(!Bayes){
        resMat <- try(do.call('rbind',res))
        if(!inherits(resMat,'try-error')) res <- resMat
    }
  attr(res,"time") <- time
  res
}


fullsim <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),#sepTs=c(TRUE,FALSE),
                    mu10=c(0,.3),#sepCs=c(TRUE,FALSE),
                    mu11=.3,#effs=c(TRUE,FALSE),
                    gumbs=c(TRUE,FALSE),
                    b1s=c(0,0.2,0.5,1),
                    ext='',
                    se=TRUE,
                    cl=NULL,
                    Bayes=FALSE,
		    start=1
                    ){

    cases=expand.grid(ns,mu01,mu10,mu11,gumbs,b1s)
    names(cases) <- c('n','mu01','mu10','mu11','gumb','b1')

    cat(nrow(cases),' conditions\n')

    save(cases,file=paste0('simResults/cases',ext,'.RData'))

    for(i in start:nrow(cases)){
    	  cat(round(i/nrow(cases)*100),'%\n')
	  facs <- cases[i,]
    	  res <- oneCase(nsim=nsim,cl=cl,facs=facs,se=se,Bayes=Bayes)
	  save(res,facs,file=paste0('simResults/sim',i,ext,'.RData'))
	  }

    return(0)
}


datNewForm <- function(dat){
    X <- as.matrix(dat[,startsWith(names(dat),'x')])
    with(dat,list(
                 Y0=Y[Z==0],
                 Y1=Y[Z==1],
                 S1=S[Z==1],
                 X0=X[Z==0,],
                 X1=X[Z==1,]
             )
         )
}


fullsimJustM <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),#sepTs=c(TRUE,FALSE),
                    mu10=c(0,.3),#sepCs=c(TRUE,FALSE),
                    mu11=.3,#effs=c(TRUE,FALSE),
                    gumbs=c(TRUE,FALSE),
                    b1s=c(0,0.2,0.5,1),
                    ext='',
                    se=TRUE,
                    cl=NULL,
                    Bayes=FALSE,
		    start=1
                    ){

    cases=expand.grid(ns,mu01,mu10,mu11,gumbs,b1s)
    names(cases) <- c('n','mu01','mu10','mu11','gumb','b1')

    if(nsim==0) return(cases)

    cat('% done:')
    list(cases=cases,
         res=lapply(
             1:nrow(cases),
             function(i){
                 cat(round(i/nrow(cases)*100))
                 replicate(nsim,effs(do.call("makeDat",cases[i,])))
             }
         )
         )
}





summarizeSimReg <- function(simRes){
    require(dplyr)
    cases <- simRes$cases%>%
        mutate(
            eff0=mu10,
            eff1=mu11-mu01,
            diff=eff1-eff0
        )

    sss <- simRes$res
    cases <- cbind(cases,
                   est=t(vapply(1:length(sss),function(i) rowMeans(sss[[i]][,1,]),numeric(3))))

    cover <- function(x,b) vapply(rownames(x),function(y) mean((x[y,1,]-2*sqrt(x[y,2,])<b[y])&(x[y,1,]+2*sqrt(x[y,2,])>b[y])),numeric(1))

    cases <- cbind(cases,trueVar=t(vapply(1:length(sss),function(i) apply(sss[[i]][,1,],1,var),numeric(3))))
    cases <- cbind(cases,estVar=t(vapply(1:length(sss),function(i) rowMeans(sss[[i]][,2,]^2),numeric(3))))
    cases <- cbind(cases,coverage=t(vapply(1:length(sss),function(i) cover(sss[[i]],cases[i,c('eff0','eff1','diff')]),numeric(3))))



    cases <- cbind(cases,biasEst=do.call("cbind",setNames(lapply(c('eff0','eff1','diff'),function(x) cases[,paste0('est.',x)]-cases[,x]),c('eff0','eff1','diff'))))
    cases <- cbind(cases,biasVar=do.call("cbind",setNames(lapply(c('eff0','eff1','diff'),function(x) cases[,paste0('estVar.',x)]-cases[,paste0('trueVar.',x)]),c('eff0','eff1','diff'))))

#    cases <- cbind(cases,coverage=t(vapply(1:length(sss),setNames(lapply(c('eff0','eff1','diff'),function(x) cover(sss,paste0('estVar.',x)]-cases[,paste0('trueVar.',x)]),c('eff0','eff1','diff'))))

    cases
}

