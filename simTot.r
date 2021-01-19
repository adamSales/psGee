library(tidyverse)
library(rstan)
library(geex)
library(parallel)
library(pbapply)

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





makeDat <- function(n,mu00=0,mu01=0.2,mu10=0,mu11=0.5,b1=1,gumb=FALSE){

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

    Yc <- 0.5*(x1+x2+x3)+mu01*S
    Y <- Yc+(mu11-mu01)*S*Z+error


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

simStan2step <- function(n,...){
    dat <- makeDat(n,...)
    two <- twoStep(dat)
    MLE <- mle(dat)

    out <- setNames(c(two[c('S0','S1','eff0','eff1')],MLE),
             c('true0','true1','mom0','mom1','mle0','mle1')
             )
    #if(any(abs(out)>2)) return(list(out,dat))
    out
}

bias <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    colMeans(cbind(sss1[,3:6]-sss1[,c(1,2,1,2)],sss2))
}

rmse <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    c(
        sqrt(colMeans((sss1[,3:6]-sss1[,c(1,2,1,2)])^2)),
        colMeans(sss2)
    )
}


summ <- function(sss){
    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]
    rbind(colMeans(sss1),
          apply(sss1,2,sd),
          sqrt(row1Means((sss-sss[rep(1:2,3),])^2))
          )
    }


fullsim <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),#sepTs=c(TRUE,FALSE),
                    mu10=c(0,.3),#sepCs=c(TRUE,FALSE),
                    mu11=c(.3,.3),#effs=c(TRUE,FALSE),
                    gumbs=c(TRUE,FALSE),
                    b1s=c(0,0.2,0.5,1),
                    cl=NULL
                    ){

    cases=expand.grid(ns,mu01,mu10,mu11,gumbs,b1s)
    names(cases) <- c('n','mu01','mu10','mu11','gumb','b1')

    cases%>%
        rowwise()%>%
        mutate(res=list(pbreplicate(nsim,simStan2step(n=n,mu00=0,mu01=mu01,mu10=mu10,mu11=mu11,gumb=gumb,b1=b1),cl=cl)))
}


## if(!exists('cl')){
##     cl <- makeCluster(6)
##     clusterEvalQ(cl,library(tidyverse))
##     clusterEvalQ(cl,library(rstan))
##     clusterEvalQ(cl,library(geex))

##     clusterExport(cl,c('simStan2step','makeDat','mle','twoStep','stanMod'))
## }

res <- fullsim(500,cl=50)
save(res,file='simulation.RData')

 for(i in 1:99){
 load(paste0('simResults/sim',i,'.RData'))
 resList[[i]] <- res
}

biases <- sapply(resList,bias)
rmses <- sapply(resList,rmse)


plotRes <- function(rrr,fac,S){
    www <- (max(rrr[fac,])-min(rrr[fac,]))*0.02
    plot(rrr[fac,]-www,rrr[paste0('mom',S),],ylim=range(c(rrr[paste0('mom',S),],rrr[paste0('mle',S),])),xlim=c(min(rrr[fac,])-2*www,max(rrr[fac,])+2*www),pch=16,xlab=fac)
    points(rrr[fac,]+www,rrr[paste0('mle',S),],col='red',pch=16)
    legend('topright',legend=c('MOM','MLE'),col=c('black','red'),pch=16)
    abline(h=0,lty=2)
}


loadRes <- function(){
    load('simResults/cases.RData')

    results <- list()
    for(i in 1:nrow(cases)){
        if(i %% 10==0) cat(round(i/nrow(cases)*100), '% ')
        load(paste0('simResults/sim',i,'.RData'))
        stopifnot(identical(facs,cases[i,]))
        results[[i]] <- as.data.frame(res)
        results[[i]]$n <- facs$n
        rm(res,facs)
    }

    do.call('rbind',results)

}


results%>%
    filter(n==1000,gumb==0,mu01==0.3,mu10==0.3,mu11==0.3)%>%
    group_by(b1,mu01,mu10,mu11)%>%
    summarize(across(c(true0,true1),mean),across(c(mom0,mle0),~mean(.-true0)),across(c(mom1,mle1),~mean(.-true1)))
