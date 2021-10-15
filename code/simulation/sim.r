library(tidyverse)
library(rstan)
library(geex)
library(parallel)
library(pbapply)


source('mom2step.r')
source('momStan.r')



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


summ <- function(sss){
    rbind(rowMeans(sss),
          apply(sss,1,sd),
          sqrt(rowMeans((sss-sss[rep(1:2,3),])^2))
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


if(!exists('cl')){
    cl <- makeCluster(6)
    clusterEvalQ(cl,library(tidyverse))
    clusterEvalQ(cl,library(rstan))
    clusterEvalQ(cl,library(geex))

    clusterExport(cl,c('simStan2step','makeDat','mle','twoStep','stanMod'))
}

res <- fullsim(5)
#save(res,file='simulation.RData')
