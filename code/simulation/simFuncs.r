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
    fit <- stan('code/ps.stan',data=sdat,...)

    summary(fit, par=c('eff0','eff1','effDiff'))$summary
}



### mu00=0
makeDat <- function(n,mu01=0.2,mu10=0,mu11=0.5,b1=1,errDist='norm',intS=FALSE,intZ=FALSE,...){

    x1 <- rnorm(2*n)
    x2 <- rnorm(2*n)
    x3 <- if(errDist=='norm') rnorm(2*n) else runif(2*n,-sqrt(12)/2,sqrt(12)/2)

    x1 <- x1-mean(x1)
    x2 <- x2-mean(x2)
    x3 <- x3-mean(x3)

    psTrue <- plogis(b1*(x1+x3)-b1*x2)

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)

    error <- if(errDist=='norm') rnorm(2*n,0,sqrt(0.5))
             else if(errDist=='mix') c(rnorm(3*n/2,-1/3,sqrt(1/6)),rnorm(n/2,1,sqrt(1/6)))
             else runif(2*n,-sqrt(6)/2,sqrt(6)/2)
    error <- error-mean(error)

    Yc <- sqrt(1/6)*(x1+x2+x3)+mu01*S+error
    if(intS) Yc <- Yc-sqrt(1/6)/4*(x1+x2)+S*sqrt(1/6)/2*(x1+x2)

    #x0 <- mean(x1[S==0]+x2[S==0])

    Yt <- Yc+mu10+(mu11-mu01-mu10)*S
    if(intZ) Yt <- Yt+sqrt(1/6)/2*(x1+x2)

    Y <- ifelse(Z==1,Yt,Yc)

    ## if(!norm){
    ##     Y <- round(Y)
    ##     #Y[Y< -4] <- 4
    ##     #Y[Y> 4] <- 4
    ## }

    dat <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2,Strue=S)
    attr(dat,'trueEffs') <- c(
        S0=mean(Yt[S==0])-mean(Yc[S==0]),
        S1=mean(Yt[S==1])-mean(Yc[S==1])
    )
    attr(dat,'facs') <- unlist(as.list(match.call())[-1])

    dat

}


simOneBayes <- function(dat){

    mest <- effs(dat)
    BAYES <- bayes(dat,chains=2,iter=3000,warmup=1000)

    list(
        true=attr(dat,'trueEffs'),
        mest=mest,
        bayes=BAYES,
        facs=attr(dat,'facs')
    )
}


oneCase <- function(nsim,ext,ncores, facs){ #n,mu00,mu01,mu10,mu11,gumb,b1,cl){

    print(Sys.time())

    datasets <- mclapply(1:nsim,function(i) do.call("makeDat",facs),mc.cores=ncores)
    save(datasets,file=paste0('simData/dat',ext,'.RData'))

    time <- system.time(
        res <-
            mclapply(
                datasets,
                function(dat) try(simOneBayes(dat)),
                mc.cores=ncores
            )
    )
    print(time)

  attr(res,"time") <- time
  res
}

fullsim <- function(nsim,
                    ns=c(100,500,1000),
                    mu01=c(0,.3),#sepTs=c(TRUE,FALSE),
                    mu10=c(0,.3),#sepCs=c(TRUE,FALSE),
                    mu11=.3,#effs=c(TRUE,FALSE),
                    errDist=c('norm','mix','unif'),
                    b1s=c(0,0.2,0.5),
                    ext='',
                    se=TRUE,
                    ncores=parallel::detectCores(),
		    start=1
                    ){

    cases=expand.grid(ns,mu01,mu10,mu11,errDist,b1s,intS=c(TRUE,FALSE),intZ=c(TRUE,FALSE))
    names(cases) <- c('n','mu01','mu10','mu11','norm','b1')

    cat(nrow(cases),' conditions\n')

    save(cases,file=paste0('simResults/cases',ext,'.RData'))

    for(i in start:nrow(cases)){
    	  cat(round(i/nrow(cases)*100),'%\n')
	  facs <- cases[i,]
    	  res <- oneCase(nsim=nsim,ext=paste0(i,ext),
                         ncores=ncores,facs=facs)
	  save(res,facs,file=paste0('simResults/sim',i,ext,'.RData'))
    }

    return(0)
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
                 replicate(nsim,effs(do.call("makeData",cases[i,])))
             }
         )
         )
}




xSim <- function(){
     n <- 10000
 mu01=0.2
 mu10=0
 mu11=0.5
 b1=1
 errDist='norm'
 x1 <- rnorm(2*n)
 x2 <- rnorm(2*n)
 x3 <- if(errDist=='norm') rnorm(2*n) else runif(2*n,-sqrt(12)/2,sqrt(12)/2)
 x1 <- x1-mean(x1)
 x2 <- x2-mean(x2)
 x3 <- x3-mean(x3)
 psTrue <- plogis(b1*(x1+x3)-b1*x2)
 S <- rbinom(2*n,1,psTrue)
 Z <- rep(c(1,0),n)
 x12 <- x1+x2
     c(mean(x12),
       mean(x12[S==1]),
       mean(x12[S==0]))
}
