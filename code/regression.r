library(sandwich)

print.vcv <- function(x,...) print(x[1:NROW(x),1:NCOL(x)],...)

Attach <- function(.list,.names=names(.list),env=environment()){
    if(is.null(.names))
        .names <- paste0(
            as.character(
                as.list(match.call())$.list)[1],
            seq(length(.list)))

    for(i in 1:length(.list)){
        assign(.names[i],.list[[i]], envir=parent.frame())
    }

}

pointEst <- function(data,covForm=~x1+x2,int=FALSE){
    psMod <- glm(update(covForm,S~.),data=data,family=binomial,subset=Z==1)

    ps <- predict(psMod,data,type='response')
    if('Sp'%in%names(data)) warning('replacing Sp')
    data <- within(data,Sp <- ifelse(Z==1,S,ps))

    outMod <- lm(
        update(covForm,if(int) Y~Z*(Sp+.) else Y~Z*Sp+.),
        data=data
    )

    list(psMod=psMod,outMod=outMod)
}

EFps <- function(psMod,data){
    efPStmp <- estfun(psMod)
    efPS <- matrix(0,nrow(data),ncol(efPStmp))
    efPS[as.numeric(rownames(efPStmp)),] <- efPStmp
    efPS
}

sandwichMats <- function(psMod,outMod,int=any(grepl(":x",names(coef(outMod))))){
    data <- model.frame(outMod)
    ### estimating equations
    efPS <- EFps(psMod,data)
    efOut <- estfun(outMod)

    out <- list(
        a11inv = bread(psMod)/sum(data$Z),
        a22inv = bread(outMod)/nrow(data),
        a21 = A21(psMod,outMod),
        b11=meat(psMod)*nrow(model.frame(psMod)),
        b22 = meat(outMod)*nrow(data)#,adjust=TRUE)
    )
    out$b12 <-
        if(int){
            matrix(0,nrow(b11),ncol(b22))
        } else crossprod(efPS,efOut) #/sum(data$Z==0)

    out
}

bigM <- function(m11,m12,m22)
    rbind(
        cbind(m11,m12),
        cbind(t(m12),m22)
    )


vcvPS <- function(psMod,outMod,int=any(grepl(":x",names(coef(outMod))))){
    Attach(sandwichMats(psMod=psMod,outMod=outMod,int=int))

    A <- rbind(
        cbind(solve(a11inv),matrix(0,nrow(a11inv),ncol(a22inv))),
        cbind(a21,solve(a22inv))
    )
    colnames(A) <- rownames(A)
    ##bigM(solve(a11inv),t(a21),solve(a22inv))

    B <- bigM(b11,b12,b22)

    n <- nrow(model.frame(outMod))

    vcvFull <- solve(A)%*%B%*%t(solve(A))
    main <- a22inv%*%b22%*%t(a22inv)
    vcv <- vcvFull[-(1:nrow(b11)),-(1:nrow(b11))]
    dimnames(vcv) <- dimnames(main)
    correction <- vcv-main

    attr(vcv,"full") <- vcvFull
    attr(vcv,"bread") <- A
    attr(vcv,"meat") <- B
    attr(vcv,"main") <- main
    attr(vcv,"correction") <- correction

    class(vcv) <- "vcv"

    vcv
}

estEff1bread <- function(dat){
    Attach(pointEst(dat))
    vcv <- vcvPS(psMod,outMod)
    c(est=coef(outMod)['Z'],
      var=diag(solve(attr(vcv,'bread')))['Z']/(nrow(model.frame(outMod))^2))
}


oldVCV <- function(psMod,outMod,int=any(grepl(":x",names(coef(outMod))))){

    data <- model.frame(outMod)
### estimating equations

    efPStmp <- estfun(psMod)
    efPS <- matrix(0,nrow(data),ncol(efPStmp))
    efPS[as.numeric(rownames(efPStmp)),] <- efPStmp
    rm(efPStmp)

    efOut <- estfun(outMod)

    a11inv <- bread(psMod)/sum(data$Z)
    a22inv <- bread(outMod)/nrow(data)
    a21 <- A21(psMod,outMod)
    b22 <- meat(outMod)*nrow(data)#,adjust=TRUE)

    if(!int) b12 <- crossprod(efPS,efOut) #/sum(data$Z==0)

    main <- sandwich(outMod)#,bread.=a22inv*nrow(data),meat.=b22)
    vcv1 <- sandwich(psMod)#,adjust=TRUE)

    correction <- a21%*%vcv1%*%t(a21)
    if(!int) correction <- correction-
                 a21%*%a11inv%*%b12-
                 t(b12)%*%t(a11inv)%*%t(a21)

    correction <- a22inv%*%correction%*%t(a22inv)#/nrow(data)#sum(data$Z==0)

    vcv <- main+correction
    attr(vcv,"main") <- main
    attr(vcv,"correction") <- correction
    vcv
}


est <- function(data,covForm=~x1+x2,int=FALSE){

    Attach(pointEst(data=data,covForm=covForm,int=int))

    vcv <- vcvPS(psMod,outMod)

    list(outMod=outMod,psMod=psMod,vcv=vcv)
}


effs <- function(data,covForm=~x1+x2,int=FALSE){
    ests <- est(data=data,covForm=covForm,int=int)

    estimates <- with(as.list(coef(ests$outMod)),
                      list(
        #mu00=`(Intercept)`
        #mu01=mu00+Sp
                          eff0=Z,
        #mu10=mu00+eff0
        #mu11=mu00+Z+Sp+`Z:Sp`
                          eff1=Z+`Z:Sp`,
                          diff=`Z:Sp`))
    vcv <- ests$vcv
    ddd <- diag(vcv)
    vars <- list(
        #mu00=ddd['(Intercept)']
        #mu01=mu00+vcv['Sp','Sp']+2*vcv['(Intercept)','Sp']
        eff0 <- ddd['Z'],
        #mu10 <- mu00+eff0+2*vcv['(Intercept)','Z']
        #mu11 <- mu00+eff0+ddd['Sp']+ddd['Z:Sp']+
         #   2*(sum(vcv['(Intercept)',c('Z','Sp','Z:Sp')])+
         #      sum(vcv['Z',c('Sp','Z:Sp')])+vcv['Sp','Z:Sp'])
        eff1 <- ddd['Z']+ddd['Z:Sp']+2*vcv['Z','Z:Sp'],
        diff <- ddd['Z:Sp']
    )
    cbind(
        estimates=unlist(estimates),
        SE=sqrt(unlist(vars))
    )
}

A21 <- function(psMod,outMod){
    Z <- model.frame(outMod)$Z
    Y0 <- model.frame(outMod)$Y[Z==0]
    aX <- predict(psMod,model.frame(outMod),type='link')[Z==0]
    q <- family(psMod)$mu.eta(aX)      ## dp/dalpha=qX'
    u <- 2*family(psMod)$linkinv(aX)*q   ## dp^2/dalpha=uX'

    X0 <- model.matrix(update(formula(psMod),Sp~.),data=model.frame(outMod))[Z==0,]


    est <- list(beta=coef(outMod)[colnames(X0)[-1]],
                eta=coef(outMod)['Sp'],
                mu0=coef(outMod)['(Intercept)']
                )

    a21=rbind(
       est$eta*t(q)%*%X0,
       -t((Y0-est$mu0-X0[,-1]%*%est$beta)*q-est$eta*u)%*%X0,
        est$eta*t(X0[,-1])%*%diag(as.vector(q))%*%X0
    )#/nrow(X0)
    rownames(a21)[1:2] <- c('(Intercept)','Sp')

    t(vapply(
        names(coef(outMod)),
        function(n)
            if(n %in% rownames(a21)) a21[n,] else rep(0,ncol(a21)),
        numeric(ncol(a21)))
      )
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





### mu00=0
makeData <- function(n,mu01=0.2,mu10=0,mu11=0.5,b1=1,gumb=FALSE,...){

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


    data <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2,Strue=S)
    attr(data,'trueEffs') <- c(
        S0=mean(Y[Z==1&S==0])-mean(Y[Z==0&S==0]),
        S1=mean(Y[Z==1&S==1])-mean(Y[Z==0&S==1])
    )
    data

}

summarizeSim <- function(simRes){
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

    cases <- cbind(cases,trueVar=t(vapply(1:length(sss),function(i) apply(sss[[i]][,1,],1,var),numeric(3))))
    cases <- cbind(cases,estVar=t(vapply(1:length(sss),function(i) rowMeans(sss[[i]][,2,]^2),numeric(3))))

    cases <- cbind(cases,biasEst=do.call("cbind",setNames(lapply(c('eff0','eff1','diff'),function(x) cases[,paste0('est.',x)]-cases[,x]),c('eff0','eff1','diff'))))
    cases <- cbind(cases,varEst=do.call("cbind",setNames(lapply(c('eff0','eff1','diff'),function(x) cases[,paste0('estVar.',x)]-cases[,paste0('trueVar.',x)]),c('eff0','eff1','diff')))

    cases
}
