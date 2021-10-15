library(sandwich)

### I want a bunch of attributes for the VCV matrix, but I don't
### want them to print every time!
print.vcv <- function(x,...) print(x[1:NROW(x),1:NCOL(x)],...)

### convenience function for functions that return lists
### why cant R be more like python in this regard??
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

### computes/returns regression models
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

### estimating equations for stage 1 model (PS)
EFps <- function(psMod,data){
    efPStmp <- estfun(psMod)
    efPS <- matrix(0,nrow(data),ncol(efPStmp))
    efPS[as.numeric(rownames(efPStmp)),] <- efPStmp
    efPS
}

### estimates bread and meat for sandwich variance
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

### combining smaller matrices into bigger one
bigM <- function(m11,m12,m22)
    rbind(
        cbind(m11,m12),
        cbind(t(m12),m22)
    )

### computes sandwich vcov for regressions
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



### wrapper function for estimating regressions+vcov
est <- function(data,covForm=~x1+x2,int=FALSE){

    Attach(pointEst(data=data,covForm=covForm,int=int))

    vcv <- vcvPS(psMod,outMod)

    list(outMod=outMod,psMod=psMod,vcv=vcv)
}

### estimates effects of interest, starting from est() output
effsFromFit <- function(ests){
        estimates <- with(as.list(coef(ests$outMod)),
                      list(
                          eff1=Z+`Z:Sp`,
                          diff=`Z:Sp`))
    vcv <- ests$vcv
    ddd <- diag(vcv)
    vars <- list(
        eff0 <- ddd['Z'],
        eff1 <- ddd['Z']+ddd['Z:Sp']+2*vcv['Z','Z:Sp'],
        diff <- ddd['Z:Sp']
    )
    cbind(
        estimates=unlist(estimates),
        SE=sqrt(unlist(vars))
    )
}

### estimates effects of interest, starting from data
effs <- function(data,covForm=~x1+x2,int=FALSE){
    ests <- est(data=data,covForm=covForm,int=int)

    effsFromFit(ests)
}

### estimates lower left of A (bread) matrix in A^{-1}BA^{-t}
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

