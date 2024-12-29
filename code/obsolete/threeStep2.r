library(sandwich)


Est <- function(Y0,X0,ps,...){

    ## system of eq's
    ## Yi-mu0-p_i*eta-X_i'beta
    ## piYi-p_i*mu0-pi^2*eta-p_i*X_i'beta
    ## X_iY_i-X_i*mu0-X_ip_i*eta-X_i*X_i'beta

    n <- length(ps)
    A <- lapply(1:n,function(i)
        rbind(
            c(mu0=1,eta=ps[i],X0[i,]),
            c(ps[i],ps[i]^2,ps[i]*X0[i,]),
            cbind(X0[i,],X0[i,]*ps[i],X0[i,]%*%t(X0[i,]))
        )
        )

    B <- lapply(1:n, function(i) c(Y0[i],ps[i]*Y0[i],Y0[i]*X0[i,]))

    est <- solve(Reduce("+",A),Reduce("+",B))
    estList <- list(mu0=est['mu0'],eta=est['eta'],beta=est[!names(est)%in%c('mu0','eta')])

    list(est=estList,ef=do.call('rbind',lapply(1:n,function(i) t(B[[i]]-A[[i]]%*%est))))
}

## Y - mu0-eta*ps-x'beta
## psY-ps mu0-eta*ps^2-ps x'beta
## XY-Xmu0-eta*ps*X-XX'beta


A22 <- function(ps,X0){
    rbind(
        c(length(ps),sum(ps),colSums(X0)),
        c(sum(ps),sum(ps^2),t(ps)%*%X0),
        cbind(colSums(X0),t(X0)%*%ps,crossprod(X0))
    )/length(ps)
}

A21 <- function(smod,X0,Y0,est){
    aX <-cbind(1,X0)%*%coef(smod)
    q <- family(smod)$mu.eta(aX)      ## dp/dalpha=qX'
    u <- 2*family(smod)$linkinv(aX)*q   ## dp^2/dalpha=uX'
    X0tilde <- cbind(1,X0)

    rbind(
        t(q)%*%X0tilde,
        -t((Y0-est$mu0-X0%*%est$beta)*q-est$eta*u)%*%X0tilde,
        est$eta*t(X0)%*%diag(as.vector(q))%*%X0tilde
    )/nrow(X0)
}

VCV <- function(a21,a22,vcv1,b22,n){
    a22inv <- solve(a22)

    main <- a22inv%*%b22%*%t(a22inv)/n
    correction <- a22inv%*%a21%*%vcv1%*%t(a21)%*%t(a22inv)/n
    out <- main+correction
    attr(out,"main") <- main
    attr(out,"correction") <- correction
    out
}

estTot <- function(Y0,X0,S1,X1,...){
    n <- length(Y0)

### stage 1: logit model
    smod <- glm(S1~X1,family=binomial)
    ps <- plogis(cbind(1,X0)%*%coef(smod)) ## principal scores
### stage 2: point estimate and estimating eq's
    est <- Est(Y0=Y0,X0=X0,ps=ps)
### A matrices
    #a11inv <- bread(smod)
    a22 <- A22(ps=ps,X0=X0)
    a21 <- A21(smod=smod,X0=X0,Y0=Y0,est=est$est)

### B matrices
    #b11 <- meat(smod,adjust=TRUE)
    b22 <- crossprod(est$ef)/(n-nrow(a22))

    vcv <- VCV(
               a21=a21,
               a22=a22,
               vcv1=sandwich(smod,adjust=TRUE),
               b22=b22,
        n=n)

    colnames(vcv) <- rownames(vcv) <-
        strsplit(names(unlist(est$est)),"\\.")|>
        vapply(function(x) x[2],'a')

    list(est=est$est,vcv=vcv)
}

effectFit <- function(dat){

    ols1 <- lm(Y~S+x1+x2,data=dat,subset=Z==1)
    vcv1 <- sandwich(ols1,adjust=TRUE)

    ests2 <- do.call("estTot",datNewForm(dat))

    list(ols1=ols1,vcv1=vcv1,ests2=ests2)
}

effectCompute <- function(ols1,ests2){

    ## compute effects of interest w/ SE
    estimates <- list()
    estimates <- within(estimates,{
        mu00 <- ests2$est$mu0
        diff0 <- ests2$est$eta
        mu01 <- mu00+diff0
        mu10 <- coef(ols1)['(Intercept)'] ### x1 & x2 are centered
        diff1 <- coef(ols1)['S']
        mu11 <- mu10+diff1
        eff0 <- mu10-mu00
        eff1 <- mu11-mu01
        diff <- diff1-diff0
    })

    unlist(lapply(estimates,unname))
}


SEcompute <- function(vcv1,vcv2){
    SEs <- list()
    vars <- list()
    vars <- within(vars,{
        mu00 <- vcv2['mu0','mu0']
        diff0 <- vcv2['eta','eta']
        mu01 <- mu00+diff0+2*vcv2['mu0','eta']
        mu10 <- vcv1['(Intercept)','(Intercept)']
        diff1 <- vcv2['eta','eta']
        mu11 <- mu00+diff0+2*vcv2['mu0','eta']
        eff0 <- mu00+mu10
        eff1 <- mu01+mu11
        diff <- diff0+diff1
    }

    )
    sqrt(unlist(lapply(vars,unname)))
}

effects <- function(dat){
    fit <- effectFit(dat)
    ests <- effectCompute(fit$ols1,fit$ests2)
    ses <- SEcompute(fit$vcv1,fit$ests2$vcv)
    cbind(ests,ses)
    }

### quick simulation to see if code is working
simTest <- function(n=1000,B=1000){
    source('code/simFuncs.r')
    sss <<- replicate(B,
              makeDat(1000)|>
              datNewForm()|>
              {\(d){do.call("estTot",d)}}()|>
              {\(x){cbind(unlist(x[[1]]),x[[2]])}}()
    )

    ests <- sss[,1,]
    vcvs <- sss[,-1,]
    print("Point estimate bias:")
    bias <- sweep(ests,1,c(mu0=0,eta=0.2,x1=0.5,x2=0.5))
    rownames(bias) <- c('mu0','eta','x1','x2')
    ttests <- apply(bias,1,t.test)
    sigBias <- sapply(ttests,function(x) x$p.value)
    print(rbind(rowMeans(bias),sigBias))

    print("Empirical covariance matrix of estimates:")
    print(cov(t(ests)))
    print("Average estimated covariance matrix:")
    print(apply(vcvs,1:2,mean))
    invisible(sss)
}
