library(sandwich)


### step 1: S~X logit (for Z=1)
### step 2: Y~X+S OLS (for Z=1)
### step 3: estimating equations:
## Y-(p mu1+(1-p)mu0+xBeta)
## Yp-p(p mu1+(1-p)mu0+xBeta)
## XY-X(p mu1+(1-p)mu0+xBeta)



pointEst <- function(Y0,S1,X1,X0,...){
    ## step 1
    Smod <- glm(S1~X1,family=binomial)

    ps <- plogis(coef(Smod)[1]+X0%*%coef(Smod)[-1])

    ## system of eq's
    ## Y=mu1*p+mu0*(1-p)+Xbeta
    ## Yp=mu1*p^2+mu0*p*(1-p)+beta'X*p
    ## XY=mu1*p*X+mu0*(1-p)*X+beta*X*X

    n <- length(ps)
    A <- lapply(1:n,function(i)
        rbind(
            c(mu1=ps[i],mu0=1-ps[i],X0[i,]),
            c(ps[i]^2,ps[i]*(1-ps[i]),ps[i]*X0[i,]),
            cbind(X0[i,]*ps[i],X0[i,]*(1-ps[i]),X0[i,]%*%t(X0[i,]))
        )
        )

    B <- lapply(1:n, function(i) c(Y0[i],ps[i]*Y0[i],Y0[i]*X0[i,]))

    est <- solve(Reduce("+",A),Reduce("+",B))

    list(est=est,ef=do.call('rbind',lapply(1:n,function(i) t(B[[i]]-A[[i]]%*%est))))
}


A21 <- function(logit1,Y0,mu01,mu00,beta0,X0){
### derivative of step 3 EEs as functions of logit coefs
    p <- plogis(coef(logit1)[1]+X0%*%coef(logit1)[-1])
    pq <- p*(1-p)
    Xtilde <- cbind(1,X0)

    mat <- rbind(
        -t(pq)%*%Xtilde*(mu01-mu00),
        t(pq*(Y0-X0%*%beta0-mu00))%*%Xtilde-2*(mu01-mu00)*t(p*pq)%*%Xtilde,
        -(mu01-mu00)*t(sweep(X0,1,pq,"*"))%*%Xtilde
    )
    mat
}

A22 <- function(logit1,X0){
    p <- plogis(coef(logit1)[1]+X0%*%coef(logit1)[-1])

    -rbind(
         c(sum(p),sum(1-p),colSums(X0)),
         c(sum(p^2),sum(p^2-p),t(p)%*%X0),
         cbind(t(X0)%*%p,t(X0)%*%(p*(1-p)),crossprod(X0))
     )
}

ef2 <- function(Y0,X0,mu1,mu0,logit1,beta0){
    p <- plogis(coef(logit1)[1]+X0%*%coef(logit1)[-1])
    xb <- X0%*%beta0
    cbind(
        Y0-p*mu1-(1-p)*mu0-xb,
        Y0*p-p^2*(mu1-mu0)-p*mu0-p*xb,
        sweep(X0,1,Y0-p*mu1-(1-p)*mu0-xb,"*")
    )
}

B12 <- function(logit1,Y0,X0,mu1,mu0,beta0){
    EF1 <- estfun(logit1)
    EF0 <- ef2(Y0,X0,mu1,mu0,logit1,beta0)
    t(rbind(EF1,matrix()%*%

B22 <- function(Y0,X0,mu1,mu0,logit1,beta0)
    crossprod(ef2(Y0,X0,mu1,mu0,logit1,beta0))

vcv2 <- function(dat,est,ef,logit1){
    a11inv <- bread(logit1)
    a21 <- A21(logit1=logit1,Y0=dat$Y0,est$mu1,est$mu0,est$beta0,dat$X0)
    a22 <- A22(logit1,dat$X0)
    b11 <- meat(logit1)
    b12 <- B12(logit1,dat$Y0,dat$X0,est$mu1,est$mu0,est$beta0)
    b22 <- B22(dat$Y0,dat$X0,est$mu1,est$mu0,logit1,est$beta0)

    n <- length(dat$Y0)
    solve(a22)%*%(b22+a21%*%a11inv%*%b11%*%t(a11inv)%*%t(a21))/n

}
