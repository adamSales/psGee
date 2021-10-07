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

    a <- rbind(
        c(mu1=sum(ps), mu0=sum(1-ps), colSums(X0)),
        c(mu1=sum(ps^2), mu0=sum(ps*(1-ps)), t(ps)%*%X0),
        cbind(t(X0)%*%ps,t(X0)%*%(1-ps),crossprod(X0))
    )

    b <- c(sum(Y0),sum(ps*Y0),t(Y0)%*%X0)

    solve(a,b)
}

bdiag <- function(m1,m2)
    cbind(
        rbind(m1,matrix(0,nrow(m1),ncol(m2))),
        rbind(matrix(0,nrow(m2),ncol(m1)),m2)
        )


A11 <- function(logit1,ols1)
    bdiag(bread(logit1),bread(ols1))

meat12 <- function(logit1,ols1,adjust=TRUE){

    psi <- cbind(estfun(logit1),estfun(ols1))
    k <- NCOL(psi)
    n <- NROW(psi)
    rval <- crossprod(as.matrix(psi))/n
    if (adjust)
        rval <- n/(n - k) * rval
    rownames(rval) <- colnames(rval) <- colnames(psi)
    rval
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
    ## append 0s for first OLS regression
    cbind(mat,matrix(0,nrow(mat),ncol(X)+2))
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

B12 <- function(logit1,ols1,Y0,X0,mu1,mu0,beta0)
    t(cbind(estfun(logit1),estfun(ols1)))%*%ef2(Y0,X0,mu1,mu0,logit1,beta0)

B22 <- function(Y0,X0,mu1,mu0,logit1,beta0)
    crossprod(ef2(Y0,X0,mu1,mu0,logit1,beta0))

vcv2 <- function(dat,est,logit1,ols1){
    a11inv <- A11(logit1,ols1)
    a21 <- A21(logit1,dat$Y0,est$mu01,est$mu00,est$beta0,dat$X0)
    a22 <- A22(logit1,dat$X0)
    b11 <- cbind(meat(logit1),meat(ols1))
    b12 <- B12(logit1,ols1,dat$Y0,dat$X0,est$mu1,est$mu0,est$beta0)
    b22 <- B22(dat$Y0,dat$X0,est$mu1,est$mu0,logit1,est$beta0)

    n <- length(dat$Y)
    solve(a22)%*%(b22-a21%*%a11inv%*%b12-t(b12)%*%t(a11inv)%*%t(a21)+a21%*%a11inv%*%b11%*%t(a11inv)%*%t(a21))/n

