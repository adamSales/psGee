### these are basically functions I wrote to find/diagnose errors and confirm estimates
PS <- function(S1,X1,X0,...){
 ## step 1
    Smod <- glm(S1~X1,family=binomial)

    plogis(coef(Smod)[1]+X0%*%coef(Smod)[-1])
}

pointEst <- function(dat){
    ps <- with(dat,PS(S1,X1,X0))
    Est(dat$Y0,dat$X0,ps)$est
}

estVCV2step <- function(dat){
    ps <- do.call("PS",dat)
    n <- length(ps)
    est <- Est(Y0=dat$Y0,X0=dat$X0,ps=ps)
    a22 <- A22(ps=ps,X0=dat$X0)
    a22inv <- solve(a22/n)
    b22 <- crossprod(est[[2]])/n
    VCV <- a22inv%*%b22%*%t(a22inv)/n
    cbind(est=est[[1]],vcv=VCV)
}


simSimp <- function(n){
    ##mu00=0
    ##mu01=0.2

    x1 <- rnorm(n)
    x2 <- rnorm(n)

    x1 <- x1-mean(x1)
    x2 <- x2-mean(x2)

    ps <- plogis(x1+x2)

    S <- rbinom(n,1,ps)

    error <- rnorm(n,0,0.2)
    error <- error-mean(error)

    Y <- 0.5*(x1+x2)+0.2*S+error

    X0 <- cbind(x1,x2)

    eee <- Est(Y,X0,ps)

    a22 <- A22(ps,X0)/n
    a22inv <- solve(a22)
    b22 <- crossprod(eee[[2]])/n


    cbind(est=eee[[1]],vcv=a22inv%*%b22%*%t(a22inv)/n)
}

estFun2 <- function(data){
    function(theta){
        c(
            theta[1]-data$ps,
            theta[2]-data$ps^2,
            theta[4]*theta[1]+theta[3]*(1-theta[1])+theta[5]*data$x1+theta[6]*data$x2-data$Y,
            theta[4]*theta[2]+theta[3]*(theta[1]-theta[2])+theta[5]*theta[1]*data$x1+theta[6]*theta[1]*data$x2-data$Y*data$ps,
            theta[4]*theta[1]*data$x1+theta[3]*(1-theta[1])*data$x1+theta[5]*data$x1^2+theta[6]*data$x2*data$x1-data$Y*data$x1,
            theta[4]*theta[1]*data$x2+theta[3]*(1-theta[1])*data$x2+theta[5]*data$x1*data$x2+theta[6]*data$x2^2-data$Y*data$x2
        )
    }
}

estFun3 <- function(data){
    function(theta){
        c(
            theta[1]-data$ps,
            theta[2]-data$ps^2,
            theta[4]*theta[1]+theta[3]+theta[5]*data$x1+theta[6]*data$x2-data$Y,
            theta[4]*theta[2]+theta[3]*theta[1]+theta[5]*theta[1]*data$x1+theta[6]*theta[1]*data$x2-data$Y*data$ps,
            theta[4]*theta[1]*data$x1+theta[3]*data$x1+theta[5]*data$x1^2+theta[6]*data$x2*data$x1-data$Y*data$x1,
            theta[4]*theta[1]*data$x2+theta[3]*data$x2+theta[5]*data$x1*data$x2+theta[6]*data$x2^2-data$Y*data$x2
        )
    }
}

estFun4 <- function(data){
    function(theta){
        c(
            theta[2]*data$ps+theta[1]+theta[3]*data$x1+theta[4]*data$x2-data$Y,
            theta[2]*data$ps^2+theta[1]*data$ps+theta[3]*data$ps*data$x1+theta[4]*data$ps*data$x2-data$Y*data$ps,
            theta[2]*data$ps*data$x1+theta[1]*data$x1+theta[3]*data$x1^2+theta[4]*data$x2*data$x1-data$Y*data$x1,
            theta[2]*data$ps*data$x2+theta[1]*data$x2+theta[3]*data$x1*data$x2+theta[4]*data$x2^2-data$Y*data$x2
        )
    }
}

ef1 <- function(data){
    function(theta){
        ## ols z=1
        xb1 <- with(data,(cbind(1,S,x1,x2)%*%theta[1:4])[,1])

        ## logit z=1
        xb3 <- with(data,(cbind(1,x1,x2)%*%theta[5:7])[,1])
        ps <- plogis(xb3)

        out <- c(
        ### regression for treatment group
            data$Y-xb1,
            data$S*(data$Y-xb1),
            data$x1*(data$Y-xb1),
            data$x2*(data$Y-xb1),
### logistic regression
            data$S-ps,
            data$x1*(data$S-ps),
            data$x2*(data$S-ps)
        )

        out
    }
}


mmm1 <- m_estimate(ef1,dat,root_control = setup_root_control(start = rep(0.1,7)))



geexNoInt <- function(data){
    function(theta){
        xmat0 <- with(data,cbind(1,x1,x2))
        ps <- plogis(xmat0%*%theta[7:9])
        r0 <- ifelse(data$Z==1,data$S-ps,0)
        Sp <- ifelse(data$Z==1,data$S,ps)
        ZSp <- data$Z*Sp

        xmat <- cbind(1,data$Z,Sp,ZSp,data$x1,data$x2)
        xb <- xmat%*%theta[1:6]
        r <- data$Y-xb
        rbind(
            t(as.vector(r)*xmat),
            t(as.vector(r0)*xmat0)
        )
    }
}

