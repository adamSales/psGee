mest <- function(data){
    form <- Y~state+grade+race+sex+frl+xirt+esl
    data <- droplevels(data)
    X <- model.matrix(form,data=data)[,-1]
    p <- ncol(X)

    X <- scale(X)
    newDat <- data.frame(cbind(Y=data$Y,S=ifelse(data$everCP,1,0),Z= data$treatment,scl=data$schoolid2,X))

    estFun <- function(data){

        XX <- as.matrix(select(data,-Y,-S,-Z,-scl))

        YY <- data$Y

        ZZ <- data$Z
        SS <- data$S

        p <- ncol(XX)
        function(theta){
            ## YY model (ZZ=1): YY=a10+a11SS+a2XX
            ##         (ZZ=0): YY={a00 or a01}+a2XX
            ## SS model (ZZ=1): logit(SS)=b0+b1XX
            ## effects: a10-a00 and a10+a11-a01

            a10 <- theta[1]
            a11 <- theta[2]
            a00 <- theta[3]
            a01 <- theta[4]
            eff0 <- theta[5] ## a10-a00
            eff1 <- theta[6] ## a10+a11-a01
            effDiff <- theta[7] ## a11+a00-a01
            a2 <- theta[8:(p+7)]
            b0 <- theta[p+8]
            b1 <- theta[(p+9):(2*p+8)]

            ## ols z=1
            ## print(length(SS))
            ## print(dim(XX))
            xb1 <- (cbind(1,SS,XX)%*%c(a10,a11,a2))[,1]

            ## xb for z=0
            xb2 <- (XX%*%a2)[,1] ## (a2 is same in trt groups)

            ## logit z=1
            xb3 <- (cbind(1,XX)%*%c(b0,b1))[,1]
            ps <- plogis(xb3)

            c(
### regression for treatment group
                sum(ifelse(ZZ==1,YY-xb1,0)),
                sum(ifelse(ZZ==1,SS*(YY-xb1),0)),
                as.vector(apply(XX,2,function(x)
                    sum(ifelse(ZZ==1,x*(YY-xb1),x*(YY-xb2-a01*ps-a00*(1-ps)))))),
### logistic regression
                sum(ifelse(ZZ==1,SS-ps,0)),
                as.vector(apply(XX,2,function(x) sum(ifelse(ZZ==1,x*(SS-ps),0)))),
### mixture model in control group
                sum(ifelse(ZZ==0,a01*ps+a00*(1-ps)-YY+xb2,0)),
                sum(ifelse(ZZ==0,a01*ps^2+a00*(ps-ps^2)-(YY-xb2)*ps,0)),
                eff0-(a10-a00),
                eff1-(a11+a10-a01),
                effDiff-(eff1-eff0)
            )
        }
    }

    t1 <- Sys.time()
    res <- try(m_estimate(estFun,newDat,units="scl",
                          root_control = setup_root_control(start = rep(0.1,2*p+8))))
    if(inherits(res,'try-error'))
        res <- try(m_estimate(estFun,newDat,units="scl",root_control = setup_root_control(start = rep(1,2*p+8))))

    t2 <- Sys.time()

    attr(res,"time") <- t2-t1

    res
}

ests <- function(res,data){

    form <- Y~state+grade+race+sex+frl+xirt+esl
    data <- droplevels(data)
    XX <- model.matrix(form,data=data)[,-1]


    est <- coef(res)
    se <- sqrt(diag(vcov(res)))


    out <- cbind(est,se)
    rownames(out) <- c('a10','a11','a00','a01','eff0','eff1','effDiff',paste0('YYreg:',colnames(XX)),'b0',paste0('Sreg:',colnames(XX)))
    colnames(out) <- c('est','se')

    out
}
