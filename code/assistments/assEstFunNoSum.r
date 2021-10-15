estFun <- function(data){
    XX <- as.matrix(select(data,-Y,-S,-Z))

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
            ifelse(ZZ==1,YY-xb1,0),
            ifelse(ZZ==1,SS*(YY-xb1),0),
            as.vector(apply(XX,2,function(x)
                ifelse(ZZ==1,x*(YY-xb1),x*(YY-xb2-a01*ps-a00*(1-ps))))),
### logistic regression
            ifelse(ZZ==1,SS-ps,0),
            as.vector(apply(XX,2,function(x) ifelse(ZZ==1,x*(SS-ps),0))),
### mixture model in control group
            ifelse(ZZ==0,a01*ps+a00*(1-ps)-YY+xb2,0),
            ifelse(ZZ==0,a01*ps^2+a00*(ps-ps^2)-(YY-xb2)*ps,0)#,
#            length(ZZ)*(eff0-(a10-a00)),
#            length(ZZ)*(eff1-(a11+a10-a01)),
#            length(ZZ)effDiff-(eff1-eff0)
        )
    }
}
