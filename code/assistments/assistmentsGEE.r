


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
            sum(ZZ*(YY-xb1)),
            sum(ZZ*SS*(YY-xb1)),
            (t(XX)%*%(YY-(ZZ*xb1+(1-ZZ)*(xb2+a01*ps+a00*(1-ps)))))[,1],
### logistic regression
            sum(ZZ*(SS-ps)),
            ZZ*(t(XX)%*%(SS-ps))[,1],
### mixture model in control group
            sum((1-ZZ)*(a01*ps+a00*(1-ps)-YY+xb2)),
            sum((1-ZZ)*(a01*ps^2+a00*(ps-ps^2)-(YY-xb2)*ps)),
            length(ZZ)*(eff0-(a10-a00)),
            length(ZZ)*(eff1-(a11+a10-a01)),
            length(ZZ)*(effDiff-(eff1-eff0))
        )
    }
}

system.time(
    result1 <- m_estimate(estFun,newDat,root_control=setup_root_control(start=runif(2*ncol(X)+8,-.2,.2)))
)

est1 <- coef(result1)
vcv1 <- vcov(result1)
sand1 <- result1@sandwich_components

save(ate1,ate2,ate3,est1,vcv1,sand1,estFun,newDat,file='assistmentsResult.RData')
try(system('drive push --quiet assistmentsResult.RData'))

mod1 <- lm(Y~.-Z,data=newDat,subset=Z==1)
mod2 <- glm(S~.-Z-Y,data=newDat,subset=Z==1,family=binomial)

system.time(
    result2 <- m_estimate(estFun,newDat,root_control = setup_root_control(
                                           start = c(coef(mod1)[1],
                                                     coef(mod1)['S'],
                                                     rep(.1,5), #a0*, effects, effDiff
                                                     coef(mod1)[-c(1,which(names(coef(mod1))=='S'))],
                                                     coef(mod2))
                                       )
                         )
)

est2 <- coef(result2)
vcv2 <- vcov(result2)
sand2 <- result2@sandwich_components
save(ate1,ate2,ate3,est1,vcv1,sand1,est2,vcv2,sand2,estFun,newDat,file='assistmentsResult.RData')

try(system('drive push --quiet assistmentsResult.RData'))
