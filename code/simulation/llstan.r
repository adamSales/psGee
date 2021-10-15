attach(sdat)

a00 <- runif(1,-2,2)
a01 <- runif(1,-2,2)
#a10 <-  runif(1,-2,2)
#a11 <- runif(1,-2,2)
#a2 <- runif(p,-2,2)
mm <- lm(Ytrt~St+xt)
a10 <- coef(mm)[1]
a11 <- coef(mm)[2]
a2 <- coef(mm)[-c(1,2)]

b0 <- runif(1,-2,2)

b1 <- runif(p,-2,2)
sclEffY <- runif(nscl,-2,2)

sclEffS <- runif(nscl,-2,2)
sigSclY <- runif(1,0,2)
sigSclS <- runif(1,0,2)
sigt <- runif(1,0,2)
sigc <- runif(1,0,2)

xb1 <- a01+a11*St+xt%*%a2+sclEffY[sclt]
xb2 <- xc%*%a2+sclEffY[sclc]
xb3 <- b0 +xt%*%b1+sclEffS[sclt]
xb4 <- b0 +xc%*%b1+sclEffS[sclc]

probt <- plogis(xb3)
probc <- plogis(xb4)


ll <- c(aprior=sum(dnorm(c(a01,a10,a01,a11),log=TRUE)),
        sclEffY=sum(dnorm(sclEffY,0,sigSclY,log=TRUE)),
        sclEffS=sum(dnorm(sclEffS,0,sigSclS,log=TRUE)),
        St=sum(dbinom(St,1,probt,log=TRUE)),
        Ytrt=sum(dnorm(Ytrt,xb1,sigt,log=TRUE))
        )

llc <- 0
for(i in 1:nctl)
    llc <- llc+
        log(exp(log(probc[i])+dnorm(Yctl[i],a01+xb2[i],sigc,log=TRUE))+
            exp(log(1-probc[i])+dnorm(Yctl[i],a00+xb2[i],sigc,log=TRUE)))




