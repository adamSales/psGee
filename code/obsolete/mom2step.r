twoStep <- function(data,justEst=TRUE){
    mod1 <- lm(Y~x1+x2+S,data=data,subset=Z==1)
    lmod <- glm(S~x1+x2,data=data,subset=Z==1,family=binomial)

    dat1 <- data.frame(
        r=data$Y[data$Z==0]-#coef(mod1)['x1']*data$x1[data$Z==0]-coef(mod1)['x2']*data$x2[data$Z==0],
                                        predict(mod1,subset(data,Z==0)),
        ps=predict(lmod,subset(data,Z==0),type='response')
    )

    estFun2 <- function(data){
        function(theta){
            c(
                theta[1]-data$ps,
                theta[2]-data$ps^2,
                theta[4]*theta[1]+theta[3]*(1-theta[1])-data$r,
                theta[4]*theta[2]+theta[3]*(theta[1]-theta[2])-data$r*data$ps
            )
        }
    }

    if(justEst){
        bbb <- create_basis(estFun2,dat1)
        mycontrol <- new('geex_control', .root = setup_root_control(start = rep(.5,4)))
        bbb@.control <- mycontrol
        theta <- estimate_GFUN_roots(bbb)$root
        theta <- c(theta,
                   -theta[3],#+mean(data$Y[data$Z==1&data$S==0]),
                   coef(mod1)['S']-theta[4]
                   )
        names(theta) <- c('Epi','Epi2','mu0','mu1','eff0','eff1')

        theta <- c(theta,attr(data,"trueEff"))

        return(theta)
    }

    m_estimate(estFun2,data=dat1,root_control=setup_root_control(start = rep(.5,4)))
}





## #clusterExport(cl, c('twoStep','makeDat','sim2Step'))

## res2000 <- pbreplicate(500,sim2Step(2000),cl=cl)
## res500 <- pbreplicate(500,sim2Step(500),cl=cl)
## resNoSep  <- pbreplicate(500,sim2Step(2000,mu00=0,mu01=0,mu10=0.5,mu11=0.5),cl=cl)
## res0  <- pbreplicate(1000,sim2Step(500,mu00=0,mu01=0,mu10=0,mu11=0),cl=cl)
## res0.1000  <- pbreplicate(1000,sim2Step(1000,mu00=0,mu01=0,mu10=0,mu11=0),cl=cl)

##                                         #rowMeans(res2000)
