library(geex)
library(pbapply)
### try it out

### under simple assumption
n <- 1000

mu0 <- 0
mu1 <- 0.2

ps <- runif(n)

S <- rbinom(n,1,ps)

Y <- ifelse(S==1,mu1,mu0)+rnorm(n)

data <- data.frame(ps,Y)

### parameter vector
## mu0, mu1

estFun1 <- function(data){
    function(theta){
        c(theta[2]*data$ps+theta[1]*(1-data$ps)-data$Y,
          theta[2]*data$ps^2+theta[1]*(data$ps-data$ps^2)-data$Y*data$ps)
    }
}

bbb <- create_basis(estFun,data)
mycontrol <- new('geex_control', .root = setup_root_control(start = c(1, 1)))
bbb@.control <- mycontrol

estimate_GFUN_roots(bbb)

mu0hat <- mean(ps)*mean(Y*ps)/((mean(Y)-1)*mean(ps^2)+mean(ps)^2)

results <- m_estimate(estFun,data=data,    root_control = setup_root_control(start = c(1,1)))


sim <- function(){
    mu0 <- 0
    mu1 <- 0.2

    ps <- runif(n)

    S <- rbinom(n,1,ps)

    Y <- ifelse(S==1,mu1,mu0)+rnorm(n)

    data <- data.frame(ps,Y)

    results <- m_estimate(estFun,data=data,    root_control = setup_root_control(start = c(1,1)))

    c(coef(results), sqrt(diag(vcov(results))))
}


res <- replicate(1000,sim())


#### bigger

n <- 200

x1 <- rnorm(2*n)
x2 <- rnorm(2*n)
x3 <- rnorm(2*n)

psTrue <- plogis(x1+x2+x3)

S <- rbinom(2*n,1,psTrue)

Z <- rep(c(1,0),n)

mu10 <- 0
mu00 <- 0
mu11 <- .5
mu01 <- 0.2

Y <- 0.5*(x1+x2+x3)+ifelse(Z==1,ifelse(S==1,mu11,mu10),ifelse(S==1,mu01,mu00))+rnorm(2*n,0,0.01)

data <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2)


estFun <- function(data){

    function(theta){
        xb1 <- with(data,(cbind(1,S,x1,x2)%*%theta[1:4])[,1]) ## ols z=1
        xb2 <- with(data,(cbind(x1,x2)%*%theta[12:13])[,1]) ## just xb for z=0 (b is same in trt groups)
        xb3 <- with(data,(cbind(1,x1,x2)%*%theta[5:7])[,1]) ## logit z=1

        ps <- plogis(xb3)

        c(
### regression for treatment group
            ifelse(data$Z==1,data$Y-xb1,0),
            ifelse(data$Z==1,data$S*(data$Y-xb1),0),
            ifelse(data$Z==1,data$x1*(data$Y-xb1),0),
            ifelse(data$Z==1,data$x2*(data$Y-xb1),0),
### logistic regression
            ifelse(data$Z==1,data$S-ps,0),
            ifelse(data$Z==1,data$x1*(data$S-ps),0),
            ifelse(data$Z==1,data$x2*(data$S-ps),0),
### mixture model in control group
            ifelse(data$Z==0,theta[9]*ps+theta[8]*(1-ps)-data$Y+xb2,0),
            ifelse(data$Z==0,data$x1*(theta[9]*ps+theta[8]*(1-ps)-data$Y+xb2),0),
            ifelse(data$Z==0,data$x2*(theta[9]*ps+theta[8]*(1-ps)-data$Y+xb2),0),
            ifelse(data$Z==0,theta[9]*ps^2+theta[8]*(ps-ps^2)-(data$Y-xb2)*ps,0),
            theta[10]-(theta[1]-theta[8]),
            theta[11]-(theta[1]+theta[2]-theta[9])
        )
    }
}

system.time(res <- m_estimate(estFun,data,root_control = setup_root_control(start = rep(1,13))))

cbind(coef(res),sqrt(diag(vcov(res))))

bbb <- create_basis(estFun,data)
mycontrol <- new('geex_control', .root = setup_root_control(start = rep(1,13)))
bbb@.control <- mycontrol

system.time(ccc <- estimate_GFUN_roots(bbb))


sim2 <- function(debug=FALSE){
    n <- 200

    x1 <- rnorm(2*n)
    x2 <- rnorm(2*n)
    x3 <- rnorm(2*n)

    psTrue <- plogis(x1+x2+x3)

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)

    mu10 <- 0
    mu00 <- 0
    mu11 <- .5
    mu01 <- 0.2

    Y <- 0.5*(x1+x2+x3)+ifelse(Z==1,ifelse(S==1,mu11,mu10),ifelse(S==1,mu01,mu00))+rnorm(2*n,0,0.01)

    data <- data.frame(Y,Z,S=ifelse(Z==1,S,0),x1,x2)

    bbb <- create_basis(estFun,data)
    bbb@.control <- mycontrol

    out <-     estimate_GFUN_roots(bbb)$root

    if(debug)
        return(list(out,data))
    out

}


library(parallel)
cl <- makeCluster(6)

clusterEvalQ(cl,library(geex))

clusterExport(cl, c('estFun','sim2','mycontrol'))

system.time(rrr <- parLapply(cl,1:100,function(i) sim2()))

rrr <- do.call('rbind',rrr)

