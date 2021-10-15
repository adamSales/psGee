

sim1 <- function(n){

    x1 <- rnorm(2*n)
    x2 <- rnorm(2*n)
    x3 <- rnorm(2*n)

    x1 <- x1-mean(x1)
    x2 <- x2-mean(x2)
    x3 <- x3-mean(x3)

    psTrue <- plogis(x1+x2+x3)

    S <- rbinom(2*n,1,psTrue)

    Z <- rep(c(1,0),n)



    psDat <- data.frame(x1,x2,S,Z)
    psMod <- glm(S~x1+x2,data=psDat,subset=Z==1,family=binomial)
    ps <- predict(psMod,subset(psDat,Z==0),type='response')

    x1 <- x1[Z==0]
    x2 <- x2[Z==0]

    Y <- 0.5*(x1+x2+x3[Z==0])+.2*S[Z==0]+rnorm(n,0,.2)

    a <- rbind(
        c(sum(ps), sum(1-ps), sum(x1), sum(x2)),
        c(sum(ps^2), sum(ps*(1-ps)), sum(ps*x1), sum(ps*x2)),
        c(sum(x1*ps), sum(x1*(1-ps)), sum(x1*x1), sum(x1*x2)),
        c(sum(x2*ps), sum(x2*(1-ps)), sum(x2*x1), sum(x2*x2)))

    b <- c(sum(Y),sum(ps*Y),sum(x1*Y),sum(x2*Y))

    solve(a,b)
}
