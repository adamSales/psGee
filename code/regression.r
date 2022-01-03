library(sandwich)

### I want a bunch of attributes for the VCV matrix, but I don't
### want them to print every time!
print.vcv <- function(x,...) print(x[1:NROW(x),1:NCOL(x)],...)

### convenience function for functions that return lists
### why cant R be more like python in this regard??
Attach <- function(.list,.names=names(.list),env=environment()){
    if(is.null(.names))
        .names <- paste0(
            as.character(
                as.list(match.call())$.list)[1],
            seq(length(.list)))

    for(i in 1:length(.list)){
        assign(.names[i],.list[[i]], envir=parent.frame())
    }

}

AUC = function(probs, true_Y){
    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    val = unlist(probsSort$x)
    idx = unlist(probsSort$ix)

    roc_y = true_Y[idx];
    stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)

    auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(auc)
}
AUCmod <- function(mod){
    AUC(mod$linear,mod$y)
}

datNames <- function(data,trt='Z',out='Y',use='S',block=NULL){
  data$Y <- as.numeric(data[[out]])
  data$Z <- as.numeric(data[[trt]])
  data$S <- as.numeric(data[[use]])
  
  data$S[data$Z==0] <- NA
  
  if(!is.null(block)){
    data$block <- as.factor(data[[block]])
    data <- data%>%
      group_by(block)%>%
      mutate(
        Zadj=Z-mean(Z),
        Yadj=Y-mean(Y)
      )
  } else{
    data <- mutate(data,
                   Z=Z,
                   Yadj=Y)
  }
  data
}

### computes/returns regression models
pointEst <- function(data,covFormU=~x1+x2,covFormY=covFormU,int=FALSE,psMod=NULL){
  
  data$S[data$Z==0] <- NA
  
    if(is.null(psMod)) psMod <- glm(
      update(covFormU,S~.),
      data=data,family=binomial,
      subset=!is.na(S))
    else covForm <- formula(psMod)[c(1,3)]

    attr(psMod,'auc') <- AUCmod(psMod)

    ps <- predict(psMod,data,type='response')
    if('Sp'%in%names(data)) warning('replacing Sp')
    data <- within(data,Sp <- ifelse(Z==1&!is.na(S),S,ps))

    outForm <- 
      update(covFormY,
        if(int) Y~Z*(Sp+.) else Y~Z*Sp+.
        )
    
    if(!is.null(data$block))
      outForm <- update(outForm,.~.+block)
    
    outMod <- lm(outForm,data=data)

    list(psMod=psMod,outMod=outMod)
}

### estimating equations for stage 1 model (PS)
EFps <- function(psMod,data){
    efPStmp <- estfun(psMod)
    efPS <- matrix(0,nrow(data),ncol(efPStmp))
    efPS[as.numeric(rownames(efPStmp)),] <- efPStmp
    efPS
}

### estimates bread and meat for sandwich variance
sandwichMats <- function(psMod,outMod,data,clust=NULL,int=any(grepl(":x",names(coef(outMod))))){
    #data <- model.frame(outMod)
    ### estimating equations
    efPS <- EFps(psMod,data)
    efOut <- estfun(outMod)

    out <- list(
        a11inv = bread(psMod),#/sum(data$Z),
        a22inv = bread(outMod),#/nrow(data),
        a21 = A21(psMod,outMod,data)/nrow(data),
        b11=(if(is.null(clust)) meat(psMod) else meatCL(psMod,cluster=clust[data$Z==1])),#*nrow(model.frame(psMod)),
        b22 = (if(is.null(clust)) meat(outMod) else meatCL(outMod,cluster=clust))#*nrow(data)#,adjust=TRUE)
    )
    out <- within(out,b12 <-
        if(int){
            matrix(0,nrow(b11),ncol(b22))
        } else if(is.null(clust)) crossprod(efPS,efOut)/nrow(data)
        else crossprod(apply(efPS,2L,rowsum,clust),
                       apply(efOut,2L,rowsum,clust))#/sum(data$Z==0)
        )
    out
}

### combining smaller matrices into bigger one
bigM <- function(m11,m12,m22)
    rbind(
        cbind(m11,m12),
        cbind(t(m12),m22)
    )

### computes sandwich vcov for regressions
vcvPS <- function(psMod,outMod,data,clust=NULL,int=any(grepl(":x",names(coef(outMod))))){
    Attach(
        sandwichMats(
            psMod=psMod,
            outMod=outMod,
            data=data,
            clust=clust,
            int=int))

    A <- rbind(
        cbind(solve(a11inv),matrix(0,nrow(a11inv),ncol(a22inv))),
        cbind(a21,solve(a22inv))
    )
    colnames(A) <- rownames(A)
    ##bigM(solve(a11inv),t(a21),solve(a22inv))

    B <- bigM(b11,b12,b22)

    n <- nrow(model.frame(outMod))

    vcvFull <- solve(A)%*%B%*%t(solve(A))/nrow(data)
    main <- a22inv%*%b22%*%t(a22inv)
    vcv <- vcvFull[-(1:nrow(b11)),-(1:nrow(b11))]
    dimnames(vcv) <- dimnames(main)
    correction <- vcv-main

    attr(vcv,"full") <- vcvFull
    attr(vcv,"bread") <- A
    attr(vcv,"meat") <- B
    attr(vcv,"main") <- main
    attr(vcv,"correction") <- correction

    class(vcv) <- "vcv"

    vcv
}



### wrapper function for estimating regressions+vcov
est <- function(data,covFormU=~x1+x2,covFormY=covFormU,psMod=NULL,clust=NULL,int=FALSE,
                trt='Z',out='Y',use='S',block=NULL){

  data <- datNames(data,trt=trt,out=out,use=use,block=block)
  
    Attach(pointEst(data=data,covFormU=covFormU,covFormY=covFormY,int=int,psMod=psMod))

    vcv <- vcvPS(psMod,outMod,data=data,int=int,clust=clust)

    list(outMod=outMod,psMod=psMod,vcv=vcv)
}

### estimates effects of interest, starting from est() output
effsFromFit <- function(ests){
    estimates <- with(as.list(coef(ests$outMod)),
                      list(
                          eff0=Z,
                          eff1=Z+`Z:Sp`,
                          diff=`Z:Sp`))
    vcv <- ests$vcv
    ddd <- diag(vcv)
    vars <- list(
        eff0 <- ddd['Z'],
        eff1 <- ddd['Z']+ddd['Z:Sp']+2*vcv['Z','Z:Sp'],
        diff <- ddd['Z:Sp']
    )
    out <-
        cbind(
            estimates=unlist(estimates),
            SE=sqrt(unlist(vars))
        )
    attr(out,'auc') <- attr(ests$psMod,'auc')
    out
}

### estimates effects of interest, starting from data
effs <- function(data,covFormU=~x1+x2,covFormY=covFormU,int=FALSE,
                 trt='Z',out='Y',use='S',block=NULL){
    ests <- est(data=data,covFormU=covFormU,covFormY=covFormY,
                int=int,trt=trt,out=out,use=use,block=block)

    effsFromFit(ests)
}

### estimates lower left of A (bread) matrix in A^{-1}BA^{-t}
### FIX FOR MISSING S IN TRT GROUP
A21 <- function(psMod,outMod,data){
    risp <- data$Z==0|is.na(data$S)#,1,0)
    Z <- data$Z[risp]
    Y0 <- data$Y[risp]
    aW <- predict(psMod,data,type='link')[risp]
    q <- family(psMod)$mu.eta(aW)      ## dp/dalpha=qX'
    u <- 2*family(psMod)$linkinv(aW)*q   ## dp^2/dalpha=uX'

    W0 <- model.matrix(update(formula(psMod),Y~.),data=data)[risp,]
    
    X0 <- model.matrix(outMod)[risp,
                               -c(which(names(coef(outMod))%in%c('Z','Sp')),
                                  grep('Z\\:|Sp\\:|\\:Z|\\:Sp',names(coef(outMod))))]

    est <- list(beta=coef(outMod)[colnames(X0)[-1]],
                betaZ=coef(outMod)['Z'],
                betaZr=coef(outMod)['Z:Sp'],
                eta=coef(outMod)['Sp'],
                mu0=coef(outMod)['(Intercept)']
                )

    a21=rbind(
       t(est$eta*q+est$betaZr*Z)%*%W0, # 1 x p1
       -t((Y0-est$mu0-est$betaZ*Z-X0[,-1]%*%est$beta)*q-est$eta*u)%*%W0,
       t(Z*(est$eta*q+est$betaZr))%*%W0, # 1 x p1
       -t(((Y0-est$mu0-est$betaZ*Z-X0[,-1]%*%est$beta)*q-est$eta*u)*Z)%*%W0,
       t((est$eta*q+est$betaZr*Z)*X0)%*%W0 # p2 x p1
    )#/nrow(X0)
    ## rownames(a21)[1:2] <- c('(Intercept)','Sp')

    ## t(vapply(
    ##     names(coef(outMod)),
    ##     function(n)
    ##         if(n %in% rownames(a21)) a21[n,] else rep(0,ncol(a21)),
    ##     numeric(ncol(a21)))
    ##   )
    a21
}

