bsSchool <- function(dat){
  datS <- split(dat,dat$SchIDPre)
  do.call("rbind",
          lapply(datS,function(x) x[sample(1:nrow(x),nrow(x),replace=TRUE),]))
}


#### PSW
psw1 <- function(dat,psMod){

  dat0=subset(dat,Z==0)
  dat1=subset(dat,Z==1)
  dat0$ps=predict(psMod,dat0,type='response')


  muc0=with(dat0,sum(Y*(1-ps))/sum(1-ps))
  muc1=with(dat0,sum(Y*ps)/sum(ps))

  mut0=with(dat1,mean(Y[S==0]))
  mut1=with(dat1,mean(Y[S==1]))

  c(eff0=mut0-muc0,
    eff1=mut1-muc1,
    diff=mut1-muc1-mut0+muc0
    )
}

print.psw <- function(x, ...) print(x$coef,...)


psw <- function(alt,psMod,B=5000,verbose=TRUE){
  dat <- getDat(alt)

  est=psw1(dat,psMod)

  if(verbose) step <- if(B>10) round(B/10) else 1
  bs <- matrix(nrow=B,ncol=3)
  for(i in 1:B){
    if(verbose) if(i%%step==0) cat(round(i/B*100),' ')
    datStar <- bsSchool(dat)
    psModStar <- update(psMod,data=subset(datStar,Z==1))
    bs[i,] <- psw1(datStar,psModStar)
  }
  if(verbose) cat('\n')
  out <-
    list(
      coef=cbind(est=est,se=apply(bs,2,sd)),
      bs=bs)
  class(out) <- 'psw'
  out
}
