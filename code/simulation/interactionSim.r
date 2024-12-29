sss0=function(n){
  x=rnorm(n)
  p=plogis(.5*x)
  S=rbinom(n,1,p)
  Y=0.3+0.3*S+rnorm(n)

  vy=var(Y)

  y0=mean(Y[S==0])

  beta=cov(p,Y)/var(p)
  p0=mean(p)
  c(mu0=y0,mu0hat=mean(Y)-beta*p0,vy=vy)
}

sss=function(n){
  x=rnorm(n)
  ps=plogis(.5*x)
  S=rbinom(n,1,ps)
  Y=.3-0.5*x+0.3*S+rnorm(n)
  vy=var(Y)
#  mod=lm(Y~x+ps)
  y0=mean(Y[S==0])
  x0=mean(x[S==0])
  rm(S);gc()
                                        #  c(coef(mod),coef(mod)[1]+coef(mod)[2]*mean(x[S==0]),mean(Y[S==0]))
  vx=var(x)
  xy=cov(x,Y)
  xp=cov(x,ps)
  xb=mean(x)
  rm(x); gc()
  py=cov(ps,Y)
  Yb=mean(Y)
  rm(Y); gc()
  vp=var(ps)
  pb=mean(ps)

  m1=(py-xp*xy/vx)/(vp-xp^2/vx)
  b=(xy*vp-py*xp)/(vx*vp-xp^2)
  m0=Yb-b*xb-m1*pb


  c(mu0=y0,mu0hat=m0+x0*b,m1,b,m0,vy=vy)

}


sssInt=function(n){
  x=rnorm(n)
  p=plogis(.5*x)
  S=rbinom(n,1,p)
  y=.3-0.5*x+0.3*S+.3*S*x+rnorm(n)
  vy=var(y)
#  mod=lm(y~x+ps)
  y0=mean(y[S==0])
  x0=mean(x[S==0])
  rm(S);gc()
                                        #  c(coef(mod),coef(mod)[1]+coef(mod)[2]*mean(x[S==0]),mean(y[S==0]))
  px=mean(p*x)
  p2x=mean(p^2*x)
  px2=mean(p*x^2)
  p2x2=mean(p^2*x^2)
  x2=mean(x^2)
  xy=mean(x*y)
  pxy=mean(p*x*y)
  x=mean(x)
  gc()

  p2=mean(p^2)
  py=mean(p*y)
  p=mean(p)
  gc()

  y=mean(y)
  gc()

  beta=solve(
    rbind(c(1,x,p,px),
          c(p,px,p2,p2x),
          c(x,x2,px,px2),
          c(px,px2,p2x,p2x2)),
    c(y,py,xy,pxy))

  c(mu0=y0,mu0hat=beta[1]+beta[2]*x0,beta,vy=vy)
}

intSim=function(nrep,simFun){

  lapply(2:6, function(ex){cat(ex,' '); replicate(nrep,simFun(as.integer(10^ex)))})
}


summarizeSim=function(res){
  names(res)<- paste0('n=1e',2:6)
  err=sapply(res,function(x) x['mu0hat',]-x['mu0',],simplify=FALSE)

  bias=sapply(err,function(x) unlist(t.test(x)[c('estimate','conf.int')]))

  vv=sapply(err,var)

  vvRatio=vv*10^c(2:6)/vapply(res,function(x) mean(x['vy',]),1.0)

  list(err=err,bias=bias,vv=vv,vvRatio=vvRatio)
}

plotSim=function(res,title=''){
  summs=summarizeSim(res)

  boxplot(summs$err,ylab='Estimation Error',main=title)

  if(require('ggplot2')&require('dplyr')){
    biasFig=t(summs$bias)%>%
      as.data.frame()%>%
      `names<-`(c('est','ciL','ciH'))%>%
      rownames_to_column('n')%>%
      ggplot(aes(n,est,ymin=ciL,ymax=ciH))+geom_point()+geom_errorbar(width=0)+ggtitle(title,subtitle='Bias')

    print(biasFig)
  }

  nrep=ncol(res[[1]])
  q1=qchisq(.025,nrep-1)
  q2=qchisq(.975,nrep-1)
  Hmisc::errbar(2:6,summs$vvRatio,summs$vvRatio*(nrep-1)/q1,summs$vvRatio*(nrep-1)/q2,
                ylab='n*Var(muc0hat)/Var(Yc)',
                xlab='log_10(n)',main=title)
}
