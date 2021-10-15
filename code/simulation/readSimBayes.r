library(tidyverse)
library(rstan)
library(geex)
library(parallel)
library(pbapply)

oneRun <- function(res1,facs){
    out <- c(
        facs,
        setNames(res1$true,c('true0','true1')),
        res1$mest,
        setNames(res1$bayes[1,c('mean','sd','Rhat')],c('mean0','sd0','Rhat0')),
        setNames(res1$bayes[2,c('mean','sd','Rhat')],c('mean1','sd1','Rhat1')),
        diffTrue=res1$true['S1']-res1$true['S0'],
        diffMest=res1$mest['eff1']-res1$mest['eff0'],
        diffBayes=res1$bayes['effDiff','mean']
    )

    as_tibble(rbind(out))
}

loadRes <- function(){
    files <- list.files('./simResults','Bayes.RData')
    results <- list()
    for(ff in files[startsWith(files,'sim')]){
        load(paste0('simResults/',ff))
        results[[ff]] <- map_dfr(res,oneRun,facs=facs)%>%map_dfc(unlist)
    }
    results
}

bias <- function(res,s){
    out=c(
        mest=mean(res[[paste0('true',s)]]-res[[paste0('eff',s)]]),
        mest.se=sd(res[[paste0('true',s)]]-res[[paste0('eff',s)]])/sqrt(nrow(res)),
        bayes=mean(res[[paste0('true',s)]]-res[[paste0('mean',s)]]),
        bayes.se=sd(res[[paste0('true',s)]]-res[[paste0('mean',s)]])/sqrt(nrow(res))
    )
    names(out) <- paste0(names(out),s)
    out
}

biasDiff <- function(res){
    out <- sweep(res[,c('diffMest.eff1','diffBayes')],2,res$diffTrue.S1,'-')
    c(colMeans(out),apply(out,2,sd)/sqrt(nrow(out)))
}

getBiases <- function(res){


}

summs <- function(res){
    facs <- res[1,map_dbl(results,var)==0]


bias <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    colMeans(cbind(sss1[,3:6]-sss1[,c(1,2,1,2)],sss2))
}

rmse <- function(sss){

    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]

    c(
        sqrt(colMeans((sss1[,3:6]-sss1[,c(1,2,1,2)])^2)),
        colMeans(sss2)
    )
}


summ <- function(sss){
    sss1 <- sss[,1:6]
    sss2 <- sss[,7:11]
    rbind(colMeans(sss1),
          apply(sss1,2,sd),
          sqrt(row1Means((sss-sss[rep(1:2,3),])^2))
          )
    }


plotRes <- function(rrr,fac,S){
    www <- (max(rrr[fac,])-min(rrr[fac,]))*0.02
    plot(rrr[fac,]-www,rrr[paste0('mom',S),],ylim=range(c(rrr[paste0('mom',S),],rrr[paste0('mle',S),])),xlim=c(min(rrr[fac,])-2*www,max(rrr[fac,])+2*www),pch=16,xlab=fac)
    points(rrr[fac,]+www,rrr[paste0('mle',S),],col='red',pch=16)
    legend('topright',legend=c('MOM','MLE'),col=c('black','red'),pch=16)
    abline(h=0,lty=2)
}


getSumms <- function(res,S,EST){
    #cases <- results%>%select(starts_with('mu'),n,b1,gumb)%>%distinct()
    true <- res[[paste0('true',S)]]
    est <- res[[paste0(EST,S)]]
    se <- res[[paste0(EST,S,'se')]]

    cbind(
        res[1,which(sapply(res,n_distinct)==1)],
        data.frame(
            EST=EST,
            S=S,
            getEst=mean(is.finite(est)),
            getSE=mean(is.finite(se)),
            bias=mean(est-true),
            rmse=sqrt(mean((est-true)^2)),
            biasV=mean(se^2)-var(est-true),
            avgSE=mean(se),
            reject=mean(est>=2*se)+mean(est<=2*se),
            cover=mean(true> est-2*se & true< est+2*se)
        )
    )
}




allSums <- function(results){
    res1 <- results[[1]]
    eee <- names(select(res1,ends_with('0se')))%>%gsub('0se','',.)

    summs <- list()
    for(EST in eee)
        for(S in c(0,1))
            summs[[paste(EST,S)]] <-
                map_dfr(results,getSumms,EST=EST,S=S)

    do.call('rbind',summs)
}

    ## I <- length(results)

##     load('simResults/casesreverse.RData')

##     for(i in 1:nrow(cases)){#length(fn)){
##         if(i %% 10==0) cat(round(i/nrow(cases)*100),'%')#length(fn)*100), '% ')
##         load(paste0('simResults/sim',i,'reverse.RData'))#fn[i]))
##         stopifnot(identical(facs,cases[i,]))
##         results[[i+I]] <- as.data.frame(res)
##         results[[i+I]]$n <- facs$n
##         rm(res,facs)
##     }

##     do.call('rbind',results)

## }

### Bias
err <- map(results,
            ~.x%>%mutate(errMom0=mom0-true0,errMom1=mom1-true1,
                         errMle0=mle0-true0,errMle1=mle1-true1)%>%
                select(starts_with('err'),mu01:n))

## bias b1 vs n
err%>%
    map_dfr(~with(.x,tibble(errMom=c(errMom0,errMom1),errMle=c(errMle0,errMle1),b1=rep(b1,2),n=rep(n,2),S=rep(c(0,1),each=nrow(.x)))))%>%
    pivot_longer(-c(b1,n,S),names_to='estimator',values_to="err",names_prefix='err')%>%
    group_by(n,b1,S,estimator)%>%
    summarize(bias=mean(err,na.rm=TRUE),ci=t.test(err)$conf.int)%>%
    mutate(hl=rep(c('l','h'),n()/2))%>%
    pivot_wider(names_from="hl",values_from="ci")%>%
    ungroup()%>%
    mutate(x=as.numeric(as.factor(n))+ifelse(estimator=='Mle',-.1,+.1),
           S=ifelse(S==0,'Prin. Strat. 0','Prin. Strat. 1'),
           b1=paste("beta=",b1),
           Estimator=ifelse(estimator=="Mom","GEE","MLE"))%>%
    ggplot(aes(x,bias,ymin=l,ymax=h,color=Estimator,group=Estimator))+
    geom_point()+
    geom_line()+
    geom_errorbar(width=0)+
    geom_hline(yintercept=0)+
    facet_grid(S~b1)+
    scale_x_continuous(name="Sample Size",breaks=1:3,labels=c(100,500,1000))
ggsave("biasBetaN.jpg",width=6.5,height=4)



## RMSE
err%>%
    map_dfr(~with(.x,tibble(errMom=c(errMom0,errMom1),errMle=c(errMle0,errMle1),b1=rep(b1,2),n=rep(n,2),S=rep(c(0,1),each=nrow(.x)))))%>%
    pivot_longer(-c(b1,n,S),names_to='estimator',values_to="err",names_prefix='err')%>%
    group_by(n,b1,S,estimator)%>%
    summarize(RMSE=sqrt(mean(err^2,na.rm=TRUE)))%>%
    ungroup()%>%
    mutate(x=as.numeric(as.factor(n))+ifelse(estimator=='Mle',-.1,+.1),
           S=ifelse(S==0,'Prin. Strat. 0','Prin. Strat. 1'),
           b1=paste("beta=",b1),
           Estimator=ifelse(estimator=="Mom","GEE","MLE"))%>%
    ggplot(aes(x,RMSE,color=Estimator,group=Estimator))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=0)+
    facet_grid(S~b1)+
    scale_x_continuous(name="Sample Size",breaks=1:3,labels=c(100,500,1000))
ggsave("rmseBetaN.jpg",width=6.5,height=4)

## coverage
cover <- map(results,
           ~.x%>%mutate(covMom0=abs(true0-mom0)/mom0se<2,covMom1=abs(true1-mom1)/mom1se<2,
                        covMle0=abs(true0-mle0)/mle0se<2,covMle1=abs(true1-mle1)/mle1se<2)%>%
                select(starts_with('cov'),mu01:n))



cover%>%
    map_dfr(~with(.x,tibble(covMom=c(covMom0,covMom1),covMle=c(covMle0,covMle1),b1=rep(b1,2),n=rep(n,2),S=rep(c(0,1),each=nrow(.x)))))%>%
    pivot_longer(-c(b1,n,S),names_to='estimator',values_to="cov",names_prefix='cov')%>%
    group_by(n,b1,S,estimator)%>%
    summarize(Coverage=mean(cov,na.rm=TRUE))%>%
    ungroup()%>%
    mutate(x=as.numeric(as.factor(n))+ifelse(estimator=='Mle',-.1,+.1),
           S=ifelse(S==0,'Prin. Strat. 0','Prin. Strat. 1'),
           b1=paste("beta=",b1),
           Estimator=ifelse(estimator=="Mom","GEE","MLE"))%>%
    ggplot(aes(x,Coverage,color=Estimator,group=Estimator))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=0.95,linetype="dotted")+
    facet_grid(S~b1)+
    scale_x_continuous(name="Sample Size",breaks=1:3,labels=c(100,500,1000))
ggsave("coverBetaN.jpg",width=6.5,height=4)

#### now compare with the MUs:
## bias mu01 vs mu10
err%>%
    map_dfr(~with(.x,tibble(errMom=c(errMom0,errMom1),errMle=c(errMle0,errMle1),b1=rep(b1,2),n=rep(n,2),mu01=rep(mu01,2),mu10=rep(mu10,2),S=rep(c(0,1),each=nrow(.x)))))%>%
    filter(b1==0.2)%>%
    pivot_longer(starts_with('err'),names_to='estimator',values_to="err",names_prefix='err')%>%
    group_by(mu10,mu01,S,n,estimator)%>%
    summarize(
        eff0=ifelse(mu10[1]==0,'Eff0 = 0','Eff0 = 0.3'),
        eff1=ifelse(mu01[1]==0.3,'Eff1 = 0','Eff1 = 0.3'),
        eff=paste(eff0,eff1,sep='\n'),
        bias=mean(err,na.rm=TRUE),ci=t.test(err)$conf.int)%>%
    mutate(hl=rep(c('l','h'),n()/2))%>%
    pivot_wider(names_from="hl",values_from="ci")%>%
    ungroup()%>%
    mutate(x=as.numeric(as.factor(n))+ifelse(estimator=='Mle',-.1,+.1),
           S=ifelse(S==0,'Prin. Strat. 0','Prin. Strat. 1'),
           Estimator=ifelse(estimator=="Mom","GEE","MLE"))%>%
    ggplot(aes(x,bias,ymin=l,ymax=h,color=Estimator,group=Estimator))+
    geom_point()+
    geom_line()+
    geom_errorbar(width=0)+
    geom_hline(yintercept=0)+
    facet_grid(S~eff)+
    scale_x_continuous(name="Sample Size",breaks=1:3,labels=c(100,500,1000))+
    ggtitle('beta=0.2')
ggsave("biasMuN.jpg",width=6.5,height=4)



## RMSE
err%>%
    map_dfr(~with(.x,tibble(errMom=c(errMom0,errMom1),errMle=c(errMle0,errMle1),b1=rep(b1,2),n=rep(n,2),mu01=rep(mu01,2),mu10=rep(mu10,2),S=rep(c(0,1),each=nrow(.x)))))%>%
    filter(b1==0.2)%>%
    pivot_longer(starts_with('err'),names_to='estimator',values_to="err",names_prefix='err')%>%
    group_by(mu10,mu01,S,n,estimator)%>%
    summarize(
        eff0=ifelse(mu10[1]==0,'Eff0 = 0','Eff0 = 0.3'),
        eff1=ifelse(mu01[1]==0.3,'Eff1 = 0','Eff1 = 0.3'),
        eff=paste(eff0,eff1,sep='\n'),
        RMSE=sqrt(mean(err^2,na.rm=TRUE)))%>%
    ungroup()%>%
    mutate(x=as.numeric(as.factor(n))+ifelse(estimator=='Mle',-.1,+.1),
           S=ifelse(S==0,'Prin. Strat. 0','Prin. Strat. 1'),
           Estimator=ifelse(estimator=="Mom","GEE","MLE"))%>%
    ggplot(aes(x,RMSE,color=Estimator,group=Estimator))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=0)+
    facet_grid(S~eff)+
    scale_x_continuous(name="Sample Size",breaks=1:3,labels=c(100,500,1000))+
    ggtitle('beta=0.2')
ggsave("rmseMuN.jpg",width=6.5,height=4)

## coverage
cover%>%
    map_dfr(~with(.x,tibble(covMom=c(covMom0,covMom1),covMle=c(covMle0,covMle1),b1=rep(b1,2),n=rep(n,2),mu01=rep(mu01,2),mu10=rep(mu10,2),S=rep(c(0,1),each=nrow(.x)))))%>%
    filter(b1==0.2)%>%
    pivot_longer(starts_with('cov'),names_to='estimator',values_to="Coverage",names_prefix='cov')%>%
    group_by(mu10,mu01,S,n,estimator)%>%
    summarize(
        eff0=ifelse(mu10[1]==0,'Eff0 = 0','Eff0 = 0.3'),
        eff1=ifelse(mu01[1]==0.3,'Eff1 = 0','Eff1 = 0.3'),
        eff=paste(eff0,eff1,sep='\n'),
        Coverage=mean(Coverage,na.rm=TRUE))%>%
    ungroup()%>%
    mutate(x=as.numeric(as.factor(n))+ifelse(estimator=='Mle',-.1,+.1),
           S=ifelse(S==0,'Prin. Strat. 0','Prin. Strat. 1'),
           Estimator=ifelse(estimator=="Mom","GEE","MLE"))%>%
    ggplot(aes(x,Coverage,color=Estimator,group=Estimator))+
    geom_point()+
    geom_line()+
    geom_hline(yintercept=0.95,linetype="dotted")+
    facet_grid(S~eff)+
    scale_x_continuous(name="Sample Size",breaks=1:3,labels=c(100,500,1000))+
    ggtitle('beta=0.2')
ggsave("coverMuN.jpg",width=6.5,height=4)




results%>%
    bind_rows()%>%
    filter(b1>0,mu11>0)%>%#head()#,abs(mom0-true0)<5,abs(mle0-true0)<5)%>%
    group_by(n,mu01,mu10,b1)%>%
    summarize(across(c(mom0,mle0),list(est=~mean(.-true0,na.rm=TRUE),ci=~t.test(.-true0)$conf.int)))%>%#,across(c(mom1,mle1),list(est=~mean(.-true1),ci=~t.test(.-true1)$conf.int)))
    ungroup()%>%
    mutate(ci=rep(c('L','H'),n()/2))%>%
    select(-starts_with("true"))%>%#head()#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_wider(names_from=ci,values_from=ends_with('_ci'),names_sep='')%>%
    pivot_longer(-c(b1,n,mu01,mu10),names_to=c('Estimator','.value'),names_sep='_')%>%
    mutate(b1=b1+log(n,10)/50,sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,est,linetype=Estimator,color=factor(n),shape=Estimator,group=paste(Estimator,n),ymin=ciL,ymax=ciH))+geom_point()+geom_line()+geom_hline(yintercept=0,linetype=2)+geom_errorbar(width=0)+facet_grid(sepC~sepT)+
    ggtitle('Bias','Normal')


results%>%
    filter(gumb==1,b1>0,mu11>0)%>%
    group_by(n,b1,mu01,mu10)%>%
    summarize(across(c(true0,true1),mean),across(c(mom0,mle0),list(est=~mean(.-true0),ci=~t.test(.-true0)$conf.int)))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
    ungroup()%>%
    mutate(ci=rep(c('L','H'),n()/2))%>%
    select(-starts_with("true"))%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_wider(names_from=ci,values_from=ends_with('_ci'),names_sep='')%>%
    pivot_longer(-c(b1,n,mu01,mu10),names_to=c('Estimator','.value'),names_sep='_')%>%
    mutate(b1=b1+log(n,10)/50,sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,est,linetype=Estimator,color=factor(n),shape=Estimator,group=paste(Estimator,n),ymin=ciL,ymax=ciH))+geom_point()+geom_line()+geom_hline(yintercept=0,linetype=2)+geom_errorbar(width=0)+facet_grid(sepC~sepT)+
    ggtitle('Bias','Gumbel')+
    ylim(-.1,.1)


results%>%
    filter(gumb==0,b1>0,mu11>0)%>%
    group_by(n,b1,mu01,mu10)%>%
    summarize(across(c(mom0,mle0),list(est=~sqrt(mean((.-true0)^2)),ci=~sqrt(t.test((.-true0)^2)$conf.int))),across(c(true0,true1),mean))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
    ungroup()%>%
    mutate(ci=rep(c('L','H'),n()/2))%>%
    select(-starts_with("true"))%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_wider(names_from=ci,values_from=ends_with('_ci'),names_sep='')%>%
    pivot_longer(-c(b1,n,mu01,mu10),names_to=c('Estimator','.value'),names_sep='_')%>%
    mutate(b1=b1+log(n,10)/50,sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,est,linetype=Estimator,color=factor(n),shape=Estimator,group=paste(Estimator,n)))+#,ymin=ciL,ymax=ciH))+
    geom_point()+
    geom_line()+#geom_hline(yintercept=0,linetype=2)+
                                        #    geom_errorbar(width=0)+
    facet_grid(sepC~sepT)+
    ggtitle('RMSE','Normal') +
    coord_cartesian(ylim=c(0,0.5))


results%>%
    filter(gumb==1,b1>0,mu11>0)%>%
    group_by(n,b1,mu01,mu10)%>%
    summarize(across(c(mom0,mle0),list(est=~sqrt(mean((.-true0)^2)),ci=~sqrt(t.test((.-true0)^2)$conf.int))),across(c(true0,true1),mean))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
    ungroup()%>%
    mutate(ci=rep(c('L','H'),n()/2))%>%
    select(-starts_with("true"))%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_wider(names_from=ci,values_from=ends_with('_ci'),names_sep='')%>%
    pivot_longer(-c(b1,n,mu01,mu10),names_to=c('Estimator','.value'),names_sep='_')%>%
    mutate(b1=b1+log(n,10)/50,sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,est,linetype=Estimator,color=factor(n),shape=Estimator,group=paste(Estimator,n)))+#,ymin=ciL,ymax=ciH))+
    geom_point()+
    geom_line()+#geom_hline(yintercept=0,linetype=2)+
                                        #    geom_errorbar(width=0)+
    facet_grid(sepC~sepT)+
    ggtitle('RMSE','Gumbel') +
    coord_cartesian(ylim=c(0,0.5))

                                        #     ylim(0,.25)



results%>%
    filter(gumb==0,b1>0,mu11==0,mu01==0.3,mu10==0.3)%>%
    group_by(n,b1,mu01,mu10)%>%
    summarize(across(c(true0,true1),mean),across(c(mom0,mle0),list(est=~mean(.-true0),ci=~t.test(.-true0)$conf.int)))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
    ungroup()%>%
    mutate(ci=rep(c('L','H'),n()/2))%>%
    select(-starts_with("true"))%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_wider(names_from=ci,values_from=ends_with('_ci'),names_sep='')%>%
    pivot_longer(-c(b1,n,mu01,mu10),names_to=c('Estimator','.value'),names_sep='_')%>%
    mutate(b1=b1+log(n,10)/50)%>%#,sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,est,linetype=Estimator,color=factor(n),shape=Estimator,group=paste(Estimator,n),ymin=ciL,ymax=ciH))+geom_point()+geom_line()+geom_hline(yintercept=0,linetype=2)+geom_errorbar(width=0)#+#facet_grid(sepC~sepT)+
#    ggtitle('Bias','Normal')

results%>%
    filter(gumb==0,n==500,b1>0,mu11>0)%>%
    mutate(b1=factor(b1),across(c(mom0,mle0),~.-true0))%>%
    select(-starts_with("true"),-mle1,-mom1)%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_longer(c(mom0,mle0),names_to='Estimator',values_to='Estimate')%>%
    mutate(sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,Estimate,fill=Estimator))+geom_boxplot()+geom_hline(yintercept=0,linetype=2)+
    coord_cartesian(ylim=c(-1,1))+
    facet_grid(sepC~sepT)+
    ggtitle('n=500','Normal')


results%>%
    filter(gumb==0,n==1000,b1>0,mu11>0)%>%
    mutate(b1=factor(b1),across(c(mom0,mle0),~.-true0))%>%
    select(-starts_with("true"),-mle1,-mom1)%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_longer(c(mom0,mle0),names_to='Estimator',values_to='Estimate')%>%
    mutate(sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,Estimate,fill=Estimator))+geom_boxplot()+geom_hline(yintercept=0,linetype=2)+
    coord_cartesian(ylim=c(-1,1))+
    facet_grid(sepC~sepT)+
    ggtitle('n=1000','Normal')

results%>%
    filter(gumb==1,n==1000,b1>0,mu11>0)%>%
    mutate(b1=factor(b1),across(c(mom0,mle0),~.-true0))%>%
    select(-starts_with("true"),-mle1,-mom1)%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
    pivot_longer(c(mom0,mle0),names_to='Estimator',values_to='Estimate')%>%
    mutate(sepC=ifelse(mu01>0,'sepC','no sepC'),sepT=ifelse(mu10==0,'sepT','no sepT'))%>%
    ggplot(aes(b1,Estimate,fill=Estimator))+geom_boxplot()+geom_hline(yintercept=0,linetype=2)+
    coord_cartesian(ylim=c(-.1,.1))+
    facet_grid(sepC~sepT)+
    ggtitle('n=1000','Normal')


