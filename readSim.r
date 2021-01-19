library(tidyverse)
library(rstan)
library(geex)
library(parallel)
library(pbapply)

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


loadRes <- function(){
    load('simResults/cases.RData')

    results <- list()
    #fn <- list.files('./simResults','sim[[:alnum:]]+.RData')
    for(i in 1:nrow(cases)){#length(fn)){
        if(i %% 10==0) cat(round(i/nrow(cases)*100),'%')#length(fn)*100), '% ')
        load(paste0('simResults/sim',i,'.RData'))#fn[i]))
        stopifnot(identical(facs,cases[i,]))
        results[[i]] <- as.data.frame(res)
        results[[i]]$n <- facs$n
        rm(res,facs)
    }

    I <- length(results)

    load('simResults/casesreverse.RData')

    for(i in 1:nrow(cases)){#length(fn)){
        if(i %% 10==0) cat(round(i/nrow(cases)*100),'%')#length(fn)*100), '% ')
        load(paste0('simResults/sim',i,'reverse.RData'))#fn[i]))
        stopifnot(identical(facs,cases[i,]))
        results[[i+I]] <- as.data.frame(res)
        results[[i+I]]$n <- facs$n
        rm(res,facs)
    }

    do.call('rbind',results)

}


results%>%
    filter(gumb==0,b1>0,mu11>0,abs(mom0-true0)<5,abs(mle0-true0)<5)%>%
    group_by(n,b1,mu01,mu10)%>%
    summarize(across(c(true0,true1),mean),across(c(mom0,mle0),list(est=~mean(.-true0),ci=~t.test(.-true0)$conf.int)))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
    ungroup()%>%
    mutate(ci=rep(c('L','H'),n()/2))%>%
    select(-starts_with("true"))%>%#b1,n,starts_with('mom0'),starts_with('mle0'),ci)%>%
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
    summarize(across(c(true0,true1),mean),across(c(mom0,mle0),list(est=~sqrt(mean((.-true0)^2)),ci=~sqrt(t.test((.-true0)^2)$conf.int))))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
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
    summarize(across(c(true0,true1),mean),across(c(mom0,mle0),list(est=~sqrt(mean((.-true0)^2)),ci=~sqrt(t.test((.-true0)^2)$conf.int))))%>%#,across(c(mom1,mle1),~mean(.-true1)))#%>%
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


