library(tidyverse)
library(broom)
library(kableExtra)

source('code/simulation/readSimFuncs.r')
#### after simulation results have been loaded and pre-processed...


print(load('simResults/fullResults.RData'))
print(load('simResults/pswResults.RData'))

if(!is.data.frame(results)) results <- reduce(results,bind_rows)

facsDF<-distinct(results,n,mu01,mu10,mu11,b1,errDist,intS,intZ,run)
stopifnot(all.equal(facsDF$run,1:length(pswResults)))
### transform pswResults
pswResults <-
  map(seq(length(pswResults)),
      function(i){
        res<-pswResults[[i]]
        facs<-attr(res,'facs')
        if(!all(as.data.frame(facsDF[i,names(facs)])==facs)){
          print(i)
          stop()
        }
        res<-cbind(res,facs)
        res$run<-i
        res$effDiff=res$eff1-res$eff0
        res=pivot_longer(res,starts_with('eff'),names_to='eff',names_prefix='eff',values_to='est')
        samp<-unique(na.omit(results$samp[results$run==i]))
        res$samp <- if(length(samp)==nrow(res)) samp else NA
        res
      })%>%
  reduce(bind_rows)

pswResults$estimator<-'psw'

### rhats
rhats <- results%>%
  filter(eff=='1',estimator=='bayes')%>%
  group_by(n,mu01,errDist,b1,intS,intZ)%>%
    summarize(Rhat1.1=mean(rhat<1.1,na.rm=TRUE),Rhat1.01=mean(rhat<1.01,na.rm=TRUE))%>%arrange(Rhat1.1)

hist(rhats$Rhat1.1)
hist(rhats$Rhat1.01)

lm(Rhat1.1~I(n/500)+mu01+errDist+b1+intS+intZ,data=rhats)

##### check some stuff out wrt

### how does AUC very w b1 and n?
results%>%ggplot(aes(as.factor(b1),auc))+geom_jitter(alpha=0.2)+geom_violin(draw_quantiles = 0.5)

results%>%ggplot(aes(as.factor(n),auc))+geom_jitter(alpha=0.2)+geom_violin(draw_quantiles = 0.5)+facet_wrap(~b1,nrow=1)

#### how do estimates vary with auc (simple case)?
results%>%
  filter(n==500,errDist=='norm',!intS,!intZ,eff=='1')%>%
  group_by(b1)%>%
  mutate(AUC=round(mean(auc,na.rm=TRUE),1))%>%
  ungroup()%>%
  ggplot(aes(auc,est))+
  geom_point(aes(color=factor(b1)),alpha=0.2)+
  geom_boxplot(aes(AUC,est,group=AUC),outlier.shape = NA)+
  geom_smooth()+
  geom_hline(aes(yintercept=mu11-mu01))+
  facet_grid(estimator~mu01)+ylim(-1,1)



results<-bind_rows(results,pswResults)%>%
  group_by(run,eff)%>%
  mutate(pop=na.omit(pop)[1],samp=ifelse(is.na(samp),pop,samp))%>%
  ungroup()


pd <- results %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=c(0,0.2,0.5),
                labels=c(bquote(alpha==0),bquote(alpha==0.2),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='M-Est',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction")
  ) %>%
  filter(estimator=='PSW'|rhat<1.1)

### boxplots for when M-Estimation assumptions hold+normal erros




bp(pd,
   subset=b1 > 0& errDist=='norm'& !intS& !intZ&eff!='Diff'&n==500&PE=='Stratum 1',
   title="Normal Residuals, No Interactions",
   facet = B1~ M1,
   ylim=c(-1,1),
   Labeller=label_parsed
   )+
  labs(subtitle=bquote("Stratum 1 Principal Effects; n=500,"~alpha>0))
ggsave('simFigs/normalOutcomesNoInteractionPosB1.jpg',height=4,width=6)

bp(pd,
   subset=b1 > 0& errDist=='unif'& !intS& !intZ&eff!='Diff'&n==500&PE=='Stratum 1',
   title="Uniform Residuals, No Interactions",
   facet = B1~ M1,
   ylim=c(-1,1),
   Labeller=label_parsed
   )+
  labs(subtitle=bquote("Stratum 1 Principal Effects; n=500,"~alpha>0))
ggsave('simFigs/unifOutcomesNoInteractionPosB1.jpg',height=4,width=6)

bp(pd,
   subset=b1 > 0& errDist=='norm'& n==500&PE=='Stratum 1'&b1==.5&mu01==0.3,
   title="Interactions with S and Z, Normal Residuals",
   facet = interactionZ~interactionS,
   ylim=NULL
   )+
  labs(subtitle=bquote("Stratum 1 Principal Effects; "~alpha==0.5~","~mu[C]^1-mu[C]^0==0.3~","~ n==500))
ggsave('simFigs/InteractionPosB1.jpg',height=4,width=6)



challengeAssumptions <- list(
### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& !intS& !intZ&errDist!='norm',
   title="Stratum 1, No Interaction, b1>0",
   facet=AUCff+N~errDist+ M1),


### now the same, but with interaction in S

bp(pd,
   subset=b1 > 0& errDist=='norm'& intS& !intZ&eff!='Diff',
   title="Normal Outcomes, Interaction in S, b1>0",
   facet = AUCff+N~PE + M1
) ,

### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& intS& !intZ&errDist!='norm',
   title="Stratum 1, Interaction in S, b1>0",
   facet=AUCff+N~errDist+ M1),



### now the same, but with interaction in Z

bp(pd,
   subset=b1 > 0& errDist=='norm'& !intS& intZ&eff!='Diff',
   title="Normal Outcomes, Interaction in Z, b1>0",
   facet = AUCff+N~PE + M1
) ,

### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& !intS& intZ&errDist!='norm',
   title="Stratum 1, Interaction in Z, b1>0",
   facet=AUCff+N~errDist+ M1),

### both interactions
### now the same, but with interaction in S

bp(pd,
   subset=b1 > 0& errDist=='norm'& intS& intZ&eff!='Diff',
   title="Normal Outcomes, Both Interactions , b1>0",
   facet = AUCff+N~PE + M1
) ,

### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& intS& intZ&errDist!='norm',
   title="Stratum 1, Both Interaction, b1>0",
   facet=AUCff+N~errDist+ M1)
)

#challengeAssumptions <- do.call(marrangeGrob, c(challengeAssumptions, list(nrow = 2, ncol = 2)))

pdf("simFigs/challengeAssumptionsBoxplots.pdf", width = 6.5, height = 9,onefile = TRUE)

lapply(challengeAssumptions,print)
dev.off()

#######################################################
### summary stats
#######################################################

#######################################################
#### bias
#######################################################
bias <- results%>%
  filter(estimator=='psw'|rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,intS,intZ,eff,estimator)%>%
  summarize(bias=mean(est-pop,na.rm=TRUE),
            ttest=tidy(t.test(est-pop)))%>%
  ungroup()%>%
  bind_cols(.$ttest)%>%
  select(-ttest)

pdf('simFigs/biasFigs.pdf',height = 9,width=6.5)

### normal errors, no interaction
bias%>%
  filter(errDist=='norm',!intS,!intZ)%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')%>%print()

# with error bars:
bias%>%
  filter(errDist=='norm',!intS,!intZ,b1>0)%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width=0)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')%>%print()

### non-normal errors, no interaction
bias%>%
  filter(errDist!='norm',eff=='1',!intS,!intZ)%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(errDist~mu01+n,scales='free')+
  ggtitle('non-normal errors, no interaction stratum 1')%>%print()


#### interactions, normal errors
bias%>%
  filter(errDist=='norm',eff=='1')%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'normal errors, stratum 1, interactions')%>%print()


#### interactions, mixture errors
bias%>%
  filter(errDist=='mix',eff=='1')%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'mixture errors, stratum 1, interactions')%>%print()


#### interactions, uniform errors
bias%>%
  filter(errDist=='unif',eff=='1')%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'uniform errors, stratum 1, interactions')%>%print()

dev.off()


#######################################################
#### rmse
#######################################################
 rmse <- pd%>%
  filter(rhat<1.1|estimator=='PSW')%>%
  group_by(n,mu01,errDist,b1,intS,intZ,interactionS,interactionZ,eff,estimator)%>%
   summarize(rmse=sqrt(mean(errP^2,na.rm=TRUE)))%>%ungroup()

rmse%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1>0)%>%
  ggplot(aes(b1,rmse,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+#geom_hline(yintercept=0.95)+
  facet_grid(interactionS~errDist+interactionZ)

sink('writeUps/rmseTab.tex')
cbind(
rmse%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1==0)%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `X:Z\nInt.?`=ifelse(intZ,'Yes','No'),
         `X:S\nInt.?`=ifelse(intS,'Yes','No'),
         estimator=ifelse(estimator=='Mixture','Mix.',estimator),
         rmse)%>%
pivot_wider(names_from=estimator,values_from=rmse),
rmse%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1==.2)%>%
  pivot_wider(names_from=estimator,values_from=rmse)%>%
select(`M-Est`,`Mix.`=Mixture,PSW),
rmse%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1==.5)%>%
  pivot_wider(names_from=estimator,values_from=rmse)%>%
select(`M-Est`,`Mix.`=Mixture,PSW))%>%
  kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
  add_header_above(c(" " = 3, "$\\\\alpha=0$" = 3, "$\\\\alpha=0.2$" = 3, "$\\\\alpha=0.5" = 3),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")
sink()



#######################################################
#### coverage
#######################################################
 coverage <- pd%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,interactionS,interactionZ,intS,intZ,eff,estimator)%>%
   summarize(coverage=mean(CInormL<=pop & CInormU>=pop,na.rm=TRUE))%>%ungroup()

coverage%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3)%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+geom_hline(yintercept=0.95)+
  facet_grid(interactionS~errDist+interactionZ)

sink('writeUps/coverageTab.tex')
cbind(
  coverage%>%
  filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1==0)%>%
  transmute(`Residual\nDist.`=c(norm='Normal',unif='Uniform')[errDist],
         `X:Z\nInt.?`=ifelse(intZ,'Yes','No'),
         `X:S\nInt.?`=ifelse(intS,'Yes','No'),
         estimator=ifelse(estimator=='Mixture','Mix.',estimator),coverage)%>%
pivot_wider(names_from=estimator,values_from=coverage)
,
coverage%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1==.2)%>%
  pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`M-Est`,`Mix.`=Mixture),
coverage%>%filter(n==500,eff==1,errDist!='mix',mu01==0.3,b1==.5)%>%
  pivot_wider(names_from=estimator,values_from=coverage)%>%
select(`M-Est`,`Mix.`=Mixture))%>%
  kbl('latex',booktabs=TRUE,col.names=linebreak(names(.)),escape=FALSE,digits=2)%>%
  add_header_above(c(" " = 3, "$\\\\alpha=0$" = 2, "$\\\\alpha=0.2$" = 2, "$\\\\alpha=0.5" = 2),escape=FALSE)%>%
  collapse_rows(columns=1,latex_hline="major",valign="middle")
sink()

pdf('simFigs/coverage.pdf',height=9,width=6.5)

### normal errors, no interaction
coverage%>%
  filter(errDist=='norm',!intS,!intZ)%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0.95)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')%>%print()


### non-normal errors, no interaction
coverage%>%
  filter(errDist!='norm',eff=='1',!intS,!intZ)%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(errDist~mu01+n,scales='free')+
  ggtitle('non-normal errors, no interaction stratum 1')%>%print()


#### interactions, normal errors
coverage%>%
  filter(errDist=='norm',eff=='1')%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'normal errors, stratum 1, interactions')%>%print()


#### interactions, mixture errors
coverage%>%
  filter(errDist=='mix',eff=='1')%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'mixture errors, stratum 1, interactions')%>%print()


#### interactions, uniform errors
coverage%>%
  filter(errDist=='unif',eff=='1')%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'uniform errors, stratum 1, interactions')%>%print()

dev.off()

wSE <- results%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,intS,intZ,eff,estimator)%>%
  summarize(SE2=mean(se^2,na.rm=TRUE),Vrep=var(est,na.rm = TRUE),
            vBias=SE2-Vrep,se=sqrt(SE2),sdRep=sqrt(Vrep),
            seBias=se-sdRep)


### normal errors, no interaction
SE%>%
  filter(b1>0,errDist=='norm',!intS,!intZ,b1>0)%>%
  ggplot(aes(b1,seBias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')

### non-normal errors, no interaction
SE%>%
  filter(b1>0,errDist!='norm',eff=='1',!intS,!intZ)%>%
  ggplot(aes(b1,seBias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(errDist~mu01+n,scales='free')+
  ggtitle('non-normal errors, no interaction stratum 1')


#### interactions, normal errors
SE%>%
  filter(b1>0,errDist=='norm',eff=='1')%>%
  ggplot(aes(b1,seBias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'normal errors, stratum 1, interactions')


#### interactions, mixture errors
SE%>%
  filter(b1>0,errDist=='mix',eff=='1')%>%
  ggplot(aes(b1,seBias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'mixture errors, stratum 1, interactions')


#### interactions, uniform errors
SE%>%
  filter(b1>0,errDist=='unif',eff=='1')%>%
  ggplot(aes(b1,seBias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'uniform errors, stratum 1, interactions')


#######################################################
#### error rates
#######################################################
rate <- results%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,intS,intZ,eff,estimator)%>%
  summarize(rate=mean(abs(est)>=2*se,na.rm=TRUE),
            powerORlevel=ifelse(pop==0,'Level',
                                paste0('Power (eff=',pop,')'))
  )

pdf('simFigs/rates.pdf',height=6.5,width=6.5)

### normal errors, no interaction
rate%>%
  filter(errDist=='norm',!intS,!intZ,eff!='0')%>%
  ggplot(aes(b1,rate,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0.05)+
  facet_grid(powerORlevel~eff+n,scales='free')+
  ggtitle( 'normal errors, no interaction')%>%print()


### non-normal errors, no interaction
rate%>%
  filter(errDist!='norm',eff=='1',!intS,!intZ)%>%
  ggplot(aes(b1,rate,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(powerORlevel~errDist+n,scales='free')+
  ggtitle('non-normal errors, no interaction stratum 1')%>%print()


#### interactions, normal errors
rate%>%
  filter(errDist=='norm',eff=='1')%>%
  ggplot(aes(b1,rate,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(powerORlevel~paste('intS=\n',intS,'\nintZ=\n',intZ)+n,scales='free')+
  ggtitle( 'normal errors, stratum 1, interactions')%>%print()


#### interactions, mixture errors
rate%>%
  filter(errDist=='mix',eff=='1')%>%
  ggplot(aes(b1,rate,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(powerORlevel~paste('intS=\n',intS,'\nintZ=\n',intZ)+n,scales='free')+
  ggtitle( 'mixture errors, stratum 1, interactions')%>%print()


#### interactions, uniform errors
rate%>%
  filter(errDist=='unif',eff=='1')%>%
  ggplot(aes(b1,rate,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept=0.05)+
  facet_grid(powerORlevel~paste('intS=\n',intS,'\nintZ=\n',intZ)+n,scales='free')+
  ggtitle( 'uniform errors, stratum 1, interactions')%>%print()

dev.off()



pd100 <- res100 %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=c(0,0.2,0.5),
                labels=c(bquote(alpha==0),bquote(alpha==0.2),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='M-Est',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction")
  ) %>%
  filter(estimator=='PSW'|rhat<1.1)
