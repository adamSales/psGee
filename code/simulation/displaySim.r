library(tidyverse)

source('code/simulation/readSimFuncs.r')
#### after simulation results have been loaded and pre-processed...


print(load('simResults/fullResults.RData'))

if(!is.data.frame(results)) results <- reduce(results,bind_rows)

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






pd <- results %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = as.factor(b1),
    AUCff=paste0('AUCf=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=paste0('disp=',mu01),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors')
  ) %>%
  filter(rhat<1.1)

### boxplots for when M-Estimation assumptions hold+normal erros


bp(pd,
   subset=b1 > 0& errDist=='norm'& !intS& !intZ&eff!='Diff',
   title="Normal Outcomes, No Interaction, b1>0",
   facet = AUCf+N~PE + M1
   ) 



### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& !intS& !intZ&errDist!='norm',
   title="Stratum 1, No Interaction, b1>0",
   facet=AUCf+N~errDist+ M1)


### now the same, but with interaction in S

bp(pd,
   subset=b1 > 0& errDist=='norm'& intS& !intZ&eff!='Diff',
   title="Normal Outcomes, Interaction in S, b1>0",
   facet = AUCf+N~PE + M1
) 

### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& intS& !intZ&errDist!='norm',
   title="Stratum 1, Interaction in S, b1>0",
   facet=AUCf+N~errDist+ M1)



### now the same, but with interaction in Z

bp(pd,
   subset=b1 > 0& errDist=='norm'& !intS& intZ&eff!='Diff',
   title="Normal Outcomes, Interaction in Z, b1>0",
   facet = AUCf+N~PE + M1
) 

### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& !intS& intZ&errDist!='norm',
   title="Stratum 1, Interaction in Z, b1>0",
   facet=AUCf+N~errDist+ M1)

### both interactions
### now the same, but with interaction in S

bp(pd,
   subset=b1 > 0& errDist=='norm'& intS& intZ&eff!='Diff',
   title="Normal Outcomes, Both Interactions , b1>0",
   facet = AUCf+N~PE + M1
) 

### now just stratum 1, but different error dists

bp(pd,
   subset=b1 > 0& eff=='1'& intS& intZ&errDist!='norm',
   title="Stratum 1, Both Interaction, b1>0",
   facet=AUCf+N~errDist+ M1)


#######################################################
### summary stats
#######################################################

#######################################################
#### bias
#######################################################
bias <- results%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,intS,intZ,eff,estimator)%>%
  summarize(bias=mean(est-pop,na.rm=TRUE),
            ttest=tidy(t.test(est-pop)))%>%
  ungroup()%>%
  bind_cols(.$ttest)%>%
  select(-ttest)



### normal errors, no interaction
bias%>%
  filter(errDist=='norm',!intS,!intZ)%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')

# with error bars:
bias%>%
  filter(errDist=='norm',!intS,!intZ,b1>0)%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width=0)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')

### non-normal errors, no interaction
bias%>%
  filter(errDist!='norm',eff=='1',!intS,!intZ)%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(errDist~mu01+n,scales='free')+
  ggtitle('non-normal errors, no interaction stratum 1')


#### interactions, normal errors
bias%>%
  filter(errDist=='norm',eff=='1')%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'normal errors, stratum 1, interactions')


#### interactions, mixture errors
bias%>%
  filter(errDist=='mix',eff=='1')%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'mixture errors, stratum 1, interactions')


#### interactions, uniform errors
bias%>%
  filter(errDist=='unif',eff=='1')%>%
  ggplot(aes(b1,bias,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'uniform errors, stratum 1, interactions')


#######################################################
#### coverage
#######################################################
coverage <- res%>%
  filter(rhat<1.1)%>%
  group_by(n,mu01,errDist,b1,intS,intZ,eff,estimator)%>%
  summarize(coverage=mean(CInormL<=pop & CInormU>=pop,na.rm=TRUE))


### normal errors, no interaction
coverage%>%
  filter(errDist=='norm',!intS,!intZ)%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept = 0.95)+
  facet_grid(eff~mu01+n,scales='free')+
  ggtitle( 'normal errors, no interaction')


### non-normal errors, no interaction
coverage%>%
  filter(errDist!='norm',eff=='1',!intS,!intZ)%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept =0.95)+
  facet_grid(errDist~mu01+n,scales='free')+
  ggtitle('non-normal errors, no interaction stratum 1')


#### interactions, normal errors
coverage%>%
  filter(errDist=='norm',eff=='1')%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept =0.95)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'normal errors, stratum 1, interactions')


#### interactions, mixture errors
coverage%>%
  filter(errDist=='mix',eff=='1')%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept =0.95)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'mixture errors, stratum 1, interactions')


#### interactions, uniform errors
coverage%>%
  filter(errDist=='unif',eff=='1')%>%
  ggplot(aes(b1,coverage,color=estimator,fill=estimator,group=estimator))+
  geom_point()+geom_line()+
  geom_hline(yintercept =0.95)+
  facet_grid(paste('intS=',intS,'\nintZ=',intZ)~mu01+n,scales='free')+
  ggtitle( 'uniform errors, stratum 1, interactions')


#######################################################
#### standard errors
#######################################################
SE <- results%>%
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




