library(tidyverse)
library(broom)
library(kableExtra)
library(gridExtra)

source('code/simulation/readSimFuncs.r')

#### after simulation results have been loaded and pre-processed...
## source('code/simulation/readSim.r')
load('simResults/fullResults.RData')
load('simResults/resultsNs.RData')
load('simResults/resultsB1s.RData')


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


resultsB1s%>%
  group_by(b1)%>%
  mutate(meanAUC=mean(auc))%>%
  ggplot(aes(as.factor(b1),auc))+geom_boxplot()+geom_point(aes(y=meanAUC))+geom_smooth(se=FALSE)
                                        #  scale_y_continuous('Avg. AUC',seq(.5,1,.1))
resultsB1s%>%
  group_by(b1)%>%
  summarize(meanAUC=mean(auc))

resultsB1s%>%
  group_by(b1)%>%
  summarize(meanAUC=mean(auc))%>%ggplot(aes(b1,meanAUC))+geom_point()+geom_smooth(method='lm')



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
    B1 = factor(b1,levels=c(0,0.2,0.5),
                labels=c(bquote(alpha==0),bquote(alpha==0.2),bquote(alpha==0.5))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              ##labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
              labels=c(bquote(mu[c]^1==0),bquote(mu[c]^1==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='M-Est',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z\ninteraction","No Z\ninteraction"),
    interactionS=ifelse(intS,"S\ninteraction","No S\ninteraction")
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



norm<-bp(pd,
   subset=b1 > 0& errDist=='norm'& !intS& !intZ&eff!='Diff'&n==500&PE=='Stratum 1',
   title="Normal Residuals",
   facet = B1~ M1,
   ylim=c(-1,1),
   Labeller=label_parsed,labSize=2
   )+
  labs(y=bquote("Estimation Error for "~tau^1),subtitle="No Interactions",x=NULL)+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.subtitle = element_text(size=10),legend.pos="none")+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

#ggsave('simFigs/normalOutcomesNoInteractionPosB1.jpg',height=4,width=6)

unif<-bp(pd,
   subset=b1 > 0& errDist=='unif'& !intS& !intZ&eff!='Diff'&n==500&PE=='Stratum 1',
   title="Uniform Residuals",
   facet = B1~ M1,
   ylim=c(-1,1),
   Labeller=label_parsed,labSize=2
   )+
  labs(subtitle="No Interactions",x=NULL,y=NULL)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.subtitle = element_text(size=10))+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))
#ggsave('simFigs/unifOutcomesNoInteractionPosB1.jpg',height=4,width=6)

int<-bp(pd,
   subset=b1 > 0& errDist=='norm'& n==500&PE=='Stratum 1'&b1==.5&mu01==0.3,
   title="S&Z Interactions",
   facet = interactionZ~interactionS,
   ylim=NULL
   )+
  labs(subtitle=bquote(~alpha==0.5~","~mu[C]^1==0.3~"Norm. Resid"),x=NULL,y=NULL)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.subtitle = element_text(size=8))+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))

#ggsave('simFigs/InteractionPosB1.jpg',height=4,width=6)

pdf("simFigs/boxplots.pdf",width=6.4,height=4)
grid.arrange(norm,unif,int,nrow=1)
dev.off()

int.2<-bp(pd,
   subset=b1 > 0& errDist=='norm'& n==500&PE=='Stratum 1'&b1==.2&mu01==0.3,
   title="S&Z Interactions",
   facet = interactionZ~interactionS,
   ylim=c(-1.5,1.5)
   )+
  labs(subtitle=bquote(~alpha==0.5~","~mu[C]^1==0.3~"Norm. Resid"),x=NULL,y=NULL)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.subtitle = element_text(size=8))+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))



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


### n=100 results
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

###############################################################################
### results across ns
###############################################################################
 pdns <- resultsNs %>%
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
    interactionS=ifelse(intS,"S interaction","No\nS interaction"),
    N=paste0("n=",n)
  ) %>%
  filter(rhat<1.1)


bp(pdns,subset=eff==1,facet=~N,title="Estimation Error by N",ylim=c(-1,1))+
  labs(subtitle=bquote("Stratum 1 Principal Effects; Normal Residuals, No Interactions, "~alpha>0~", "~mu[C]^1-mu[C]^0==0.3))
ggsave('simFigs/byN.jpg',width=6,height=3)


rmseN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
  summarize(rmse=sqrt(mean(errP^2,na.rm=TRUE)))%>%ungroup()

rmseN%>%filter(eff==1)%>%select(-eff)%>%
  pivot_wider(names_from=n,names_prefix="n=",values_from=rmse)

rmseN%>%filter(eff==1)%>%select(-eff)%>%
  mutate(N1000=factor(ifelse(n==1000,"n=1000","99<n<501"),levels=c("99<n<501","n=1000")))%>%
  ggplot(aes(n,rmse,color=estimator,group=estimator,fill=estimator))+
  geom_point()+geom_line()+facet_grid(cols=vars(N1000),scales="free_x",space="free")

 coverage <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
   summarize(coverage=mean(CInormL<=pop & CInormU>=pop,na.rm=TRUE))%>%ungroup()

coverageN%>%filter(eff==1)%>%select(-eff)%>%
  pivot_wider(names_from=n,names_prefix="n=",values_from=coverage)

coverageN%>%filter(eff==1)%>%select(-eff)%>%
  mutate(N1000=factor(ifelse(n==1000,"n=1000","99<n<501"),levels=c("99<n<501","n=1000")))%>%
  ggplot(aes(n,coverage,color=estimator,group=estimator,fill=estimator))+
  geom_point()+geom_line()+geom_hline(yintercept=0.95)+
  facet_grid(cols=vars(N1000),scales="free_x",space="free")



biasN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
    summarize(bias=mean(est-pop,na.rm=TRUE),
            ttest=tidy(t.test(est-pop)))%>%
  ungroup()%>%
  bind_cols(.$ttest)%>%
  select(-ttest)

biasN%>%filter(eff==1)%>%select(-eff)%>%
  pivot_wider(names_from=n,names_prefix="n=",values_from=bias)

biasN%>%filter(eff==1)%>%select(-eff)%>%
  mutate(N1000=factor(ifelse(n==1000,"n=1000","99<n<501"),levels=c("99<n<501","n=1000")))%>%
  ggplot(aes(n,bias,color=estimator,group=estimator,fill=estimator))+
  geom_point()+geom_line()#+facet_grid(cols=vars(N1000),scales="free_x",space="free")

seN <- pdns%>%
  filter(rhat<1.1)%>%
  group_by(n,eff,estimator)%>%
  summarize(se=sd(errP,na.rm=TRUE))%>%ungroup()

full_join(seN,biasN)%>%filter(eff==1,estimator!="PSW")%>%summarize(across(c(se,bias),max),coef=round(se/bias))


bind_rows(
  seN%>%mutate(se=se/3,meas="SE")%>%rename(what=se),
  biasN%>%mutate(meas="Bias")%>%rename(what=bias))%>%
  filter(eff==1,estimator!="PSW")%>%
  ggplot(aes(n,what,color=estimator,group=paste0(estimator,meas),linetype=meas,fill=estimator,shape=meas))+
  geom_point()+geom_line()+geom_hline(yintercept=0)+#,linetype="dotted",size=2)+
  scale_y_continuous(name="Bias",sec.axis=sec_axis(trans=~.*3,name='Standard Error'))+
  scale_shape_manual(values=c(0,16))+
  annotate('text',600,.125,label=list(bquote(atop("Normal Resid., No Interactions, ",alpha==0.5~", "~mu[C]^1-mu[C]^0==0.3))),parse=TRUE)+
  scale_x_continuous(name="Sample Size Per Group",breaks=c(seq(100,500,200),1000))+
  ggtitle("Bias and Standard Error by n")+theme(legend.title=element_blank())
ggsave("simFigs/biasSEbyN.jpg",width=6,height=3)


###############################################################################
### results across B1
###############################################################################

 pdb1 <- resultsB1s %>%
  group_by(b1)%>%
  mutate(AUCf=mean(auc,na.rm=TRUE))%>%
  ungroup()%>%
  mutate(
    errP = est - pop,
    B1 = factor(b1,levels=seq(0,1,.1),
                labels=lapply(seq(0,1,.1), function(x) bquote(alpha==.(x)))),
    AUCff=paste0('AUC=',round(AUCf,1)),
    N=paste0('n=',n),
    M1=factor(mu01,levels=c("0","0.3"),
              labels=c(bquote(mu[c]^1-mu[c]^0==0),bquote(mu[c]^1-mu[c]^0==0.3))),
    PE=paste('Stratum',eff),
    dist=paste(c(mix='Mixture',unif='Unform',norm='Normal')[errDist],'Errors'),
    estimator=c(bayes='Mixture',mest='M-Est',psw='PSW')[estimator],
    interactionZ=ifelse(intZ,"Z interaction","No\nZ interaction"),
    interactionS=ifelse(intS,"S interaction","No\nS interaction"),
    N=paste0("n=",n)
  ) %>%
  filter(rhat<1.1)


bp(pdb1,subset=eff==1,facet=~B1,title=bquote("Estimation Error by "~alpha),Labeller=label_parsed,labSize=2)+
  labs(subtitle=bquote("Stratum 1; Normal Resid., No Interactions, "~n==500~", "~mu[C]^1-mu[C]^0==0.3))
ggsave('simFigs/byAlpha.jpg',width=6,height=3)


rmseB1 <- pdb1%>%
  filter(rhat<1.1)%>%
  group_by(b1,eff,estimator)%>%
  summarize(rmse=sqrt(mean(errP^2,na.rm=TRUE)))%>%ungroup()

rmseB1%>%filter(eff==1)%>%select(-eff)%>%
  pivot_wider(names_from=b1,names_prefix="b1=",values_from=rmse)%>%knitr::kable(format='markdown',digits=3)

rmseB1%>%filter(eff==1)%>%select(-eff)%>%
  ggplot(aes(b1,rmse,color=estimator,group=estimator,fill=estimator))+
  geom_point()+geom_line()

 coverageB1 <- pdb1%>%
  filter(rhat<1.1)%>%
  group_by(b1,eff,estimator)%>%
   summarize(coverage=mean(CInormL<=pop & CInormU>=pop,na.rm=TRUE))%>%ungroup()

coverageB1%>%filter(eff==1)%>%select(-eff)%>%
  pivot_wider(names_from=b1,names_prefix="n=",values_from=coverage)

coverageB1%>%filter(eff==1)%>%select(-eff)%>%
  ggplot(aes(b1,coverage,color=estimator,group=estimator,fill=estimator))+
  geom_point()+geom_line()+geom_hline(yintercept=0.95)


biasB1 <- pdb1%>%
  filter(rhat<1.1)%>%
  group_by(b1,eff,estimator)%>%
  summarize(bias=mean(errP,na.rm=TRUE))%>%ungroup()

biasB1%>%filter(eff==1)%>%select(-eff)%>%
  pivot_wider(names_from=b1,names_prefix="n=",values_from=bias)

biasB1%>%filter(eff==1)%>%select(-eff)%>%
  ggplot(aes(b1,bias,color=estimator,group=estimator,fill=estimator))+
  geom_point()+geom_line()#+facet_grid(cols=vars(N1000),scales="free_x",space="free")

seB1 <- pdb1%>%
  filter(rhat<1.1)%>%
  group_by(b1,eff,estimator)%>%
  summarize(se=sd(errP,na.rm=TRUE))%>%ungroup()

full_join(seB1,biasB1)%>%filter(eff==1,estimator!="PSW",b1>0)%>%summarize(across(c(se,bias),max),coef=round(se/bias))

bind_rows(
  seB1%>%mutate(se=se/5,meas="SE")%>%rename(what=se),
  biasB1%>%mutate(meas="Bias")%>%rename(what=bias))%>%
  filter(b1>0,eff==1,estimator!="PSW")%>%
  ggplot(aes(b1,what,color=estimator,group=paste0(estimator,meas),linetype=meas,fill=estimator,shape=meas))+
  geom_point()+geom_line()+geom_hline(yintercept=0)+#,linetype="dashed")+
  scale_y_continuous(name="Bias",sec.axis=sec_axis(trans=~.*5,name='Standard Error'))+
  scale_shape_manual(values=c(0,16))+
  annotate('text',.7,.1,label=list(bquote(atop("Normal Resid., No Interactions, ",n==500~", "~mu[C]^1-mu[C]^0==0.3))),parse=TRUE)+xlab(bquote(alpha))+ggtitle(bquote("Bias and Standard Error by "~alpha))
ggsave("simFigs/biasSEbyB1.jpg",width=6,height=3)
