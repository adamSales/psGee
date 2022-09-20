library(coefplot)

## ps models
coefplot(estimates3$BAU$psMod, predictors=names(model.frame(estimates3$BAU$psMod))[-c(1:2)],
         newNames=c(Scale.Score5='G5 Std. Test',pre.total_time_on_tasks='Log(Pretest Time on Task)',
                    ESOL1='ESOL',IEP1='IEP',
                    pre_MSE_total_score='Math Self Efficacy',pre_MA_total_scoreNATRUE='Baseline Missing',
                    fullYear5TRUE='Full Year 5th Grd',GenderM='Male',"raceEthHispanic/Latino"="Hispanic/Latinx",
                    raceEthAsian="Asian",raceEthOther="Other Race/Eth.",GIFTED1="Gifted",
                    "poly(pre_PS_tasks_total_score, 2, raw = TRUE)1"="Perceptual Sensitivity",
                    "poly(pre_PS_tasks_total_score, 2, raw = TRUE)2"="Perceptual Sensitivity^2"),
         title="Principal Score Model",xlab="Estimate",ylab=NULL)

ggsave('figure/psModCoef.jpg',width=5,height=3,units='in')


### plot estimates

ateDF <- imap_dfr(ates,~data.frame(EFF=.25,Alternative=ifelse(.y=='Dragon','DragonBox',.y),estimates=.x$coefficients['Z'],
                              ymin=.x$conf.low['Z'],ymax=.x$conf.high['Z']))

bind_rows(
  imap_dfr(estimates3,~effsFromFit(.x)%>%
                        data.frame()%>%
                        rownames_to_column("eff")%>%
                        mutate(int="None",Alternative=.y)),
    imap_dfr(estimatesInt1,~effsFromFit(.x)%>%
                        data.frame()%>%
                        rownames_to_column("eff")%>%
                        mutate(int="Pretest &\nG5 Test",Alternative=.y))
)%>%
  filter(startsWith(eff,'eff'))%>%
  mutate(
    Alternative=ifelse(Alternative=='Dragon','DragonBox',Alternative),
    EFF=ifelse(eff=='eff0',ifelse(int=="None",-.05,.05),ifelse(int=="None",.45,.55)),
    ymin=estimates-2*SE,
    ymax=estimates+2*SE
  )%>%
  ggplot(aes(EFF,estimates,color=int,ymin=ymin,ymax=ymax))+
  geom_point()+geom_errorbar(width=0)+
  geom_point(data=ateDF,aes(EFF,estimates),color='black',inherit.aes=FALSE)+
  geom_errorbar(data=ateDF,aes(EFF,ymin=ymin,ymax=ymax),width=0,color='black',inherit.aes=FALSE)+
  geom_hline(yintercept=0)+
  scale_x_continuous("Principal Stratum",breaks=c(0,.25,.5),minor_breaks=NULL,
                     labels=c("Non-\nBottom-\nOuter","ATE","Bottom-\nOuter"),limits=c(-.15,.65))+
  labs(color='S Interaction',y='Principal Effect')+
  facet_wrap(~Alternative,nrow=1)+
  theme(legend.position="top")
ggsave("figure/prinEffs.pdf",width=5,height=3,units="in")

#################
### comparing GEEPERS to other methods
bind_rows(
  imap_dfr(estimates3,~effsFromFit(.x)%>%
                        data.frame()%>%
                        rownames_to_column("eff")%>%
                        mutate(method="GEEPERS",Alternative=.y)),
  imap_dfr(pswResults,~.$coef%>%
                        data.frame()%>%
                        rownames_to_column("eff")%>%
                        mutate(method='PSW',Alternative=.y)%>%
                        rename(SE='se',estimates='est')),
#### need to replace this with function or something
  data.frame(eff=rep(paste0('eff',c(0,1)),3),
             estimates=c(.3,.16,-.11,-.07,-.23,-.28),
             SE=c(.28,.24,.23,.21,.24,.27),
             method='Mixture',
             Alternative=rep(alts,each=2))
)%>%
  filter(startsWith(eff,'eff'))%>%
  mutate(
    Alternative=ifelse(Alternative=='Dragon','DragonBox',Alternative),
    EFF=ifelse(eff=='eff0',ifelse(method=="GEEPERS",-.07,ifelse(method=='Mixture',0,.07)),ifelse(method=="GEEPERS",.43,ifelse(method=='Mixture',.5,.57))),
    ymin=estimates-2*SE,
    ymax=estimates+2*SE
  )%>%
  ggplot(aes(EFF,estimates,color=method,ymin=ymin,ymax=ymax))+
  geom_point()+geom_errorbar(width=0)+
  geom_hline(yintercept=0)+
  scale_x_continuous("Principal Stratum",breaks=c(0,.5),minor_breaks=NULL,
                     labels=c("Non-\nBottom-\nOuter","Bottom-\nOuter"),limits=c(-.15,.65))+
  labs(color='Method',y='Principal Effect')+
  facet_wrap(~Alternative,nrow=1)+
  theme(legend.position="top")
ggsave("figure/compareMethods.pdf",width=5,height=3,units="in")
