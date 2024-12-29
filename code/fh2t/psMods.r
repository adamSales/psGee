library(arm)

load('data/psdat.RData')

psdat <- droplevels(filter(psdat,SchIDPre!=7))

auc <- function(x,y)
  if(length(unique(y))==2 & length(y)==length(x)){
    wilcox.test(x~y)$statistic/(sum(y)*(sum(1-y)))
  }else
    wilcox.test(x,y)$statistic/(legnth(x)*length(y))

aucMod <- function(mod)
  auc(mod$linear,1-mod$y)


psMod1 <- glm(
  S~.-StuID-TeaIDPre-ClaIDPre-S-Y-virtual,
  data=psdat%>%filter(trt=='ASSISTments')%>%select(-trt),
  family=binomial)

aucMod(psMod1)
boxplot(fitted(psMod1)~psMod1$y)




summary(psMod2 <- step(psMod1))
aucMod(psMod2)

psMod3 <- update(psMod2,.~.-SchIDPre+TeaIDPre)
aucMod(psMod3)

### what about lumping teachers with fewer than 10 kids, within school?
psdat <- psdat%>%
  group_by(SchIDPre)%>%
  mutate(maxTch=levels(TeaIDPre)[which.max(table(TeaIDPre))])%>%
  group_by(TeaIDPre)%>%
  mutate(nTrt=sum(trt=='ASSISTments'))%>%
  ungroup()%>%
  mutate(tid2=factor(ifelse(nTrt>9,as.character(TeaIDPre),maxTch)))

psMod4 <- update(psMod3,.~.-TeaIDPre+tid2)
aucMod(psMod4)

### but is it just overfit?
### 10-fold cv
cv <- function(mod,folds=10,seed=613){
  mf=mod$data[names(mod$y),]
  set.seed(seed)
  folds <- seq(nrow(mf))%>%sample()%>%cut(10,labels=FALSE)
  cvPred <- numeric(nrow(mf))
  for(ff in 1:10)
    cvPred[folds==ff] <- predict(update(mod,data=mf[folds!=ff,]),mf[folds==ff,])
  AUC <-auc(cvPred,1-mod$y)
  print(AUC)
  invisible(list(cvPred=cvPred,folds=folds,auc=AUC))
}

cv(psMod2)
cv(psMod4)

#### binned plots for psMod2
binnedplot(fitted(psMod2),resid(psMod2,type='response'))

### by (continuous) covariate
par(mfrow=c(3,3))
model.frame(psMod2)%>%select(where(is.numeric))%>%
  iwalk(~binnedplot(.x,y=resid(psMod2,type='response'),main=.y))
par(mfrow=c(1,1))

### curvature in PS score?
binnedplot(model.frame(psMod2)$pre_PS_tasks_total_score,resid(psMod2,type='response'),main="pre_PS_tasks_total_score",ylim=c(-.5,.5))

summary(psMod5 <- update(psMod2,.~.-pre_PS_tasks_total_score+splines::ns(pre_PS_tasks_total_score,3)))

#### binned plots for psMod5
binnedplot(fitted(psMod5),resid(psMod5,type='response'))

### by (continuous) covariate
par(mfrow=c(3,3))
model.frame(psMod2)%>%select(where(is.numeric))%>%
  iwalk(~binnedplot(.x,y=resid(psMod5,type='response'),main=.y))
par(mfrow=c(1,1))

par(mfrow=c(4,5))
model.frame(psMod1)%>%select(!(S:Y))%>%select(where(~is.numeric(.)&n_distinct(.)>2))%>%
  iwalk(~binnedplot(.x,y=resid(psMod5,type='response'),main=.y))
par(mfrow=c(1,1))

model.frame(psMod1)%>%select(!(S:Y))%>%select(!where(~is.numeric(.)&n_distinct(.)>2))%>%
  mutate(resid=resid(psMod5,type='response'))%>%
  mutate(across(-resid,as.factor))%>%
  pivot_longer(-resid,names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),sdResid=sd(resid),up=meanResid+2*sdResid,down=meanResid-2*sdResid)%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down))+geom_point()+geom_errorbar(width=0)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")


model.frame(psMod1)%>%select(ends_with("NA"))%>%rowSums()%>%binnedplot(resid(psMod5,type='response'))

aucMod(psMod5)
cv5 <- cv(psMod5)

binnedplot(plogis(cv5$cvPred),psMod5$y-plogis(cv5$cvPred))

nd <- model.frame(psMod5)[rep(1,100),]
nd$pre_PS_tasks_total_score <- seq(0,16,length=100)
pred <- predict(psMod5,nd,se.fit=TRUE)
up <- pred$fit+2*pred$se.fit
down <- pred$fit-2*pred$se.fit


plot(seq(0,15,length=100),pred$fit,xlab="PS score",type='l',ylim=range(c(up,down)))
lines(seq(0,15,length=100),up,lty=2)
lines(seq(0,15,length=100),down,lty=2)

save(psMod5,file='results/psModGlm.RData')


#### hmmmm maybe want to add in more variables for the sake of the writeup
summary(psMod6 <- update(psMod5,
                         .~.#+Gender+raceEth
                         +GIFTED+pre_MA_total_score))

aucMod(psMod6)
cv6 <- cv(psMod6)

binnedplot(fitted(psMod6),resid(psMod6,type='response'))
par(mfrow=c(4,5))
model.frame(psMod1)%>%select(!(S:Y))%>%select(where(~is.numeric(.)&n_distinct(.)>2))%>%
  iwalk(~binnedplot(.x,y=resid(psMod6,type='response'),main=.y))
par(mfrow=c(1,1))

model.frame(psMod1)%>%select(!(S:Y))%>%select(!where(~is.numeric(.)&n_distinct(.)>2))%>%
  mutate(resid=resid(psMod6,type='response'))%>%
  mutate(across(-resid,as.factor))%>%
  pivot_longer(-resid,names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),sdResid=sd(resid),up=meanResid+2*sdResid,down=meanResid-2*sdResid)%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down))+geom_point()+geom_errorbar(width=0)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")


binnedplot(plogis(cv6$cvPred),psMod6$y-plogis(cv6$cvPred))

nd <- model.frame(psMod6)[rep(1,100),]
nd$pre_PS_tasks_total_score <- seq(0,16,length=100)
pred <- predict(psMod6,nd,se.fit=TRUE)
up <- pred$fit+2*pred$se.fit
down <- pred$fit-2*pred$se.fit


plot(seq(0,15,length=100),pred$fit,xlab="PS score",type='l',ylim=range(c(up,down)))
lines(seq(0,15,length=100),up,lty=2)
lines(seq(0,15,length=100),down,lty=2)

save(psMod6,file='results/psModGlm6.RData')

## how different are they, really?
ps5 <- predict(psMod5,psdat)
ps6 <- predict(psMod6,psdat)

plot(ps5,ps6,col=ifelse(psdat$Z==1,'blue','red'))
cor(ps5,ps6)

ps5 <- predict(psMod5,psdat,type='response')
ps6 <- predict(psMod6,psdat,type='response')

plot(ps5,ps6,col=ifelse(psdat$Z==1,'blue','red'))
cor(ps5,ps6)

#### OK actually maybe use a quadratic polynomial for PS instead of spline, for interpretability
summary(psMod7 <- update(psMod6,
                 .~.-splines::ns(pre_PS_tasks_total_score,3)+poly(pre_PS_tasks_total_score,2,raw=TRUE)))

aucMod(psMod7)
cv7 <- cv(psMod7)

binnedplot(fitted(psMod7),resid(psMod7,type='response'))
par(mfrow=c(4,5))
model.frame(psMod1)%>%select(!(S:Y))%>%select(where(~is.numeric(.)&n_distinct(.)>2))%>%
  iwalk(~binnedplot(.x,y=resid(psMod7,type='response'),main=.y))
par(mfrow=c(1,1))

model.frame(psMod1)%>%select(!(S:Y))%>%select(!where(~is.numeric(.)&n_distinct(.)>2))%>%
  mutate(resid=resid(psMod7,type='response'))%>%
  mutate(across(-resid,as.factor))%>%
  pivot_longer(-resid,names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),sdResid=sd(resid),up=meanResid+2*sdResid,down=meanResid-2*sdResid)%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down))+geom_point()+geom_errorbar(width=0)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")


binnedplot(plogis(cv7$cvPred),psMod7$y-plogis(cv7$cvPred))

nd <- model.frame(psMod7)[rep(1,100),]
nd$pre_PS_tasks_total_score <- seq(min(psdat$pre_PS_tasks_total_score),max(psdat$pre_PS_tasks_total_score),length=100)
pred <- predict(psMod7,nd,se.fit=TRUE)
up <- pred$fit+2*pred$se.fit
down <- pred$fit-2*pred$se.fit


plot(seq(0,15,length=100),pred$fit,xlab="PS score",type='l',ylim=range(c(up,down)))
lines(seq(0,15,length=100),up,lty=2)
lines(seq(0,15,length=100),down,lty=2)

save(psMod7,file='results/psModGlm7.RData')

## how different are they, really?
ps5 <- predict(psMod5,psdat)
ps7 <- predict(psMod7,psdat)

plot(ps5,ps7,col=ifelse(psdat$Z==1,'blue','red'))
cor(ps5,ps7)

ps5 <- predict(psMod5,psdat,type='response')
ps7 <- predict(psMod7,psdat,type='response')

plot(ps5,ps7,col=ifelse(psdat$Z==1,'blue','red'))
cor(ps5,ps7)


### what if we just tried everything?
psModAll <- update(psMod1,.~.-pre_PS_tasks_total_score+poly(pre_PS_tasks_total_score,2,raw=TRUE))

anova(psMod7,psMod1,psModAll,test='LRT')

plot(fitted(psMod7),fitted(psModAll))
abline(0,1)

binnedplot(fitted(psModAll),resid(psModAll,type='response'))

model.frame(psMod1)%>%select(-S)%>%select(!where(~is.numeric(.)&n_distinct(.)>2))%>%
  mutate(resid=resid(psModAll,type='response'))%>%
  mutate(across(-resid,as.factor))%>%
  pivot_longer(-resid,names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),sdResid=sd(resid),up=meanResid+2*sdResid,down=meanResid-2*sdResid)%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down))+geom_point()+geom_errorbar(width=0)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")

model.frame(psMod1)%>%select(-S)%>%select(where(~is.numeric(.)&n_distinct(.)>2&n_distinct(.)<20))%>%
  mutate(resid=resid(psModAll,type='response'),pred=predict(psModAll,type='response'))%>%
  pivot_longer(-c(resid,pred),names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),n=n(),sdResid=ifelse(n>1,sd(resid),sqrt(pred*(1-pred))),up=2*sdResid,down=-2*sdResid)%>%
  ungroup()%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down,size=n))+geom_point()+#geom_errorbar(width=0)+
  geom_line(aes(y=up),color='grey',size=1)+geom_line(aes(y=down),color='grey',size=1)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")

par(mfrow=c(3,5))
model.frame(psMod1)%>%select(-S,-StuID)%>%select(where(~is.numeric(.)&n_distinct(.)>=20))%>%
  iwalk(~binnedplot(.x,resid(psModAll,type='response'),sub=.y))
par(mfrow=c(1,1))



    psModRest <- update(psModAll,.~.-pre.total_math_score-Scale.Score5-Performance.Level5-pre.sub_P_score-pre.sub_F_score-pre.math_completed_num)

anova(psModRest,psModAll,test='LRT')

model.frame(psMod1)%>%select(-S)%>%select(!where(~is.numeric(.)&n_distinct(.)>2))%>%
  mutate(resid=resid(psModRest,type='response'))%>%
  mutate(across(-resid,as.factor))%>%
  pivot_longer(-resid,names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),sdResid=sd(resid),up=meanResid+2*sdResid,down=meanResid-2*sdResid)%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down))+geom_point()+geom_errorbar(width=0)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")

model.frame(psMod1)%>%select(-S,-Y)%>%select(where(~is.numeric(.)&n_distinct(.)>2&n_distinct(.)<20))%>%
  mutate(resid=resid(psModRest,type='response'),pred=predict(psModRest,type='response'))%>%
  pivot_longer(-c(resid,pred),names_to='var',values_to='level')%>%
  group_by(var,level)%>%
  summarize(meanResid=mean(resid),n=n(),sdResid=ifelse(n>1,sd(resid),sqrt(pred*(1-pred))),up=2*sdResid,down=-2*sdResid)%>%
  ungroup()%>%
  ggplot(aes(level,meanResid,ymax=up,ymin=down,size=n))+geom_point()+#geom_errorbar(width=0)+
  geom_line(aes(y=up),color='grey',size=1)+geom_line(aes(y=down),color='grey',size=1)+
  geom_hline(yintercept=0)+facet_wrap(~var,scales="free")

par(mfrow=c(3,5))
model.frame(psMod1)%>%select(-S,-StuID)%>%select(where(~is.numeric(.)&n_distinct(.)>=20))%>%
  iwalk(~binnedplot(.x,resid(psModRest,type='response'),sub=.y))
par(mfrow=c(1,1))

anova(psModRest,update(psModRest,.~.+I(pre.sub_F_score^2)),test="LRT")


save(psModAll,psModRest,file='results/biggerPSmods.RData')
