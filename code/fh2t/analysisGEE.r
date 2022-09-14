library(tidyverse)
library(arm)
library(randomForest)
library(estimatr)
source('code/regression.r')
select <- dplyr::select

load('data/psdat.RData')


load('results/psModGlm7.RData')

psdat <- filter(psdat,SchIDPre!=7) ### only one student from school 7, in FH2T

psdat$Z <- ifelse(psdat$trt=='ASSISTments',1,0)

alts <- unique(psdat$trt[psdat$Z==0])
alts <- setNames(alts,alts)




getDat <- function(alt){
  dat <- psdat%>%
    filter(trt%in%c('ASSISTments',alt))
 # dat$pred <- predict(remnantMods[[alt]],dat)

  droplevels(dat)
}


ates <- sapply(setNames(alts,alts),
                  function(alt)
                    estimatr::lm_robust(
                                Y~Z+Scale.Score5+ESOL+IEP+pre.total_time_on_tasks+pre_MSE_total_score+
                            fullYear5+pre_PS_tasks_total_score+raceEth+Gender+GIFTED+
                            pre.total_math_score+pre.total_math_scoreNA+ClaIDPre,
                            data=getDat(alt),fixed_effects=~ClaIDPre),
                  simplify=FALSE)

estimates3 <- lapply(setNames(alts,alts),
                    function(alt)
                      est(getDat(alt),
                          covFormU=formula(psMod7)[-2],
                          covFormY=~Scale.Score5+ESOL+IEP+pre.total_time_on_tasks+pre_MSE_total_score+
                            fullYear5+pre_PS_tasks_total_score+raceEth+Gender+GIFTED+
                            pre.total_math_score+pre.total_math_scoreNA+ClaIDPre
                                          ))

estimatesInt1 <- lapply(setNames(alts,alts),
                    function(alt)
                      est(getDat(alt),
                          covFormU=formula(psMod7)[-2],
                          covFormY=~Scale.Score5+ESOL+IEP+pre.total_time_on_tasks+pre_MSE_total_score+
                            fullYear5+pre_PS_tasks_total_score+raceEth+Gender+GIFTED+
                            pre.total_math_score+pre.total_math_scoreNA+ClaIDPre,
                          intS=c("Scale.Score5","pre.total_math_score")))


estimates0 <- lapply(alts,
                    function(alt)
                      est(getDat(alt),
                          covFormU=formula(psMod7)[-2])
                    )

sapply(estimates0,effsFromFit,simplify=F)

pdf('results/outcomeModels0.pdf')
iwalk(estimates0,~plot(.x$outMod,which=1,main=.y))
dev.off()

estimates1 <- lapply(setNames(alts,alts),
                    function(alt)
                      est(getDat(alt),
                          covFormU=formula(psMod7)[-2],
                          covFormY=update(formula(psMod7)[-2],.~.-SchIDPre),
                          block='ClaIDPre')
                    )

sapply(estimates1,effsFromFit,simplify=F)

pdf('results/outcomeModels1.pdf')
iwalk(estimates1,
      ~print(ggplot(aes(fitted(.x$outMod),resid(.x$outMod)))+geom_jitter(height=0.3)+ggtitle(.y)))#",main=.y)))
dev.off()

estimates2 <- lapply(setNames(alts,alts),
                    function(alt)
                      est(getDat(alt),#droplevels(filter(psdat,trt%in%c('ASSISTments',alt))),
                          covFormU=formula(psMod7)[-2],
                          covFormY=update(formula(psMod7)[-2],.~.-SchIDPre+pred+
                                                                splines::ns(pre.total_math_score,3)+
                                                                raceEth+
                                                                Gender+
                                                                GIFTED
                                          ),
                          block='ClaIDPre')
                    )


iwalk(estimates2,
      ~print(ggplot(mapping=aes(fitted(.x$outMod),resid(.x$outMod)))+geom_jitter(height=0.3)+ggtitle(.y)))#",main=.y)))

for(alt in names(estimates2))
print(getDat(alt)%>%
  select(-(StuID:Y),-Z)%>%select(where(is.numeric)&where(~n_distinct(.)>2))%>%
  imap_dfr(~data.frame(x=.x,resid=resid(estimates2[[alt]]$outMod),predictor=.y))%>%
  ggplot(aes(x,resid))+geom_jitter(height=0.3)+geom_hline(yintercept=0)+geom_smooth()+facet_wrap(~predictor,scales="free")+ggtitle(alt))


imap_dfr(estimates2,~data.frame(fitted=fitted(.x$outMod),resid=resid(.x$outMod),alt=.y))%>%
  ggplot(aes(fitted,resid))+geom_jitter(height=0.3)+geom_hline(yintercept=0)+geom_smooth()+facet_wrap(~alt,scales="free")




imap_dfr(estimates3,~data.frame(fit=fitted(.x$outMod),res=resid(.x$outMod),alt=.y))%>%
  ggplot(aes(fit,res))+geom_jitter(height=0.3)+geom_smooth()+facet_wrap(~alt)

for(alt in names(estimates3))
print(getDat(alt)%>%
  select(-(StuID:Y),-Z)%>%select(where(is.numeric)&where(~n_distinct(.)>2))%>%
  imap_dfr(~data.frame(x=.x,resid=resid(estimates3[[alt]]$outMod),predictor=.y))%>%
  ggplot(aes(x,resid))+geom_jitter(height=0.3)+geom_hline(yintercept=0)+geom_smooth()+facet_wrap(~predictor,scales="free")+ggtitle(alt))

par(mfrow=c(3,1))
for(alt in alts){
  mod <- estimates3[[alt]]$outMod
  nd <- model.frame(mod)[rep(1,100),]
  nd$pre.total_math_score <- seq(0,10,length=100)
  pred <- predict(mod,nd,se.fit=TRUE)
  up <- pred$fit+2*pred$se.fit
  down <- pred$fit-2*pred$se.fit


  plot(seq(0,15,length=100),pred$fit,xlab="PS score",type='l',ylim=range(c(up,down)),main=alt)
  lines(seq(0,15,length=100),up,lty=2)
  lines(seq(0,15,length=100),down,lty=2)
}
par(mfrow=c(1,1))

aics <- NULL
for(i in 0:3)
  aics <- rbind(aics,map_dbl(get(paste0('estimates',i),envir=.GlobalEnv),~AIC(.$outMod)))



 sapply(estimatesInt1,effsFromFit,simplify=F)

### test standard errors w bootstrap
bsSchool <- function(dat){
  datS <- split(dat,dat$SchIDPre)
  do.call("rbind",
          lapply(datS,function(x) x[sample(1:nrow(x),nrow(x),replace=TRUE),]))
}

bsInt1 <- lapply(setNames(alts,alts),
                 function(alt){
                   print(alt)
                   datAlt <- getDat(alt)
                   replicate(5000, {
                     bsDat <- datAlt[sample(1:nrow(datAlt),nrow(datAlt),replace=TRUE),]
                     if(!any(with(bsDat,table(SchIDPre,Z))==0))
                       effsFromFit(
                         est(bsDat,#bsSchool(datAlt),
                             covFormU=formula(psMod7)[-2],
                             covFormY=~Scale.Score5+ESOL+IEP+pre.total_time_on_tasks+
                               pre_MSE_total_score+
                               fullYear5+pre_PS_tasks_total_score+raceEth+Gender+GIFTED+
                               pre.total_math_score+pre.total_math_scoreNA+ClaIDPre,
                             intS=c("Scale.Score5","pre.total_math_score")))})})

bsInt1a <- lapply(bsInt1,function(x) do.call("abind",list(x,along=0)))
save(bsInt1,file='results/bsInt1.RData')

par(mfrow=c(3,3))
for(alt in alts)
  for(ee in c('eff0','eff1','diff'))
    hist(bsInt1[[alt]][ee,1,],main=paste(alt, ee))
par(mfrow=c(1,1))

lapply(alts, function(alt) cbind(effsFromFit(estimatesInt1[[alt]])[,'estimates'],
                                 rowMeans(bsInt1[[alt]][,'estimates',])))


bsSE1 <- lapply(bsInt1a,function(x) apply(x[,,1],2,sd))

lapply(alts, function(alt) cbind(effsFromFit(estimatesInt1[[alt]])[,'SE'],bsSE1[[alt]]))

bsSE1comp <- lapply(bsInt1,function(x) cbind(apply(x[,2,]^2,1,mean),apply(x[,1,],1,var)))



### robustness checks
probit=glm(
  S ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ log(pre.avg_time_on_tasksIMP)+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(dat,Z==1),family=binomial('probit'))
estimate1.1 <- est(dat, psMod=probit,
                covFormY=
                  update(formula(mod4imp)[-2],
                         .~.-virtualIMP),block = "class")
plot(estimate1.1$outMod,which=1)
effsFromFit(estimate1.1)

bayes=bayesglm(
  S ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ log(pre.avg_time_on_tasksIMP)+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(dat,Z==1),family=binomial)
estimate1.2 <- est(dat, psMod=bayes,
                covFormY=
                  update(formula(mod4imp)[-2],
                         .~.-virtualIMP),block = "class")
plot(estimate1.2$outMod,which=1)
effsFromFit(estimate1.2)


library(lspline)
### 3-df linear splines on pretest
mod5 <-     update(mod4imp,.~.-pretestIMP-Scale.Score5IMP-log(pre.avg_time_on_tasksIMP)-  pre_MA_total_scoreIMP- pre_negative_reaction_scoreIMP- pre_numerical_confindence_scoreIMP+qlspline(pretestIMP,3)+qlspline(Scale.Score5IMP,3)+qlspline(log(pre.avg_time_on_tasksIMP),3)+qlspline(  pre_MA_total_scoreIMP,3)+qlspline( pre_negative_reaction_scoreIMP,3)+qlspline( pre_numerical_confindence_scoreIMP,3))

summary(mod5)
AUCmod(mod5)

binnedplot(mod5$fitted.values,resid(mod5,type='response'))

estimate2 <- est(dat,covFormU = formula(mod5)[-2],
                covFormY=
                  update(formula(mod5)[-2],
                         .~.-virtualIMP),block = "class")
plot(estimate2$outMod,which=1)
effsFromFit(estimate2)

#### without pretests
mod6 <-     update(mod4imp,.~.-pretestIMP-Scale.Score5IMP)

summary(mod6)
AUCmod(mod6)

binnedplot(mod6$fitted.values,resid(mod6,type='response'))

estimate3 <- est(dat,covFormU = formula(mod6)[-2],
                covFormY=
                  update(formula(mod6)[-2],
                         .~.-virtualIMP),block = "class")
plot(estimate3$outMod,which=1)
effsFromFit(estimate3)

#### random forest ps model
rf=randomForest(  anyBottom ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(dat3,Z=='ASSISTments'))

psRF=predict(rf,newdata=dat)

plot(predict(mod4imp,newdata=dat,type='response'),psRF,col=ifelse(dat$Z==1,'blue','red'))
abline(0,1)

boxplot(psRF[dat$Z==1]~dat$S[dat$Z==1])

rfDat=dat
rfDat$Sp=ifelse(rfDat$Z==1,rfDat$S,psRF)
bbb=coef(lm(update(formula(mod5)[-2],
                         Y~.-virtualIMP+class+Z*Sp),data=rfDat))
with(as.list(bbb),
                      list(
                          eff0=Z,
                          eff1=Z+`Z:Sp`,
                          diff=`Z:Sp`))

rfbs=replicate(1000,
{
    ddd=dat[sample(1:nrow(dat),nrow(dat),replace=TRUE),]
    rf=randomForest(  S~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(ddd,Z==1))

    psRF=predict(rf,newdata=ddd)


    ddd$Sp=ifelse(ddd$Z==1,ddd$S,psRF)
    bbb=coef(lm(update(formula(mod5)[-2],
                       Y~.-virtualIMP+class+Z*Sp),data=ddd))
    with(as.list(bbb),
                      c(
                          eff0=Z,
                          eff1=Z+`Z:Sp`,
                          diff=`Z:Sp`))
    })




### test standard errors
bs0=replicate(1000,
              effsFromFit(
                est(
                  dat[sample(1:nrow(dat),nrow(dat),replace=TRUE),],
                  covFormU = formula(mod4imp)[-2])))

bs1=replicate(1000,
              effsFromFit(
                est(
                  dat[sample(1:nrow(dat),nrow(dat),replace=TRUE),],
                  covFormU = formula(mod4imp)[-2],
                   covFormY=update(formula(mod4imp)[-2],
                                   .~.-virtualIMP),block = "class")))


#### model checking
estEffs=function(fdat){
  fest=est(fdat,covFormU = formula(mod4imp)[-2],
                covFormY=
                  update(formula(mod4imp)[-2],
                         .~.-virtualIMP),block = "class")
  plot(fest$outMod,which=1)
  effsFromFit(fest)
}
fakeDat=read.csv('data/FakeDataForCheck.csv')

estEffs(fakeDat)


addEff=function(fdat,eff1,eff0)
  within(fdat, Y <- ifelse(Z==0,Y,
                           ifelse(S==1,Y+eff1,Y+eff0)))

### what if there were an effect of 0.2 SDs?
fakeDat2=addEff(fakeDat,0.2,0.2)
#fakeDat2$Y=fakeDat2$Y+0.2*fakeDat2$Z

estEffs(fakeDat2)

### what if the truth was like we estimated?
fakeDat3=addEff(fakeDat,0.05,-0.05)

estEffs(fakeDat3)

### randomly varying effects
fakeDat4=addEff(fakeDat,rnorm(nrow(fakeDat4),0.2,0.1),0)
estEffs(fakeDat4)
