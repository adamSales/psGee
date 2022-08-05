library(tidyverse)
library(arm)
library(randomForest)
source('code/regression.r')
select <- dplyr::select

load('data/psdatBAU.RData')
load('data/psdat.RData')

if(!exists('alt')) alt <- 'bau'

dat <- if(alt=='bau') psdatBAU else psdat

mod4imp <- glm(
  anyBottom ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ log(pre.avg_time_on_tasksIMP)+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(dat3,Z=='ASSISTments'),family=binomial)

summary(mod4imp)
AUCmod(mod4imp)

binnedplot(mod4imp$fitted.values,resid(mod4imp,type='response'))



estimate0 <- est(dat,covFormU = formula(mod4imp)[-2])
                # covFormY=
                #   update(formula(mod4imp)[-2],
                #          .~.-virtualIMP),block = "class")
plot(estimate0$outMod,which=1)
effsFromFit(estimate0)

estimate1 <- est(dat,covFormU = formula(mod4imp)[-2],
                covFormY=
                  update(formula(mod4imp)[-2],
                         .~.-virtualIMP),block = "class")
plot(estimate1$outMod,which=1)
effsFromFit(estimate1)


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
