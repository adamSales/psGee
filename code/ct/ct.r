library(tidyverse)
library(rstan)
library(geex)

source('code/ctMest.r')
source('code/ctStan.r')

#load('../data/problemLevelUsageData/probLevelData.RData')
load('../data/RANDstudyData/HSdata.RData')
hs <- dat

hs$xirt <- hs$xirt/ifelse(hs$year==1,sd(hs$xirt[hs$year==1]),sd(hs$xirt[hs$year==2]))

load('../data/sectionLevelUsageData/advanceDataWInfo.RData')


numCP <- advance%>%group_by(field_id,year)%>%summarize(ncp=sum(status=='changed placement',na.rm=TRUE))

hs <- left_join(hs,numCP)

hs <- hs%>%group_by(pair,year)%>%
    mutate(na=mean(is.na(ncp[treatment==1])))%>%
    filter(na<0.5)%>%
    ungroup()


hs$everCP <- hs$ncp>0
hs$everCP[hs$treatment==0] <- NA
hs <- filter(hs,treatment==0|!is.na(ncp))
                                        #hs$everCP[hs$treatment==1&is.na(hs$ncp)] <- FALSE

## model cp?
library(lme4)
mod1 <- glmer(everCP~state+grade+race+sex+frl+xirt+esl+(1|schoolid2),data=hs,subset=year==1,family=binomial)
mod2 <- glmer(everCP~state+grade+race+sex+frl+xirt+esl+(1|schoolid2),data=hs,subset=year==2,family=binomial)

## X is a design matrix without an intercept
## columns should have mean 0

hs1 <- filter(hs,year==1)
hs2 <- filter(hs,year==2)

form <- ~state+grade+race+sex+frl+xirt+esl
hs1 <- droplevels(hs1)
hs1 <- within(hs1,{
    S=ifelse(hs1$everCP,1,0)
    Z= treatment
    scl=schoolid2})

res1 <- est(hs1,covForm=form,clust=hs1$scl)
effs1 <- effsFromFit(res1)

hs2 <- droplevels(hs2)
hs2 <- within(hs2,{
    S=ifelse(hs2$everCP,1,0)
    Z= hs2$treatment
    scl=schoolid2})

res2 <- est(hs2,covForm=form,clust=hs2$scl)
effs2 <- effsFromFit(res2)

arm::binnedplot(res1$psMod$fitted,resid(res1$psMod,type='response'))
plot(res1$outMod)

arm::binnedplot(res2$psMod$fitted,resid(res2$psMod,type='response'))
plot(res2$outMod)


res1 <- mest(hs1)
est1 <- ests(res1,hs1)

res2 <- mest(hs2)
est2 <- ests(res2,hs2)

stan1 <- bayes(hs1)
stan2 <- bayes(hs2)

save(res1,est1,res2,est2,stan1,stan2,hs1,hs2,mest,bayes,file='ctResults.RData')


simCT1 <- function(newS,newY,outMod,data){
    data$S <- newS
    data <- within(data,{
                   S <- newS
                   Y <- newY+coef(outMod)['Sp']*S+coef(outMod)['Z:Sp']*Z*S
                   })
    cform <- formula(psModML)[-2]|>update(.~.-(1|scl))
    est(data,covForm=cform,clust=data$scl)|>effsFromFit()
}

simCT <- function(nsim,data,outMod,re.form=~0){
    psModML <- glmer(
        S~state + grade + race + sex + frl + xirt + esl + (1 | scl),
        family=binomial,
        data=subset(data,Z==1))
    outModML <- lmer(
        Y ~ state + grade + race + sex + frl + xirt + esl + Z + (1 | scl),
        data=data)

    newSs <- simulate(psModML,nsim=nsim,re.form=re.form,newdata=data,allow.new.levels=TRUE)
    newYs <- simulate(outModML,nsim=nsim,re.form=re.form)


    print(system.time(out <- mapply(simCT1,newSs,newYs,MoreArgs=list(outMod=outMod,data=data),SIMPLIFY='array')))

    out
}

res1Sim <- simCT(5,data=hs1,outMod=res1$outMod)
