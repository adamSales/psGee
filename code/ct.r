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


res1 <- mest(hs1)
est1 <- ests(res1,hs1)

res2 <- mest(hs2)
est2 <- ests(res2,hs2)

stan1 <- bayes(hs1)
stan2 <- bayes(hs2)

save(res1,est1,res2,est2,stan1,stan2,hs1,hs2,mest,bayes,file='ctResults.RData')
