library(tidyverse)
library(arm)
library(sandwich)
library(loop.estimator)
source('code/regression.r')
select <- dplyr::select

load('data/psdat.RData')

### posttest
dat3$z=ifelse(dat3$Z=='ASSISTments',1,0)
### contrast assistments with each of the other conditions
postMods=
    lapply(c('BAU','Dragon','FH2T')%>%setNames(.,.),function(alt)
           lm(  post.total_math_score ~ z + class+pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+  EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=dat3,subset=Z=='ASSISTments'|Z==alt))

postMods$pooledTrt=update(postMods$BAU,subset=Z!='BAU')

vcvs=lapply(postMods,vcovHC,type='HC1')

effs=sapply(postMods,function(x) coef(x)['z'])
ses=sqrt(sapply(vcvs,function(x) x['z','z']))
Ts=effs/ses
pvals=2*pnorm(-abs(Ts))

cbind(effs,ses,Ts,pvals)
