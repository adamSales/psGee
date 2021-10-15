library(rstan)

print(load('ctResults.RData'))
##  [1] "res1"  "est1"  "res2"  "est2"  "stan1" "stan2" "hs1"   "hs2"   "mest"
# [10] "bayes"


pdat <- expand.grid(Year=c("Year 1","Year 2"),Rt=c('Not Reassigned','Reassigned'),Estimator=c('M-Est','Bayes'))

bayes1 <- summary(stan1,par=c('eff0','eff1'))$summary
bayes2 <- summary(stan2,par=c('eff0','eff1'))$summary

pdat$Estimate=c(
    est1['eff0','est'],
    est2['eff0','est'],
    est1['eff1','est'],
    est2['eff1','est'],
    bayes1['eff0','mean'],
    bayes2['eff0','mean'],
    bayes1['eff1','mean'],
    bayes2['eff1','mean'])

pdat$SE <-
c(
    est1['eff0','se'],
    est2['eff0','se'],
    est1['eff1','se'],
    est2['eff1','se'],
    bayes1['eff0','sd'],
    bayes2['eff0','sd'],
    bayes1['eff1','sd'],
    bayes2['eff1','sd'])

pdat$lo=pdat$Estimate-2*pdat$SE
pdat$hi=pdat$Estimate+2*pdat$SE

pdat%>%
    ggplot(aes(Estimator,Estimate,color=Estimator,ymin=lo,ymax=hi))+
    geom_point(size=3)+geom_errorbar(size=2,width=0)+
    geom_hline(yintercept=0,linetype='dotted')+
    facet_grid(Year~Rt)
ggsave('writeUps/ctResults.jpg')
