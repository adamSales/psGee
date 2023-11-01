library(rstan)
library(dplyr)

alts <- c('BAU','FH2T','Dragon')
names(alts) <- alts

ests <- data.frame(eff=rep(paste0('eff',c(0,1)),3),
             estimates=NA,
             SE=NA,
             method='Mixture',
             Alternative=rep(alts,each=2))

for(alt in alts){
  mod=load(paste0('psStan',alt,'.RData'))
  mod=get(mod)
  s=summary(mod,par=c('bottomOuterATE','notbottomOuterATE'))$summary
  ests$estimates[ests$Alternative==alt&ests$eff=='eff0'] <- 
    s['notbottomOuterATE','mean']
  ests$estimates[ests$Alternative==alt&ests$eff=='eff1'] <- 
    s['bottomOuterATE','mean']
  ests$SE[ests$Alternative==alt&ests$eff=='eff0'] <- 
    s['notbottomOuterATE','sd']
  ests$SE[ests$Alternative==alt&ests$eff=='eff1'] <- 
    s['bottomOuterATE','sd']
}
  
saveRDS(ests,file='stanEsts.rds')
