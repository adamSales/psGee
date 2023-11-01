library(rstan)
library(dplyr)
options(mc.cores = 4)
rstan_options(auto_write=TRUE)

load('stanStuff.RData')

psStanBAU=stan('code/fh2t/psMod.stan',data=makeSdat('BAU'))
save(psStanBAU,file='psStanBAU.RData')

psStanFH2T=stan('code/fh2t/psMod.stan',data=makeSdat('FH2T'))
save(psStanFH2T,file='psStanFH2T.RData')

psStanDragon=stan('code/fh2t/psMod.stan',data=makeSdat('Dragon'))
save(psStanDragon,file='psStanDragon.RData')