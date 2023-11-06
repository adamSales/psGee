

#######################################################################
### principal score models. For record of model selection, see "psMods.r"
#######################################################################

### optimize AIC & check model fit
psMod7 <- glm(
  S ~ SchIDPre + Scale.Score5 + ESOL + IEP + pre.total_time_on_tasks +
    pre_MSE_total_score + fullYear5 + pre_MA_total_scoreNA +
    Gender + raceEth + GIFTED +
    poly(pre_PS_tasks_total_score,2, raw = TRUE),
  data=psdat%>%filter(trt=='ASSISTments')%>%select(-trt),
  family=binomial)

### use all available covariates
psModAll <- glm(
  S ~ SchIDPre + Gender + raceEth + Performance.Level5 + EIP +
    Scale.Score5 + ESOL + GIFTED + IEP + IST + SECTION504 + SST +
    AbsentDays5 + UnexcusedDays5 + AbsentDays6 + UnexcusedDays6 +
    pre.total_math_score + pre.sub_P_score + pre.sub_F_score +
    pre.math_completed_num + pre.total_time_on_tasks + pre_MA_total_score +
    pre_MSE_total_score + pre_PS_part2E_score + pre_PS_part2NE_score +
    pre_PS_completed_num + pre_PS_total_RT_sec + fullYear5 +
    fullYear6 + noUnexcused5 + Performance.Level5NA + AbsentDays6NA +
    pre.total_math_scoreNA + pre_MA_total_scoreNA + poly(pre_PS_tasks_total_score,
                                                         2, raw = TRUE),
    data=psdat%>%filter(trt=='ASSISTments')%>%select(-trt),
  family=binomial)

### use covariates other than prior math ability

psModRest <- update(
  psModAll,
  .~.-pre.total_math_score-Scale.Score5-Performance.Level5-pre.sub_P_score-pre.sub_F_score-pre.math_completed_num)




#######################################################################
### estimate ATEs and principal effects with GEEPERs
### for record of model selection, testing, see analysisGEE.r
#######################################################################

### ATEs, using all covariates
atesAll <-  sapply(setNames(alts,alts),
                  function(alt)
                    estimatr::lm_robust(
                                update(formula(psModAll),Y~Z+.-SchIDPre),
                            data=getDat(alt),fixed_effects=~ClaIDPre),
               simplify=FALSE)


### principal effects using psMod7
estimates3 <- lapply(setNames(alts,alts),
                    function(alt)
                      est(getDat(alt),
                          covFormU=formula(psMod7)[-2],
                          covFormY=~Scale.Score5+ESOL+IEP+pre.total_time_on_tasks+pre_MSE_total_score+
                            fullYear5+pre_PS_tasks_total_score+raceEth+Gender+GIFTED+
                            pre.total_math_score+pre.total_math_scoreNA+ClaIDPre
                                          ))

### principal effects using psModAll
estimates1all <-
  lapply(alts,
         function(alt)
           est(getDat(alt),psMod=psModAll,
               covFormY=update(formula(psModAll)[-2],
                               .~.
                               -SchIDPre+ClaIDPre
                               -UnexcusedDays6+ns(UnexcusedDays6,3)
                               -pre_MA_total_score+ns(pre_MA_total_score,3)
                               -pre_MSE_total_score+I(pre_MSE_total_score< -2)+ns(pre_MSE_total_score,3)
                               -pre.total_math_score+ns(pre.total_math_score,4))
               )
         )


### principal effects using psModRest
estimates1rest <-
  lapply(alts,
         function(alt)
           est(getDat(alt),psMod=psModRest,
               covFormY=update(formula(psModAll)[-2],
                               .~.
                               -SchIDPre+ClaIDPre
                               -UnexcusedDays6+ns(UnexcusedDays6,3)
                               -pre_MA_total_score+ns(pre_MA_total_score,3)
                               -pre_MSE_total_score+I(pre_MSE_total_score< -2)+ns(pre_MSE_total_score,3)
                               -pre.total_math_score+ns(pre.total_math_score,4)
                               -Scale.Score5+ns(Scale.Score5,3)
                               )
               )
         )

save(atesAll,estimates3,estimates1all,estimates1rest,file='results/geeResults.RData')




#######################################################################
### estimate principal effects with PSW
#######################################################################

### use psModAll
if(file.exists('results/psw.RData')){
  load('results/psw.RData')
} else{
  pswResults <- lapply(setNames(alts,alts),psw,psMod=psModAll)
  save(pswResults,file='results/psw.RData')
}

#######################################################################
### estimate principal effects with Bayesian Mixture Modeling
#######################################################################

if(file.exists('results/stanEsts.rds')){
  stanResults <- readRDS('results/stanEsts.rds')
} else {
  makeSdat <- function(alt,ests=estimates1all){
    dat <- getDat(alt)

    Xout <- model.matrix(ests[[alt]]$outMod)
    Xout <- Xout[,-c(1,which(colnames(Xout)%in%c('Z','Sp','Z:Sp')))]
    sdat <- list(
      YctlY=dat$Y[dat$Z==0],
      YtrtY=dat$Y[dat$Z==1],
      XctlU=model.matrix(formula(ests[[alt]]$psMod)[-2],data=subset(dat,Z==0))[,-1],
      XtrtU=model.matrix(ests[[alt]]$psMod)[,-1],
      XctlY=Xout[dat$Z==0,],
      XtrtY=Xout[dat$Z==1,],
      bottomOuter=dat$S[dat$Z==1],
      nc=sum(1-dat$Z),
      nt=sum(dat$Z)
    )
    sdat$ncovU=ncol(sdat$XctlU)
    sdat$ncovY=ncol(Xout)

    sdat
  }

  psStanBAU=stan('code/fh2t/psMod.stan',data=makeSdat('BAU'))
  save(psStanBAU,file='results/psStanBAU.RData')

  psStanFH2T=stan('code/fh2t/psMod.stan',data=makeSdat('FH2T'))
  save(psStanFH2T,file='results/psStanFH2T.RData')

  psStanDragon=stan('code/fh2t/psMod.stan',data=makeSdat('Dragon'))
  save(psStanDragon,file='results/psStanDragon.RData')

  stanResults <- data.frame(eff=rep(paste0('eff',c(0,1)),3),
             estimates=NA,
             SE=NA,
             method='Mixture',
             Alternative=rep(alts,each=2))

  for(alt in alts){
    mod=load(paste0('psStan',alt,'.RData'))
    mod=get(mod)
    s=summary(mod,par=c('bottomOuterATE','notbottomOuterATE'))$summary
    stanResults$estimates[stanResults$Alternative==alt&stanResults$eff=='eff0'] <-
      s['notbottomOuterATE','mean']
    stanResults$estimates[stanResults$Alternative==alt&stanResults$eff=='eff1'] <-
      s['bottomOuterATE','mean']
    stanResults$SE[stanResults$Alternative==alt&stanResults$eff=='eff0'] <-
      s['notbottomOuterATE','sd']
    stanResults$SE[stanResults$Alternative==alt&stanResults$eff=='eff1'] <-
      s['bottomOuterATE','sd']
  }

  saveRDS(stanResults,file='results/stanEsts.rds')
}
