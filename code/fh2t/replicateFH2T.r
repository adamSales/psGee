library(tidyverse)
library(arm)
library(randomForest)
library(missForest)
library(estimatr)
library(splines)
library(rstan)
library(kableExtra)
library(texreg)
library(tableone)
library(xtable)
select <- dplyr::select


source('code/regression.r')
source('code/fh2t/psw.r')

auc <- function(x,y)
  if(length(unique(y))==2 & length(y)==length(x)){
    wilcox.test(x~y)$statistic/(sum(y)*(sum(1-y)))
  }else
    wilcox.test(x,y)$statistic/(legnth(x)*length(y))

aucMod <- function(mod)
  auc(mod$linear,1-mod$y)

############################################################
### make dataset
############################################################


if(file.exists('data/psdat.RData')){
  load('data/psdat.RData')
} else{
  source('code/fh2t/data.r')
}

## not enough variation in treatement assignment in school 7
psdat <- droplevels(filter(psdat,SchIDPre!=7))


psdat$Z <- ifelse(psdat$trt=='ASSISTments',1,0)

alts <- unique(psdat$trt[psdat$Z==0])
alts <- setNames(alts,alts)

psdat$S[psdat$Z==0]<- NA

### function to get data subset for 2-way comparison
getDat <- function(alt){
  dat <- psdat%>%
    filter(trt%in%c('ASSISTments',alt))
 # dat$pred <- predict(remnantMods[[alt]],dat)

  droplevels(dat)
}

############################################################
### estimate principal effects
############################################################

source('code/fh2t/estimates.r')

############################################################
### make plots
############################################################

source('code/fh2t/plots.r')

############################################################
### make tables (for appendix)
############################################################

source('code/fh2t/tables.r')
