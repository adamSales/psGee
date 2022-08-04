library(tidyverse)
library(rstan)
load('data/psdat.RData')

psdat=psdat%>%
    group_by(class)%>%
    mutate(nz=n_distinct(Z))%>%
    ungroup()%>%
    filter(nz>1)


rstan_options(auto_write = TRUE)
options(mc.cores = 4)
gc()

formU <-    ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP

formY  <- update(formU,~.-virtualIMP+as.factor(class))

makeCovMat=function(form,data){
    X=model.matrix(form,data=data)[,-1]
    dummy=apply(X,2,function(x) all(x%in%c(0,1)))

    Xdum=X[,dummy]
    Xcon=X[,!dummy]
    Xcon=scale(Xcon)
    out=cbind(Xcon,Xdum)
    attr(out,'scale')=attributes(Xcon)$`scaled:scale`
    attr(out,'center')=attributes(Xcon)$`scaled:center`
    out
}


fakeDat=
    psdat%>%
    filter(Z==1)%>%
    mutate(Z=0)%>%
    bind_rows(filter(psdat,Z==1))



                                        #checkdat <- read.csv("data/FakeDataForCheck.csv")
### Data Prep ####




XU=makeCovMat(formU,data=fakeDat)
XY=makeCovMat(formU,data=fakeDat)

XctlU=XU[fakeDat$Z==0,]
XtrtU=XU[fakeDat$Z==1,]
XctlY=XY[fakeDat$Z==0,]
XtrtY=XY[fakeDat$Z==1,]



stanDat <-  list(nc= sum(1-fakeDat$Z),#length(fakeDat[fakeDat$Z == 0,]$student_number), #
                 nt= sum(fakeDat$Z),#length(psdat[psdat$Z == 1,]$student_number), #
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √

                 YctlY=  (fakeDat[fakeDat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                #    YctlU= , this doesn't exist, right?

                 YtrtY= (fakeDat[fakeDat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 #   YtrtU=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),

                 XctlU=XctlU, # √
                 XctlY=XctlY, # √

                 XtrtU=XtrtU, # √
                 XtrtY=XtrtY, # √

                 bottomOuter=(ifelse(fakeDat[fakeDat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
                 )




checkmod2 <- stan('code/fh2t/psMod.stan',
                 data=stanDat,
                 iter = 4000)

save(checkmod2,stanDat,XU,XY,file='checkModNoClass.RData')


fakeDatT=filter(psdat,Z==1)

clsDat=function(cls){
    nc=sum(psdat$class==cls&psdat$Z==0)
    trtDat=filter(psdat,Z==1,class==cls)
    trtDat$Z <- 0
    trtDat[sample.int(nrow(trtDat),nc,replace=TRUE),]
}

clsDat=function(cls){
    nc=sum(psdat$class==cls&psdat$Z==0)
    trtDat=filter(psdat,Z==1,class==cls)
    trtDat$Z <- 0
                                        #    trtDat[sample.int(nrow(trtDat),nc,replace=TRUE),]
    trtDat
}


fakeDatC <-     map_dfr(unique(psdat$class),clsDat)

fakeDat <- bind_rows(
    fakeDatT,
    fakeDatC)


XU=makeCovMat(formU,data=fakeDat)
XY=makeCovMat(formY,data=fakeDat)

XctlU=XU[fakeDat$Z==0,]
XtrtU=XU[fakeDat$Z==1,]
XctlY=XY[fakeDat$Z==0,]
XtrtY=XY[fakeDat$Z==1,]



stanDat <-  list(nc= sum(1-fakeDat$Z),#length(fakeDat[fakeDat$Z == 0,]$student_number), #
                 nt= sum(fakeDat$Z),#length(psdat[psdat$Z == 1,]$student_number), #
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √

                 YctlY=  (fakeDat[fakeDat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                #    YctlU= , this doesn't exist, right?

                 YtrtY= (fakeDat[fakeDat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 #   YtrtU=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),

                 XctlU=XctlU, # √
                 XctlY=XctlY, # √

                 XtrtU=XtrtU, # √
                 XtrtY=XtrtY, # √

                 bottomOuter=(ifelse(fakeDat[fakeDat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
                 )




checkmod1 <- stan('code/fh2t/psMod.stan',
                 data=stanDat,
                 iter = 4000)

save(checkmod1,stanDat,XU,XY,file='checkMod.RData')

fakeDatT=filter(psdat,Z==1)

clsDat=function(cls){
    nc=sum(psdat$class==cls&psdat$Z==0)
    trtDat=filter(psdat,Z==1,class==cls)
    trtDat$Z <- 0
    trtDat[sample.int(nrow(trtDat),nc,replace=TRUE),]
}



fakeDatC <-     map_dfr(unique(psdat$class),clsDat)

fakeDat <- bind_rows(
    fakeDatT,
    fakeDatC)


XU=makeCovMat(formU,data=fakeDat)
XY=makeCovMat(formY,data=fakeDat)

XctlU=XU[fakeDat$Z==0,]
XtrtU=XU[fakeDat$Z==1,]
XctlY=XY[fakeDat$Z==0,]
XtrtY=XY[fakeDat$Z==1,]



stanDat <-  list(nc= sum(1-fakeDat$Z),#length(fakeDat[fakeDat$Z == 0,]$student_number), #
                 nt= sum(fakeDat$Z),#length(psdat[psdat$Z == 1,]$student_number), #
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √

                 YctlY=  (fakeDat[fakeDat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                #    YctlU= , this doesn't exist, right?

                 YtrtY= (fakeDat[fakeDat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 #   YtrtU=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),

                 XctlU=XctlU, # √
                 XctlY=XctlY, # √

                 XtrtU=XtrtU, # √
                 XtrtY=XtrtY, # √

                 bottomOuter=(ifelse(fakeDat[fakeDat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
                 )




checkmod3 <- stan('code/fh2t/psMod.stan',
                 data=stanDat,
                 iter = 4000)

save(checkmod3,stanDat,XU,XY,file='checkModResamp.RData')


