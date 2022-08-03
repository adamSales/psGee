library(tidyverse)
library(dplyr)
library(arm)
library(modEvA)
library(fastDummies)
library(rstan)
library(magrittr)
library(mice)
#source('code/regression.r')


rstan_options(auto_write = TRUE)
options(mc.cores = 4)
gc()

load('data/psdat.RData')
## write.csv(psdat, "data/psdat.csv")
## psdat<- read.csv("data/psdat.csv")
## colnames(psdat)
## table(psdat$Z,
##       psdat$rdm_condition)

# U --> Model predicting Bottom Outers
# Y --> Model predicting Bottom Outcomes


xtabs(~class+Z,psdat)%>%as.data.frame()%>%summarize(sum(Freq==0))
length(unique(psdat$class))

psdat=psdat%>%
    group_by(class)%>%
    mutate(nz=n_distinct(Z))%>%
    ungroup()%>%
    filter(nz>1)

length(unique(psdat$class))

table(psdat$raceIMP, psdat$Z)

formU <-    ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP

formY  <- update(formU,~.-virtualIMP+as.factor(class))

makeCovMat=function(form,data){
    X=model.matrix(form,data=data)[,-1]
    sdX=apply(X,2,sd)
    meanX=colMeans(X)
    scale(X,center=meanX,scale=sdX)
}

XU=makeCovMat(formU,data=psdat)
XY=makeCovMat(formY,data=psdat)

XctlU=XU[psdat$Z==0,]
XtrtU=XU[psdat$Z==1,]
XctlY=XY[psdat$Z==0,]
XtrtY=XY[psdat$Z==1,]



stanDat <-  list(nc= sum(1-psdat$Z),#length(psdat[psdat$Z == 0,]$student_number), #
                 nt= sum(psdat$Z),#length(psdat[psdat$Z == 1,]$student_number), #
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √

                 YctlY=  (psdat[psdat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                #    YctlU= , this doesn't exist, right?

                 YtrtY= (psdat[psdat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 #   YtrtU=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),

                 XctlU=XctlU, # √
                 XctlY=XctlY, # √

                 XtrtU=XtrtU, # √
                 XtrtY=XtrtY, # √

                 bottomOuter=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
                 )

md.pattern(stanDat$XctlU)

### Model ####
mod <- stan('code/fh2t/psMod.stan',data=stanDat, iter = 4000)

save(mod,stanDat,XU,XY,file='fittedModsStanPS.RData')


#### try it without class FE
sdatNoClass=stanDat
sdatNoClass$XctlY=sdatNoClass$XctlU
sdatNoClass$XtrtY=sdatNoClass$XtrtU
sdatNoClass$ncovY=sdatNoClass$ncovU
modNoClass <- stan('code/fh2t/psMod.stan',data=sdatNoClass, iter = 4000)

save(modNoClass,sdatNoClass,file='fittedModsStanPSnoClass.RData')

# saveRDS(mod, "model.rds")
# mod <- readRDS("model.rds")

# Warning message:
#   Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess

print(
  mod,
  pars = c(
    'alphaTBO',
    'alphaTNBO',
    'alphaCBO',
    'alphaCNBO',
    'bottomOuterATE',
    'notbottomOuterATE',
    'ATEdiff'
  ),
  probs = c(0.025, 0.975)
)

traceplot(
  mod,
  inc_warmup = T,
  pars = c(
    'alphaTBO',
    'alphaTNBO',
    'alphaCBO',
    'alphaCNBO',
    'bottomOuterATE',
    'notbottomOuterATE',
    'ATEdiff'
  )
)

traceplot(
  mod,
  inc_warmup = F,
  pars = c(
    'alphaTBO',
    'alphaTNBO',
    'alphaCBO',
    'alphaCNBO',
    'bottomOuterATE',
    'notbottomOuterATE',
    'ATEdiff'
  )
)

### Model Checking #####
## drop control group (FH2T)
## sample with replacement treatment group (asstistement) a sample the size of the control (FH2T)
## call the new sample control
## run same model

#
# checkdat_treatment <- psdat %>%
#   filter(rdm_condition == "ASSISTments")
#
# set.seed(4)
# coltrol_student_id <- as.data.frame(sample(checkdat_treatment$student_number, size = 853, replace = T))
# colnames(coltrol_student_id) <- c("student_number")
# install.packages("sampling")
# require(sampling)
# table(psdat$class[psdat$Z==0])
#
# checkdat_control<- sampling::strata(
#   data = checkdat_treatment %>%
#     dplyr::arrange(class),
#   stratanames = "class",
#   size = table(psdat$class[psdat$Z==0]),
#   method = "srswr"
# )
#
# checkdat_control<- getdata(
#                            checkdat_treatment%>%
#                              dplyr::arrange(class),
#                            checkdat_control)
# setdiff(colnames(checkdat_treatment), colnames(checkdat_control))
# colnames(checkdat_treatment)
# colnames(checkdat_control)
# table(checkdat_control$Z)
# checkdat <-  checkdat_treatment %>%
#   bind_rows(checkdat_control  %>%
#              dplyr::select(-Z) %>%
#              mutate(Z = 0) %>%
#              select(-Prob, -Stratum, -S, -anyBottom)
#           )
# table(checkdat$Z)
# table(psdat$Z)
#   table(checkdat$class, checkdat$Z)
#
# write.csv(checkdat, "FakeDataForCheck.csv")
checkdat <- read.csv("data/FakeDataForCheck.csv")
#md.pattern(checkdat)
### Data Prep ####

XU=makeCovMat(formU,data=checkdat)
XY=makeCovMat(formY,data=checkdat)


XctlU=XU[checkdat$Z==0,]
XtrtU=XU[checkdat$Z==1,]
XctlY=XY[checkdat$Z==0,]
XtrtY=XY[checkdat$Z==1,]



stanDat <-  list(nc= sum(1-checkdat$Z),#length(checkdat[checkdat$Z == 0,]$student_number), #
                 nt= sum(checkdat$Z),#length(psdat[psdat$Z == 1,]$student_number), #

#### Covariates matrixes #######
# Control Cov for U model
XctlU <- checkdat %>%
  dplyr::filter(Z == 0) %>%
  dplyr::select(
    pretestIMP,
    Scale.Score5IMP,
    MALEIMP,
    raceIMP,
    virtualIMP,
    EIPIMP,
    IEPIMP,
    ESOLIMP,
    GIFTEDIMP,
    pre.avg_time_on_tasksIMP,
    pre_MA_total_scoreIMP,
    pre_negative_reaction_scoreIMP,
    pre_numerical_confindence_scoreIMP
  ) %>%
  mutate(
    pretestIMP = scale(pretestIMP),
    Scale.Score5IMP = scale(Scale.Score5IMP),
    pre.avg_time_on_tasksIMP = scale(pre.avg_time_on_tasksIMP),
    pre_MA_total_scoreIMP = scale(pre_MA_total_scoreIMP),
    pre_negative_reaction_scoreIMP = scale(pre_negative_reaction_scoreIMP),
    pre_numerical_confindence_scoreIMP = scale(pre_numerical_confindence_scoreIMP)
  ) %>%
  left_join(
    race_dummy,
    by = "raceIMP"
  ) %>%
  dplyr::select(-raceIMP)
colnames(XctlU)
##md.pattern(XctlU)

# Control Cov for Y model
XctlY <- checkdat %>%
  filter(Z == 0) %>%
  dplyr::select(
    pretestIMP,
    Scale.Score5IMP,
    MALEIMP,
    raceIMP,
    #virtualIMP,
    EIPIMP,
    IEPIMP,
    ESOLIMP,
    GIFTEDIMP,
    pre.avg_time_on_tasksIMP,
    pre_MA_total_scoreIMP,
    pre_negative_reaction_scoreIMP,
    pre_numerical_confindence_scoreIMP,
    class
  ) %>%
  mutate(
    pretestIMP= scale(pretestIMP),
    Scale.Score5IMP= scale(Scale.Score5IMP),
    pre.avg_time_on_tasksIMP= scale(pre.avg_time_on_tasksIMP),
    pre_MA_total_scoreIMP= scale(pre_MA_total_scoreIMP),
    pre_negative_reaction_scoreIMP= scale(pre_negative_reaction_scoreIMP),
    pre_numerical_confindence_scoreIMP= scale(pre_numerical_confindence_scoreIMP)
  ) %>%
  left_join(
    class_dummy,
    by = "class"
  ) %>%
  left_join(
    race_dummy,
    by = "raceIMP"
  ) %>%
  dplyr::select(-raceIMP, -class)
colnames(XctlY)
##md.pattern(XctlY)

# Treatment Cov for U model
XtrtU <- checkdat %>%
  filter(Z == 1) %>%
  dplyr::select(
    pretestIMP,
    Scale.Score5IMP,
    MALEIMP,
    raceIMP,
    virtualIMP,
    EIPIMP,
    IEPIMP,
    ESOLIMP,
    GIFTEDIMP,
    pre.avg_time_on_tasksIMP,
    pre_MA_total_scoreIMP,
    pre_negative_reaction_scoreIMP,
    pre_numerical_confindence_scoreIMP,
    #class
  ) %>%
  mutate(
    pretestIMP = scale(pretestIMP),
    Scale.Score5IMP = scale(Scale.Score5IMP),
    pre.avg_time_on_tasksIMP = scale(pre.avg_time_on_tasksIMP),
    pre_MA_total_scoreIMP = scale(pre_MA_total_scoreIMP),
    pre_negative_reaction_scoreIMP = scale(pre_negative_reaction_scoreIMP),
    pre_numerical_confindence_scoreIMP = scale(pre_numerical_confindence_scoreIMP)
  ) %>%
  left_join(
    race_dummy,
    by = "raceIMP"
  ) %>%
  dplyr::select(-raceIMP)
colnames(XtrtU)
##md.pattern(XtrtU)

# Treatment Cov for Y model
XtrtY <- checkdat %>%
  filter(Z == 1) %>%
  dplyr::select(
    pretestIMP,
    Scale.Score5IMP,
    MALEIMP,
    raceIMP,
    #virtualIMP,
    EIPIMP,
    IEPIMP,
    ESOLIMP,
    GIFTEDIMP,
    pre.avg_time_on_tasksIMP,
    pre_MA_total_scoreIMP,
    pre_negative_reaction_scoreIMP,
    pre_numerical_confindence_scoreIMP,
    class
  ) %>%
  mutate(
    pretestIMP = scale(pretestIMP),
    Scale.Score5IMP = scale(Scale.Score5IMP),
    pre.avg_time_on_tasksIMP = scale(pre.avg_time_on_tasksIMP),
    pre_MA_total_scoreIMP = scale(pre_MA_total_scoreIMP),
    pre_negative_reaction_scoreIMP = scale(pre_negative_reaction_scoreIMP),
    pre_numerical_confindence_scoreIMP = scale(pre_numerical_confindence_scoreIMP)
  ) %>%
  left_join(
    class_dummy,
    by = "class"
  ) %>%
  left_join(
    race_dummy,
    by = "raceIMP"
  ) %>%
  dplyr::select(-raceIMP
                , -class
  )
colnames(XtrtY)
#md.pattern(XtrtY)
stanDat <-  list(nc= length(checkdat[checkdat$Z == 0,]$student_number), #
                 nt= length(checkdat[checkdat$Z == 1,]$student_number), #
>>>>>>> Stashed changes
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √

                 YctlY=  (checkdat[checkdat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                #    YctlU= , this doesn't exist, right?

                 YtrtY= (checkdat[checkdat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 #   YtrtU=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),

                 XctlU=XctlU, # √
                 XctlY=XctlY, # √

                 XtrtU=XtrtU, # √
                 XtrtY=XtrtY, # √

                 bottomOuter=(ifelse(checkdat[checkdat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
                 )




checkmod2 <- stan('code/fh2t/psMod.stan',
                 data=stanDat,
                 iter = 4000)

,
                 control = list(max_treedepth = 10))


sdatNoClass=stanDat
sdatNoClass$XtrtY=sdatNoClass$XtrtY[,-c(12:145)]
sdatNoClass$XctlY=sdatNoClass$XctlY[,-c(12:145)]
sdatNoClass$ncovY=ncol(sdatNoClass$XtrtY)

## justOutTrt <- stan('code/fh2t/psModJustOutcomeTrt.stan',data=sdatNoClass)
## save(justOutTrt,file='justOutTrt.RData')

sdatClass=sdatNoClass
checkdat$classN=as.numeric(droplevels(as.factor(checkdat$class)))
sdatClass$classT=checkdat$classN[checkdat$Z==1]
sdatClass$classC=checkdat$classN[checkdat$Z==0]
sdatClass$nclass=max(sdatClass$classT)
stopifnot(max(sdatClass$classT)==max(sdatClass$classC))
stopifnot(max(sdatClass$classT)==n_distinct(sdatClass$classT))
stopifnot(max(sdatClass$classT)==n_distinct(sdatClass$classC))


checkmodRandClass=stan('code/fh2t/psModRandClass.stan',data=sdatClass)
save(checkmodRandClass,file='checkmodRandClass.RData')


checkdathalf1=subset(checkdat,classN<67)
stanDatHalf <-  list(nc= length(checkdathalf1[checkdathalf1$Z == 0,]$student_number), #
                 nt= length(checkdathalf1[checkdathalf1$Z == 1,]$student_number), #
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √
                 YctlY=  (checkdathalf1[checkdathalf1$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                 #    YctlU= , this doesn't exist, right?
                 YtrtY= (checkdathalf1[checkdathalf1$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 YtrtU=(ifelse(checkdathalf1[checkdathalf1$Z == 1, ]$anyBottom == "TRUE", 1, 0)),
                 XctlU=as.matrix(XctlU)[checkdat$classN[checkdat$Z==0]<67,], # √
                 XctlY=as.matrix(XctlY)[checkdat$classN[checkdat$Z==0]<67,], # √
                 XtrtU=as.matrix(XtrtU)[checkdat$classN[checkdat$Z==1]<67,], # √
                 XtrtY=as.matrix(XtrtY)[checkdat$classN[checkdat$Z==1]<67,], # √
                 classT=checkdathalf1$classN[checkdathalf1$Z==1],
                 classC=checkdathalf1$classN[checkdathalf1$Z==0],
                 nclass=66,
                 bottomOuter=(ifelse(checkdathalf1[checkdathalf1$Z == 1, ]$anyBottom == "TRUE", 1, 0))
)

modHalf1=stan('code/fh2t/psModRandClass.stan',data=stanDatHalf)
save(modHalf1,file='modHalf1.RData')


                                        # saveRDS(checkmod2, "checkmod2.rds")
# mod <- readRDS("model.rds")
checkmod2
print(
  checkmod2,
  pars = c(
    'alphaTBO',
    'alphaTNBO',
    'alphaCBO',
    'alphaCNBO',
    'bottomOuterATE',
    'notbottomOuterATE',
    'ATEdiff'
  ),
  probs = c(0.025, 0.975)
)


traceplot(
  checkmod2,
  inc_warmup = F,
  pars = c(
    'alphaTBO',
    'alphaTNBO',
    'alphaCBO',
    'alphaCNBO',
    'bottomOuterATE',
    'notbottomOuterATE',
    'ATEdiff'
  )
)


names(stanDat)

summary(with(stanDat, glm(YtrtU~XtrtU, family = binomial())))
summary(with(stanDat, glm(YtrtY~XtrtY)))
