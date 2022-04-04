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

load('psdat.RData')
write.csv(psdat, "psdat.csv")
psdat<- read.csv("psdat.csv")
colnames(psdat)
table(psdat$Z,
      psdat$rdm_condition)

# U --> Model predicting Bottom Outers
# Y --> Model predicting Bottom Outcomes

### remove classes that don't have representation across both scripts
classes_by_z <- data.frame(table(psdat$class, psdat$Z)) %>%
  spread(Var2, Freq) %>%
  rename(class = Var1) %>%
  filter(
    `0` > 0 &
    `1` > 0) %>%
  dplyr::select(-`1`,
         -`0`)

psdat <- psdat %>%
  inner_join(
    classes_by_z
  )

length(unique(psdat$class))

table(psdat$raceIMP, psdat$Z)

### Data Prep ####
###### Dummy Codes ####
#  dummy codes for race
race_dummy <- fastDummies::dummy_cols(psdat$raceIMP, 
                                      remove_first_dummy = T) %>%
  distinct() %>%
  rename(raceIMP = ".data") 
names(dummy)[-1]<-sub(".data*", "", names(dummy)[-1])

#  dummy codes for class
class_dummy<-fastDummies::dummy_cols(psdat$class, 
                                     remove_first_dummy = T) %>%
  distinct() %>%
  rename(class = ".data") 
names(dummy)[-1]<-sub(".data*", "", names(dummy)[-1])

#### Covariates matrixes #######
# Control Cov for U model
XctlU <- psdat %>%
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

# Control Cov for Y model
XctlY <- psdat %>%
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
  dplyr::select(-raceIMP, -class)


# Treatment Cov for U model
XtrtU <- psdat %>%
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


# Treatment Cov for Y model
XtrtY <- psdat %>%
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

stanDat <-  list(nc= length(psdat[psdat$Z == 0,]$student_number), # 
                    nt= length(psdat[psdat$Z == 1,]$student_number), # 
                    ncovU= ncol(XctlU), # √
                    ncovY= ncol(XctlY), # √
                    
                    YctlY=  (psdat[psdat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                #    YctlU= , this doesn't exist, right?
                    
                    YtrtY= (psdat[psdat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                    YtrtU=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),
                
                    XctlU=as.matrix(XctlU), # √
                    XctlY=as.matrix(XctlY), # √
                    
                    XtrtU=as.matrix(XtrtU), # √
                    XtrtY=as.matrix(XtrtY), # √
                
                    bottomOuter=(ifelse(psdat[psdat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
                      )

md.pattern(stanDat$XctlU)  

### Model ####
mod <- stan('psMod.stan',data=stanDat, iter = 4000)
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
checkdat <- read.csv("FakeDataForCheck.csv")
md.pattern(checkdat)
### Data Prep ####


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
md.pattern(XctlU)

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
md.pattern(XctlY)

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
md.pattern(XtrtU)

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
md.pattern(XtrtY)
stanDat <-  list(nc= length(checkdat[checkdat$Z == 0,]$student_number), # 
                 nt= length(checkdat[checkdat$Z == 1,]$student_number), # 
                 ncovU= ncol(XctlU), # √
                 ncovY= ncol(XctlY), # √
                 
                 YctlY=  (checkdat[checkdat$Z == 0, ]$Y), ### IS THE OUTCOME Y OR post.total_math_score?
                 #    YctlU= , this doesn't exist, right?
                 
                 YtrtY= (checkdat[checkdat$Z == 1, ]$Y ), ### IS THE OUTCOME Y OR post.total_math_score?
                 YtrtU=(ifelse(checkdat[checkdat$Z == 1, ]$anyBottom == "TRUE", 1, 0)),
                 
                 XctlU=as.matrix(XctlU), # √
                 XctlY=as.matrix(XctlY), # √ 
                 
                 XtrtU=as.matrix(XtrtU), # √
                 XtrtY=as.matrix(XtrtY), # √
                 
                 bottomOuter=(ifelse(checkdat[checkdat$Z == 1, ]$anyBottom == "TRUE", 1, 0))
)

### Model ####
checkmod2 <- stan('psMod.stan',
                 data=stanDat,
                 iter = 4000,
                 control = list(max_treedepth = 10))
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
