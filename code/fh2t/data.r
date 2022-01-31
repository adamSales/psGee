library(tidyverse)
library(arm)
library(missForest)
source('code/regression.r')
select <- dplyr::select

Scale <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

dat <- read_csv('../../fh2t/data/data.csv')
dat$hints <- dat%>%select(contains("hint_count"))%>%rowSums(na.rm=TRUE)
dat$anyHint <- dat$hints>0

dat$bottom <- dat%>%select(contains("bottom"))%>%rowSums(na.rm=TRUE)
dat$anyBottom <- dat$bottom>0


full <- read_csv('../../fh2t/data/Assessment_merged_2021_07_16_state_assessment_N=4321 - Sheet1.csv')%>%
  mutate(
    race=student_raceEthnicityFed%>%
      factor()%>%
      fct_lump_min(300)%>%
      fct_recode(`Hispanic/Latino`="1",Asian="3",White="6")%>%
      fct_relevel('White'),
    pretest = pre.total_math_score-round(mean(pre.total_math_score,na.rm=TRUE)),
    ScaleScore5=`ScaleScore_5th grade`
  )

full <- dat%>%select(student_number,hints,bottom,anyHint,anyBottom)%>%
  right_join(full)

dat0 <- full

dat0 <- dat0%>%
  mutate(
    ScaleScore7 = Scale(Scale.Score_7th.grade),
    hasPretest=is.finite(pre.total_math_score),
    hasPosttest=is.finite(post.total_math_score),
    Z =rdm_condition
  )%>%
  filter(!is.na(Z))

teachDrop <- read_csv('../../fh2t/data/IES_school_teacher_ID_list_opt_out - teacher (1).csv')

dat0$schoolSupp <- dat0$initial_school_id
dat0$schoolSupp[dat0$initial_teacher_id%in%teachDrop$teacher_id[teachDrop$note=='S03']] <- 'S03'
dat0$schoolSupp[dat0$initial_teacher_id%in%teachDrop$teacher_id[teachDrop$note=='S07']] <- 'S07'

dat1 <- subset(dat0,!schoolSupp%in%c('S03','S07'))
dat2 <- dat1%>%
  filter(!endsWith(Z,'Resource'))
dat3 <- filter(dat2,!is.na(ScaleScore7),hasPretest)


### use final teacher/class id to impute missing initial teacher/class id
dat3 <- mutate(dat3,
               teach=ifelse(is.na(initial_teacher_id),final_teacher_id,initial_teacher_id),
               class=ifelse(is.na(initial_teacher_class),final_teacher_class,initial_teacher_class))

xtabs(~Z+anyBottom,addNA=TRUE,data=dat3)

### try modeling anyBottom for Z=ASSISTments
mod1 <- glm(anyBottom~
              pretest+ScaleScore5+MALE+race+#as.factor(class)+
              #as.factor(initial_school_id)+
              virtual+EIP+IEP+ESOL+GIFTED+AbsentDays6+MOBILE6,
            family=binomial,
            data=subset(dat3,Z=='ASSISTments'))

AUCmod(mod1)

summary(mod1)

binnedplot(mod1$fitted.values,resid(mod1,type='response'))
           
with(subset(dat3,Z=='ASSISTments'&!is.na(ScaleScore5)),
                 binnedplot(ScaleScore5,anyBottom)) 

mod2 <- glm(anyBottom~
              pretest+poly(ScaleScore5,2)+MALE+race+
              virtual+EIP+IEP+ESOL+GIFTED+AbsentDays5+MOBILE5,
            family=binomial,
            data=subset(dat3,Z=='ASSISTments'&!is.na(ScaleScore5)))

AUCmod(mod2)

summary(mod2)

binnedplot(mod2$fitted.values,resid(mod2,type='response'))

anova(mod1,mod2,test='Chisq')

### OK let's do some imputaiton
dat3%>%
  filter(Z=='ASSISTments')%>%
  dplyr::select(pretest,ScaleScore5,MALE,race,virtual,EIP,IEP,ESOL,GIFTED,AbsentDays5,MOBILE5)%>%sapply(function(x) mean(is.na(x)))

mod3 <- glm(anyBottom~
              as.factor(pretest)+ScaleScore5+MALE+race+
              virtual+EIP+IEP+ESOL+GIFTED,
            family=binomial,
            data=subset(dat3,Z=='ASSISTments'))

AUCmod(mod3)

summary(mod3)

binnedplot(mod3$fitted.values,resid(mod3,type='response'))
AIC(mod3)
AIC(mod1)

mod4 <- update(mod1,.~.+pre.avg_time_on_tasks+pre_MA_total_score+pre_negative_reaction_score+pre_numerical_confindence_score)

summary(mod4)
AUCmod(mod4)

binnedplot(mod4$fitted.values,resid(mod4,type='response'))

### let's do some imputation
xmis <- dat3 %>%
  dplyr::select(pretest,
                ScaleScore5,
                MALE,
                race,
                virtual,
                EIP,
                IEP,
                ESOL,
                GIFTED,
                AbsentDays5,
                MOBILE5,
                pre.avg_time_on_tasks,          
                  pre_MA_total_score,             
                 pre_negative_reaction_score,    
                 pre_numerical_confindence_score,
                pre_MA_avg_score
                ) %>%
  mutate(across(where(~all(na.omit(.)%in%c(0,1),na.rm = TRUE)),as.factor))%>%
  as.data.frame()

imp <- missForest(xmis, variablewise = TRUE)
names(imp$OOBerror) <- names(imp$ximp)
imp$OOBerror
names(imp$ximp) <- paste0(names(imp$ximp), 'IMP')
ximp <- imp$ximp%>%
  mutate(across(where(~is.factor(.)&all(levels(.)%in%c('0','1'))),
                ~as.numeric(as.character(.))))

dat3 <- bind_cols(dat3, ximp)

mod4imp <- glm(
  anyBottom ~ pretestIMP+ ScaleScore5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+ 
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+ 
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(dat3,Z=='ASSISTments'),family=binomial)

summary(mod4imp)
AUCmod(mod4imp)

binnedplot(mod4imp$fitted.values,resid(mod4imp,type='response'))


### let's go with mod4imp!

psdat <- dat3%>%
  filter(Z%in%c('ASSISTments','FH2T'))%>%
  mutate(Z=ifelse(Z=='ASSISTments',1,0),
         Y=ScaleScore7,
         S=as.numeric(anyBottom))

estimate0 <- est(psdat,covFormU = formula(mod4imp)[-2])
                # covFormY=
                #   update(formula(mod4imp)[-2],
                #          .~.-virtualIMP),block = "class")
plot(estimate0$outMod)
effsFromFit(estimate0)

estimate1 <- est(psdat,covFormU = formula(mod4imp)[-2],
                covFormY=
                  update(formula(mod4imp)[-2],
                         .~.-virtualIMP),block = "class")
plot(estimate1$outMod)
effsFromFit(estimate1)


### test standard errors
bs0=replicate(1000,
              effsFromFit(
                est(
                  psdat[sample(1:nrow(psdat),nrow(psdat),replace=TRUE),],
                  covFormU = formula(mod4imp)[-2])))

bs1=replicate(1000,
              effsFromFit(
                est(
                  psdat[sample(1:nrow(psdat),nrow(psdat),replace=TRUE),],
                  covFormU = formula(mod4imp)[-2],
                   covFormY=update(formula(mod4imp)[-2],
                                   .~.-virtualIMP),block = "class")))
