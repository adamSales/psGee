library(tidyverse)
library(arm)
library(missForest)
source('code/regression.r')
select <- dplyr::select

Scale <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

dat <- read_csv('data/data.csv')
dat$hints <- dat%>%select(contains("hint_count"))%>%rowSums(na.rm=TRUE)
dat$anyHint <- dat$hints>0

dat$bottom <- dat%>%select(contains("bottom"))%>%rowSums(na.rm=TRUE)
dat$anyBottom <- dat$bottom>0


full <- read_csv('data/Assessment_merged_2022_01_19_N=4,343 - Sheet1.csv')
older <- read_csv('data/Assessment_merged_2021_07_16_state_assessment_N=4321 - Sheet1.csv')

full <- older%>%
  select(student_number,pre.avg_time_on_tasks,pre_MA_total_score,pre_negative_reaction_score,pre_numerical_confindence_score)%>%
  right_join(full)



full=full%>%
  mutate(
    race=raceEthnicity%>%
      factor()%>%
      fct_lump_min(300)%>%
      fct_recode(`Hispanic/Latino`="H",Asian="A",White="W")%>%
      fct_relevel('White'),
    pretest = pre.total_math_score-round(mean(pre.total_math_score,na.rm=TRUE))
  )

full <- dat%>%select(student_number,hints,bottom,anyHint,anyBottom)%>%
  right_join(full)

dat0 <- full

dat0 <- dat0%>%
  mutate(
    ScaleScore7 = Scale(Scale.Score7),
    hasPretest=is.finite(pre.total_math_score),
    hasPosttest=is.finite(post.total_math_score),
    Z =rdm_condition
  )%>%
  filter(!is.na(Z))

teachDrop <- read_csv('data/IES_school_teacher_ID_list_opt_out - teacher (1).csv')

dat0$schoolSupp <- dat0$initial_school_id
dat0$schoolSupp[dat0$initial_teacher_id%in%teachDrop$teacher_id[teachDrop$note=='S03']] <- 'S03'
dat0$schoolSupp[dat0$initial_teacher_id%in%teachDrop$teacher_id[teachDrop$note=='S07']] <- 'S07'

dat1 <- subset(dat0,!schoolSupp%in%c('S03','S07'))
dat2 <- dat1%>%
  filter(!endsWith(Z,'Resource'))
###dat3 <- filter(dat2,!is.na(ScaleScore7),hasPretest)
dat3 <- filter(dat2,!is.na(post.total_math_score))#,hasPretest)


### use final teacher/class id to impute missing initial teacher/class id
dat3 <- mutate(dat3,
               teach=ifelse(is.na(initial_teacher_id),final_teacher_id,initial_teacher_id),
               class=ifelse(is.na(initial_teacher_class),final_teacher_class,initial_teacher_class))

xtabs(~Z+anyBottom,addNA=TRUE,data=dat3)

### try modeling anyBottom for Z=ASSISTments
mod1 <- glm(anyBottom~
              pretest+Scale.Score5+MALE+race+#as.factor(class)+
              #as.factor(initial_school_id)+
              virtual+EIP+IEP+ESOL+GIFTED+AbsentDays6+MOBILE6,
            family=binomial,
            data=subset(dat3,Z=='ASSISTments'))

AUCmod(mod1)

summary(mod1)

### binnedplot(mod1$fitted.values,resid(mod1,type='response'))

###with(subset(dat3,Z=='ASSISTments'&!is.na(Scale.Score5)),
                ### binnedplot(Scale.Score5,anyBottom))

mod2 <- glm(anyBottom~
              pretest+poly(Scale.Score5,2)+MALE+race+
              virtual+EIP+IEP+ESOL+GIFTED+AbsentDays5+MOBILE5,
            family=binomial,
            data=subset(dat3,Z=='ASSISTments'&!is.na(Scale.Score5)))

AUCmod(mod2)

summary(mod2)

###binnedplot(mod2$fitted.values,resid(mod2,type='response'))

###anova(mod1,mod2,test='Chisq')

### OK let's do some imputaiton
dat3%>%
  filter(Z=='ASSISTments')%>%
  dplyr::select(pretest,Scale.Score5,MALE,race,virtual,EIP,IEP,ESOL,GIFTED,AbsentDays5,MOBILE5)%>%sapply(function(x) mean(is.na(x)))

mod3 <- glm(anyBottom~
              as.factor(pretest)+Scale.Score5+MALE+race+
              virtual+EIP+IEP+ESOL+GIFTED,
            family=binomial,
            data=subset(dat3,Z=='ASSISTments'))

AUCmod(mod3)

summary(mod3)

###binnedplot(mod3$fitted.values,resid(mod3,type='response'))
AIC(mod3)
AIC(mod1)

mod4 <- update(mod1,.~.+pre.avg_time_on_tasks+pre_MA_total_score+pre_negative_reaction_score+pre_numerical_confindence_score)

summary(mod4)
AUCmod(mod4)

###binnedplot(mod4$fitted.values,resid(mod4,type='response'))

### let's do some imputation
xmis <- dat3 %>%
  dplyr::select(pretest,
                Scale.Score5,
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
save(imp,xmis,file='data/imputation.RData')
imp$OOBerror
names(imp$ximp) <- paste0(names(imp$ximp), 'IMP')
ximp <- imp$ximp%>%
  mutate(across(where(~is.factor(.)&all(levels(.)%in%c('0','1'))),
                ~as.numeric(as.character(.))))

dat3 <- bind_cols(dat3, ximp)

mod4imp <- glm(
  anyBottom ~ pretestIMP+ Scale.Score5IMP+ MALEIMP+ raceIMP+ virtualIMP+ EIPIMP+
  IEPIMP+ ESOLIMP+ GIFTEDIMP+ pre.avg_time_on_tasksIMP+
  pre_MA_total_scoreIMP+ pre_negative_reaction_scoreIMP+ pre_numerical_confindence_scoreIMP,
  data=subset(dat3,Z=='ASSISTments'),family=binomial)

summary(mod4imp)
AUCmod(mod4imp)

###binnedplot(mod4imp$fitted.values,resid(mod4imp,type='response'))


### let's go with mod4imp!

psdat <- dat3%>%
  filter(Z%in%c('ASSISTments','FH2T'))%>%
  mutate(Z=ifelse(Z=='ASSISTments',1,0),
         Y=post.total_math_score,
         S=as.numeric(anyBottom))

save(psdat,dat3,file="data/psdat.RData")
##

psdatBAU <- dat3%>%
  filter(Z%in%c('ASSISTments','BAU'))%>%
  mutate(Z=ifelse(Z=='ASSISTments',1,0),
         Y=post.total_math_score,
         S=as.numeric(anyBottom))

save(psdatBAU,dat3,file="data/psdatBAU.RData")
##
