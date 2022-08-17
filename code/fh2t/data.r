library(tidyverse)
library(arm)
library(missForest)
select <- dplyr::select

Scale <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

actDat0=read_csv('data/logs_for_adam.csv')

dat=read_csv('data/DATA20220202_3591.csv',na=c('','NA','#NULL!'))

probMeta=read_csv('data/ASSISTments-Table 1.csv')
probMeta=select(probMeta,!starts_with('...'))
names(probMeta)=gsub(' ','',names(probMeta))
probMeta$Hints[probMeta$Note=='3']=3

bottomOut=
  probMeta%>%
  filter(Hints!='n/a',CorrectAnswer!='n/a')%>%
  mutate(Hints=as.numeric(Hints),
         Subproblem=ifelse(Subproblem%in%c('n/a','N/A'),'A',Subproblem))%>%
  select(ProblemID,ProblemSet,ProblemNumber,ProblemOrder,Hints,Subproblem)%>%
  inner_join(mutate(actDat0,Subproblem=LETTERS[problem_part]),
             by=c("ProblemID"="problem_id","Subproblem"))%>%
  select(-'...1')%>%
  group_by(StuID)%>%
  summarize(nbo=sum(action=='answer_hint',na.rm=TRUE),nprob=n_distinct(ProblemID))

dat <- left_join(dat,bottomOut)%>%
  filter(!is.na(post.total_math_score))

cat(sum(is.na(dat$nbo)&dat$rdm_condition=='ASSISTments'),' students with NA for # bottom out hints; setting to 0\n')

dat$nbo[is.na(dat$nbo)&dat$rdm_condition=='ASSISTments'] <- 0

med <- median(dat$nbo[dat$rdm_condition=='ASSISTments'],na.rm=TRUE)

dat$S <- dat$nbo>med

psdat <- dat%>%
  select(StuID:Performance.Level5,EIP:PercentInAttendance6,starts_with("pre"),Y=post.total_math_score,S)%>%
  mutate(
    raceEth=raceEthnicityFed%>%
      factor()%>%
      fct_lump_min(200)%>%
      fct_recode(`Hispanic/Latino`="1",Asian="3",White="6")%>%
      fct_relevel('White'),
    Gender=as.factor(Gender))

imp <- psdat%>%
  select(virtual, Gender,raceEth,Performance.Level5,EIP,Scale.Score5:PercentInAttendance6,starts_with("pre"))%>%
  missForest()

names(imp$ximp) <- paste0(names(imp$xmiss)

save(dat,file='psdat.RData')
