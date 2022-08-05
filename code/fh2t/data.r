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

dat <- dat%>%
  select(StuID:Performance.Level5,EIP:PercentInAttendance6,starts_with("pre"),Y=post.total_math_score,S)

save(dat,file='psdat.RData')
