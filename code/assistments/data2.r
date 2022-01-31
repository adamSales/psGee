library(tidyverse)
library(fixest)
library(sandwich)

dat <- read_csv('assistmentsData/principal_stratification_dataset.csv')

tsidVid <- unique(dat$assigned_tsid[dat$contains_video])
tsidTxt <- unique(dat$assigned_tsid[!dat$contains_video])

rands <- dat%>%select(contains('tsid'))%>%
    distinct()%>%
    rowwise()%>%
    mutate(alts=paste(sort(c_across()),collapse=';'))

vid <- strsplit(rands$alts,';')%>%map(~c(nAlts=length(.),vids=sum(.%in%tsidVid)))%>%reduce(rbind)

rands <- cbind(rands,vid)

dat <- left_join(dat,rands)

save(dat,rands,file='assistmentsData/dataWcontrasts.RData')

dat <- mutate(dat,perVids=vids/nAlts)

byExp <- dat%>%#filter(vids>0,nstud>800)%>%
  group_by(alts)%>%summarize(n=n(),pVidstud=mean(contains_video),perVids=mean(perVids),diff=pVidstud-perVids)


byExp%>%
  ggplot(aes(n,pVidstud,col=perVids,yintercept=perVids))+geom_point()+geom_hline(aes(yintercept=perVids))

dat$npc=ifelse(is.na(dat$next_problem_correct),0,dat$next_problem_correct)
mod=feols(npc~contains_video|alts,data=dat,vcov='hetero')


dat%>%group_by(alts)%>%
  mutate(nstud=n())%>%
  filter(nstud>100)%>%
  summarize(nstud=nstud[1],eff=mean(npc[contains_video])-mean(npc[!contains_video]))%>%
  ggplot(aes(nstud,eff))+geom_point()+geom_smooth(method='lm')

dat <- dat%>%
  filter(vids>0)%>%
  group_by(alts)%>%
  mutate(
    Zadj=contains_video-mean(contains_video),
    Yadj=npc-mean(npc)
    )

mod2 <- lm(Yadj~Zadj,data=dat)

dat$S <- dat$time_on_task>dat$video_length

missMod <- feols(is.na(time_on_task)~contains_video|alts,data=dat,vcov = 'hetero')

missMod2 <- feols(is.na(time_on_task)~contains_video|alts,
                  data=subset(dat,!is.na(student_prior_median_time_on_task)),
                  vcov = 'hetero')

dat <- mutate(dat,S1=ifelse(contains_video&is.na(time_on_task),FALSE,S))

dat <- filter(dat,!is.na(student_prior_median_time_on_task))

save(dat,rands,file='assistmentsData/dataWcontrasts.RData')

dat <- dat%>%mutate(nstud=n(),nStudVid=sum(contains_video),pCorrect=mean(npc))%>%filter(nstud>=100,nStudVid>10,pCorrect>0,pCorrect<1)


