library(tidyverse)
library(rstan)
library(geex)

load('../data/problemLevelUsageData/probLevelData.RData')
load('../data/RANDstudyData/HSdata.RData')
hs <- dat
print(load('../data/RANDstudyData/HSdata.RData'))
ms <- dat

hs$xirt <- hs$xirt/ifelse(hs$year==1,sd(hs$xirt[hs$year==1]),sd(hs$xirt[hs$year==2]))

load('../data/sectionLevelUsageData/advanceDataWInfo.RData')

numProbs <- x%>%group_by(field_id,study.year,gradeLevel)%>%summarize(n=n(),nprob=n_distinct(Prob1))

hist(numProbs$n,100)
hist(numProbs$nprob,100)

numProbs%>%ggplot(aes(nprob,y=..density..))+geom_histogram(bins=100)+facet_grid(gradeLevel~study.year)

numProbs%>%
    filter(study.year<3,n<1000)%>%
    group_by(gradeLevel,study.year)%>%
    mutate(N=n())%>%
    group_by(gradeLevel,study.year,n)%>%
    summarize(dens=n()/N)%>%
    ggplot(aes(n,y=dens))+geom_point()+facet_grid(gradeLevel~study.year)

numSec <- x%>%group_by(field_id,study.year,gradeLevel)%>%summarize(nsec=n_distinct(unit,section))

hist(numSec$nsec,100)

numSec%>%ggplot(aes(nsec,y=..density..))+geom_histogram(bins=100)+facet_grid(gradeLevel~study.year)

numSec%>%
    filter(study.year<3)%>%
    group_by(gradeLevel,study.year)%>%
    mutate(N=n())%>%
    group_by(gradeLevel,study.year,nsec)%>%
    summarize(dens=n()/N)%>%
    ggplot(aes(nsec,y=dens))+geom_point()+facet_grid(gradeLevel~study.year)

totalTime <-
    x%>%filter(total_t1>0)%>%group_by(field_id,study.year,gradeLevel)%>%summarize(time=sum(total_t1,na.rm=TRUE))

hist(log(totalTime$time),100)


