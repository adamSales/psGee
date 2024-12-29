library(tidyverse)
library(parallel)

source('code/simulation/readSimFuncs.r')

#### read, process results from main simulation
print(load('simResults/pswResults.RData'))

results=loadRes(pswResults=pswResults)
save(results,file='simResults/fullResults.RData')

resultsN100=loadRes(ext2='n100')
save(resultsN100,file='simResults/resultsN100.RData')


resultsNs=loadRes(ext2='ns')

resultsB1s=loadRes(ext2='b1s')

### merge results by n
load('simResults/casesns.RData')
casesNs=cases
load('simResults/casesn100.RData')
casesN100=cases
load('simResults/cases.RData')

resultsNs=bind_rows(
  resultsNs,
  filter(resultsN100,run%in%which(apply(casesN100[,-1],1,function(x) all(x==casesNs[1,-1])))),
  filter(results,run%in%which(apply(cases[,-1],1,function(x) all(x==casesNs[1,-1])))))

save(resultsNs,file='simResults/resultsNs.RData')

### merge results by B1s
casesTot=cases
load('simResults/casesb1s.RData')
casesB1=cases

resultsB1s=bind_rows(
  resultsB1s,
  filter(results,run%in%(left_join(casesB1[1,],casesTot%>%mutate(run=1:n())%>%select(-b1))%>%pull(run)))
  )
save(resultsB1s,file='simResults/resultsB1s')
