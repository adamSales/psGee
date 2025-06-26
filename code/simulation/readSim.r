library(tidyverse)
library(parallel)

source('code/simulation/readSimFuncs.r')

#### read, process results from main simulation
print(load('simResults/pswResults.RData'))

results=loadRes(pswResults=pswResults)
save(results,file='simResults/fullResults.RData')

resultsN100=loadRes(ext2='n100_mu01is0')
#save(resultsN100,file='simResults/resultsN100.RData')


resultsNs=loadRes(ext2="ns_mu01is0")

resultsB1s=loadRes(ext2='b1ss_mu01is0')

### merge results by n
load('simResults/casesns_mu01is0.RData')
casesNs=cases
load('simResults/casesn100.RData')
casesN100=cases
load('simResults/cases.RData')

resultsNs=bind_rows(
  resultsNs,
  filter(resultsN100,mu01%in%resultsNs$mu01[1],
         errDist%in%resultsNs$errDist,
         b1%in%resultsNs$b1,
         intS%in%resultsNs$intS,
         intZ%in%resultsNs$intZ
         ),
  filter(results,mu01==resultsNs$mu01[1],
         errDist%in%resultsNs$errDist,
         b1%in%resultsNs$b1,
         intS%in%resultsNs$intS,
         intZ%in%resultsNs$intZ
         ))%>%
    filter(estimator!="psw")

save(resultsNs,file='simResults/resultsNs_mu01is0.RData')

### merge results by B1s
casesTot=cases
load('simResults/casesb1ss_mu01is0.RData')
casesB1=cases

resultsB1s=bind_rows(
  resultsB1s,
  filter(results,run%in%(left_join(casesB1[1,],casesTot%>%mutate(run=1:n())%>%select(-b1))%>%pull(run)))
  )
save(resultsB1s,file='simResults/resultsB1s_mu01is0.RData')
