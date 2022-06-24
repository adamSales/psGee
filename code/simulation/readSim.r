library(tidyverse)
library(parallel)

source('code/simulation/readSimFuncs.r')

#### read, process results from main simulation
print(load('simResults/pswResults.RData'))

results=loadRes(pswResults=pswResults)

resultsN100=loadRes(ext2='n100')

resultsNs=loadRes(ext2='ns')
