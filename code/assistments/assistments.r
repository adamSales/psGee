library(tidyverse)
library(estimatr)
library(geex)
library(rstan)

setwd(file.path(rprojroot::find_rstudio_root_file(),'CT','MOM'))

source('code/regression.r')

source('code/assistments/assistmentsData.r')

source('code/assistmentsGEE.r')

source('code/assistmentsStan.r')
