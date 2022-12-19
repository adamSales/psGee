library(tidyverse)
library(arm)
library(randomForest)
library(tableone)
library(kableExtra)
library(xtable)
source('code/regression.r')
library(huxtable)
library(texreg)
select <- dplyr::select

options(knitr.kable.NA = '')

load('data/psdatBAU.RData')
load('data/psdat.RData')
load('data/imputation.RData')

tab1dat <-
  dat3%>%
  filter(Z%in%c('ASSISTments','BAU'))%>%
  select(
    pretest , Scale.Score5 , MALE , race ,
    virtual , EIP , IEP , ESOL , GIFTED , pre.avg_time_on_tasks ,
    pre_MA_total_score , pre_negative_reaction_score ,
    pre_numerical_confindence_score,anyBottom,post.total_math_score,Z)%>%
  mutate(pre.avg_time_on_tasks=log(pre.avg_time_on_tasks),
         MALE=ifelse(MALE==1,'Male','Female'),
         virtual=ifelse(virtual==1,'Remote','In-Person'),
         anyBottom=anyBottom*1)


oob=as.data.frame(rbind(imp$OOBerror))%>%  select(
    `Has EIP`=EIP,
    `Has IEP`=IEP,
    ESOL,
    `Gifted`=GIFTED,
        Pretest=pretest,
    `Grade 5 Stand. Test`=Scale.Score5,
    Gender=MALE,
    `Race/Ethnicity`=race,
    `log(Pretest ToT)`=pre.avg_time_on_tasks,
    `Math Anxiety`=pre_MA_total_score,
    `Math Neg. React.`=pre_negative_reaction_score,
    `Numerical Conf.`=pre_numerical_confindence_score)%>%t()



tab1bin <- tab1dat%>%
  select(
    `Has EIP (n (%))`=EIP,
    `Has IEP (n (%))`=IEP,
    `ESOL (n (%))`=ESOL,
    `Gifted (n (%))`=GIFTED,
    `Bottom Out (n (%))`=anyBottom,Z)%>%
  mutate(across(.fns=as.factor))%>%
  CreateTableOne(vars=setdiff(names(.),'Z'),data=.,strata="Z")%>%
  print(missing=TRUE,test=FALSE,dropEqual=TRUE,explain=FALSE)%>%
  cbind(level="",.,
        `Imputation Error`=round(oob[pmatch(gsub(" (n (%))","",rownames(.),fixed=TRUE),rownames(oob)),],2))




tab1oth<- tab1dat%>%
  select(
    Pretest=pretest,
    `Grade 5 Stand. Test`=Scale.Score5,
    `Gender  (n (%))`=MALE,
    `Race/Eth. (n (%))`=race,
    `log(Pretest ToT)`=pre.avg_time_on_tasks,
    `Math Anxiety`=pre_MA_total_score,
    `Math Neg. React.`=pre_negative_reaction_score,
    `Numerical Conf.`=pre_numerical_confindence_score,
    Posttest=post.total_math_score,Z)%>%
  CreateTableOne(vars=setdiff(names(.),'Z'),data=.,strata="Z")%>%
  print(missing=TRUE,test=FALSE,showAllLevels=TRUE,explain=FALSE)%>%
  cbind(.,
        `Imp. Err.`=round(oob[pmatch(gsub("( \\(%\\))|( \\(mean \\(SD\\)\\))","",rownames(.)),rownames(oob)),],2))



tab1=rbind(
  tab1oth[1:9,],
  tab1bin[-c(1,nrow(tab1bin)),],
  tab1oth[10:nrow(tab1oth),],
  `Bottom Out (%)`=tab1bin[nrow(tab1bin),])

colnames(tab1)[4]='Miss. %'

#print(xtable(tab1),hline.after=17,sanitize.text.function=function(x) {x})

sink('writeUps/table1.tex')
kbl(tab1,format='latex',booktabs=TRUE)%>% #c("p{1.5in}",rep("l",5)))%>%
  group_rows("Baseline",1,17,indent=FALSE)%>%group_rows("Post-Treatment",18,19,indent=FALSE)
sink()


load('results/geeResults.RData')


sink('writeUps/outcomeRegAppendix.tex')
texreg(list(
  psModel=estimates3$BAU$psMod,
  BAU=estimates3$BAU$outMod,
  FH2T=estimates3$FH2T$outMod,
    DragonBox=estimates3$Dragon$outMod),
    file='writeUps/outcomeRegAppendix.tex',
    omit.coef = c('SchIDPre|ClaIDPre'))
