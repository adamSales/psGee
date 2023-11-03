options(knitr.kable.NA = '')

print(load('data/psdat.RData'))
print(load('data/imputations.RData'))

tab1dat <-
  impDat[rownames(psdat),c(intersect(names(impDat),names(psdat)),'nbo')]%>%
  mutate(
    Condition=psdat$trt,
    `Bottom-Outer`=as.factor(ifelse(is.na(psdat$S),0,psdat$S*1)),
    Posttest=psdat$Y,
    Gender=ifelse(Gender=="M",'Male','Female'),
    virtual=ifelse(virtual==1,'Remote','In-Person'),
    across(where(~length(dim(.))>0),~.[,1]))%>%
  select(-contains("ID",ignore.case=FALSE))%>%
  rename(
    `# Bottom-Out`=nbo,
    `Has EIP`=EIP,
    `Has IEP`=IEP,
    Modality=virtual,
    `Gifted`=GIFTED,
    Pretest=pre.total_math_score,
    `Grade 5 Stand. Test`=Scale.Score5,
    `Grade 5 Perf. Lev.`=Performance.Level5,
    `Race/Ethnicity`=raceEth,
    `log(Pretest ToT)`=pre.total_time_on_tasks,
    `Math Anxiety`=pre_MA_total_score,
    `Pretest-Procedural`=pre.sub_P_score,
    `Pretest-Flexibility`=pre.sub_F_score,
    `Pretest-# Completed`=pre.math_completed_num,
    `Math Self-Eff`=pre_MSE_total_score,
    `Perceptual Sens.`=pre_PS_tasks_total_score,
    `Perc. Sens. Pt 2E`=pre_PS_part2E_score,
    `Perc. Sens. Pt 2NE`=pre_PS_part2NE_score,
    `Perc. Sens. #Comp.`=pre_PS_completed_num,
    `log(PS Resp. Time)`=pre_PS_total_RT_sec,
    `log(Days Abs. 5th+1)`=AbsentDays5,
    `log(Days Unexc. 5th+1)`=UnexcusedDays5,
    `log(Days Abs. 6th+1)`=AbsentDays6,
    `log(Days Unexc. 6th+1)`=UnexcusedDays6)%>%
  select( everything(),`# Bottom-Out`,Posttest)







oob=as.data.frame(rbind(imp$OOBerror))%>%
  rename(
    `Has EIP`=EIP,
    `Has IEP`=IEP,
    Modality=virtual,
    `Gifted`=GIFTED,
    Pretest=pre.total_math_score,
    `Grade 5 Stand. Test`=Scale.Score5,
    `Grade 5 Perf. Lev.`=Performance.Level5,
    `Race/Ethnicity`=raceEth,
    `log(Pretest ToT)`=pre.total_time_on_tasks,
    `Math Anxiety`=pre_MA_total_score,
    `Pretest-Procedural`=pre.sub_P_score,
    `Pretest-Flexibility`=pre.sub_F_score,
    `Pretest-# Completed`=pre.math_completed_num,
    `Math Self-Eff`=pre_MSE_total_score,
    `Perceptual Sens.`=pre_PS_tasks_total_score,
    `Perc. Sens. Pt 2E`=pre_PS_part2E_score,
    `Perc. Sens. Pt 2NE`=pre_PS_part2NE_score,
    `Perc. Sens. #Comp.`=pre_PS_completed_num,
    `log(PS Resp. Time)`=pre_PS_total_RT_sec,
    `log(Days Abs. 5th+1)`=AbsentDays5,
    `log(Days Unexc. 5th+1)`=UnexcusedDays5,
    `log(Days Abs. 6th+1)`=AbsentDays6,
    `log(Days Unexc. 6th+1)`=UnexcusedDays6)%>%
  t()




tab1bin <- tab1dat%>%
  select(where(is.factor),Condition)%>%
  select(
    where(~all(levels(.)%in%c(0,1,'TRUE','FALSE'))),#,where(~all(.%in%c('TRUE','FALSE'),na.rm=TRUE)),
                                        #where(~n_distinct(.)<5),
    -starts_with("pre."),Condition)%>%
  select(everything(),`Bottom-Outer`)%>%
  CreateTableOne(vars=setdiff(names(.),'Condition'),data=.,strata="Condition")%>%
  print(missing=TRUE,test=FALSE,dropEqual=TRUE,explain=FALSE)%>%
  cbind(level='',.,`Imputation Error`=round(oob[match(rownames(.),rownames(oob)),1],2))




tab1fac <-  tab1dat[,setdiff(names(tab1dat),rownames(tab1bin))]%>%
  select(
    where(~n_distinct(.)<=5),
    -starts_with("pre."),Condition)%>%
  CreateTableOne(vars=setdiff(names(.),'Condition'),data=.,strata="Condition")%>%
  print(missing=TRUE,test=FALSE,dropEqual=TRUE,explain=FALSE,showAllLevels=TRUE)%>%
  cbind(`Imputation Error`=round(oob[match(rownames(.),rownames(oob)),1],2))



tab1oth<- tab1dat%>%
  select(
    where(~n_distinct(.)>5),starts_with("pre."),Condition)%>%
  CreateTableOne(vars=setdiff(names(.),'Condition'),data=.,strata="Condition")%>%
  print(missing=FALSE,explain=FALSE,test=FALSE)%>%
  cbind(`Imputation Error`=round(oob[match(rownames(.),rownames(oob)),1],2))




tab1a=rbind(
  tab1fac,
  tab1bin[-1,]
)%>%
  cbind(Variable=rownames(.),.)

colnames(tab1a)[c(ncol(tab1a)-1,ncol(tab1a))]=c('Miss. %','Imp. Err. (PFC)')
colnames(tab1a)[1:2]=''
xtable(tab1a,
       caption="Counts and percentages for categorical study variables, by randomized condition. Imputation error is the proportion falsly classified, as estimated using by missForest using out-of-bag observations. All variables were measured at baseline, with the exception of \"Bottom-Outer,\" the principal stratification variable.",
       label="table:tab1fac")%>%
  print(include.rownames=FALSE,file="results/tab1fac.tex",
        hline.after=c(-1,0,which(rownames(tab1a)=='Bottom-Outer')-1,nrow(tab1a)))

colnames(tab1oth)[c(ncol(tab1oth)-1,ncol(tab1oth))]=c('Miss. %','Imp. Err. (NRMSE)')

xtable(tab1oth,
       caption="Means and standard deviations for numeric study variables, by randomized condition. Imputation error is the normalized root mean squared error (NRMSE), as estimated using by missForest using out-of-bag observations. All variables were measured at baseline, with the exception of \"# Bottom-Out\" (the number of bottom-out hints requested) and Posttest.",
       label="table:tab1fac")%>%
  print(include.rownames=TRUE,file="results/tab1num.tex",
        hline.after=c(-1,0,which(rownames(tab1oth)=='# Bottom-Out')-1,nrow(tab1oth)))



load('results/geeResults.RData')


#sink('writeUps/outcomeRegAppendix.tex')
texreg(
  list(
    psModel=estimates1all$BAU$psMod,
    BAU=estimates1all$BAU$outMod,
    FH2T=estimates1all$FH2T$outMod,
    DragonBox=estimates1all$Dragon$outMod),
  file='writeUps/outcomeRegAppendix.tex',
  longtable=TRUE,
  omit.coef = c('SchIDPre|ClaIDPre'))
