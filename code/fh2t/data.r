Scale <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

dat=read_csv('data/DATA20220202_3591.csv',na=c('','NA','#NULL!'))

kirk <- read_csv('data/data_uptated_8_18.csv')

dat$nbo <- kirk$num_problem_parts_used_bottom_out_hint[match(dat$StuID,kirk$StuID)]
dat <- filter(dat,!is.na(post.total_math_score))

### set bottom out hints=NA to 0
cat(sum(is.na(dat$nbo)&dat$rdm_condition=='ASSISTments'),
    ' students with NA for # bottom out hints; setting to 0\n')
dat$nbo[is.na(dat$nbo)&dat$rdm_condition=='ASSISTments'] <- 0


## make binary w/ a median split
med <- median(dat$nbo[dat$rdm_condition=='ASSISTments'],na.rm=TRUE)

dat$S <- #dat$pbo4>0.1#
  dat$nbo>med

## get covariates, treatment, outcome data
psdat <- dat%>%
  select(StuID:Performance.Level5,EIP:PercentInAttendance6,starts_with("pre"),Y=post.total_math_score,S)%>%
  mutate(
    raceEth=raceEthnicityFed%>%
      factor()%>%
      fct_lump_min(200)%>%
      fct_recode(`Hispanic/Latino`="1",Asian="3",White="6")%>%
      fct_relevel('White'),
    Gender=as.factor(Gender))%>%
  rename(trt=rdm_condition)

### impute missing covariates w/ missForest
impDat <- psdat%>%
  select(virtual, Gender,raceEth,Performance.Level5,EIP,Scale.Score5:PercentInAttendance6,starts_with("pre",ignore.case=FALSE))%>%
  select(where(~mean(is.na(.))<0.2))%>%
  select(where(~n_distinct(.,na.rm=TRUE)>1))%>%
  map_dfc(~if(is.character(.)|n_distinct(.,na.rm=TRUE)<4) as.factor(.) else .)%>%
  as.data.frame()

### find highly-correlated covariates
Rho=cor(impDat[,sapply(impDat,class)=='numeric'],method='spearman',use='pairwise')
Rho[upper.tri(Rho,diag=TRUE)] <- NA
bigCor=which(!is.na(Rho)&abs(Rho)>0.8,arr.ind=TRUE)
sapply(1:nrow(bigCor),function(i) c(rownames(Rho)[bigCor[i,1]],colnames(Rho)[bigCor[i,2]],round(Rho[bigCor[i,1],bigCor[i,2]],2)))%>%t()

impDat <- impDat[,!names(impDat)%in%rownames(bigCor)]
Rho=cor(impDat[,sapply(impDat,class)=='numeric'],method='spearman',use='pairwise')
Rho[upper.tri(Rho,diag=TRUE)] <- NA
round(tail(sort(abs(Rho))),3)

### now for factors
cramerV = function(x, y){
  N      = length(x)
  Chi.sq = suppressWarnings(chisq.test(x, y, correct=FALSE)$statistic)
  Phi    = Chi.sq / N
  Row    = length(unique(x))
  C      = length(unique(y))
  sqrt(Phi / min(Row-1, C-1))
}

fac=which(sapply(impDat,class)=='factor')
V <- matrix(nrow=length(fac)-1,ncol=length(fac)-1,dimnames=list(names(fac)[-1],names(fac)[-length(names(fac))]))
for(i in 2:(length(fac))) for(j in 1:(i-1)) V[i-1,j] <- cramerV(impDat[,fac[i]],impDat[,fac[j]])
### keep all factors

### some logical stuff
impDat$pre.math_completed_num[is.na(impDat$pre.total_math_score)] <- 0
impDat$pre.total_time_on_tasks[is.na(impDat$pre.total_math_score)] <- NA
impDat$pre_PS_completed_num[is.na(impDat$pre_PS_tasks_total_score)] <- 0
## idea: they didn't take pretest, but it'd be interesting to know how much time they would have spent

### check on distributions of predictors
table(sapply(impDat,class))
par(mfrow=c(4,5))
impDat%>%select(where(is.numeric))%>%
  iwalk(~hist(.x,main=.y))
par(mfrow=c(1,1))

### make some transformations
impDat <- impDat%>%
  mutate(
    across(contains("Days",ignore=FALSE),~log(.+1)),
    fullYear5=as.factor(MEMBERSHIPDAYS5==max(MEMBERSHIPDAYS5,na.rm=TRUE)),
    fullYear6=as.factor(MEMBERSHIPDAYS6==max(MEMBERSHIPDAYS6,na.rm=TRUE)),
    pre.total_time_on_tasks=log(pre.total_time_on_tasks),
    pre_PS_total_RT_sec=log(pre_PS_total_RT_sec))%>%
  select(-starts_with("MEMBER"))

### does this help?
table(sapply(impDat,class))
par(mfrow=c(4,5))
impDat%>%select(where(is.numeric))%>%
  iwalk(~hist(.x,main=.y))
par(mfrow=c(1,1))

impDat$noUnexcused5 <- as.factor(impDat$UnexcusedDays5==0)

### impute missing values (takes a while)
set.seed(613)
imp <- missForest(impDat,variablewise = TRUE)

mse=which(names(imp$OOBerror)=='MSE')
names(imp$OOBerror) <- names(impDat)

r2imp=1-imp$OOBerror[mse]/sapply(impDat[,mse],var,na.rm=TRUE)

impDat$nbo <- dat$nbo

save(imp,mse,r2imp,impDat,file='data/imputations.RData')

#names(imp$ximp) <- paste0(names(imp$ximp),'Imp')

### create NA flags
naFlags <- map_dfc(impDat,is.na)%>%select(where(~sum(.)>1))%>%as.data.frame()
V <- matrix(nrow=ncol(naFlags)-1,ncol=ncol(naFlags)-1,dimnames=list(names(naFlags)[-1],names(naFlags)[-ncol(naFlags)]))
for(i in 2:(ncol(naFlags))) for(j in 1:(i-1)) V[i-1,j] <- cramerV(naFlags[,i],naFlags[,j])

naFlags <- naFlags[,!names(naFlags)%in%rownames(which(V>0.8,arr.ind=TRUE))]
names(naFlags) <- paste0(names(naFlags),'NA')

### put it all together
psdat <- cbind(psdat[,c('StuID','SchIDPre','TeaIDPre_within_school','ClaIDPre','trt','S','Y')],
               imp$ximp,
               naFlags)

psdat <- psdat%>%
  rename(TeaIDPre=TeaIDPre_within_school)%>%
  mutate(across(ends_with('IDPre'),as.factor))

psdat$TeaIDPre <- relevel(psdat$TeaIDPre,ref='78') ## biggest in ASSISTments group

### get rid of covariates with only one level
oneLev <- psdat%>%filter(trt=='ASSISTments')%>%droplevels()%>%select(-StuID,-TeaIDPre,-ClaIDPre,-trt,-S,-Y)%>%select(where(is.factor))%>%select(where(~nlevels(.)<2))%>%names()

psdat <- psdat[,-which(names(psdat)%in%oneLev)]

### the pre_PS_part variables all add up to the total score
psdat$pre_PS_part1_score <- NULL

### some factor variables have levels with very few
psdat <- psdat%>%select(where(~n_distinct(.)>2|min(table(.[psdat$trt=='ASSISTments']))>9))

psdat <- psdat%>%mutate(across(c(Scale.Score5,pre.total_time_on_tasks,pre_MSE_total_score,pre_PS_tasks_total_score),scale))

### only 1 student in school 7
psdat <- filter(psdat,SchIDPre!=7)%>%droplevels()

save(psdat,file='data/psdat.RData')
