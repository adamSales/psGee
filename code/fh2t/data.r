library(tidyverse)
library(arm)
library(missForest)
select <- dplyr::select

Scale <- function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)

dat=read_csv('data/DATA20220202_3591.csv',na=c('','NA','#NULL!'))

kirk <- read_csv('data/data_uptated_8_18.csv')

dat$nbo <- kirk$num_problem_parts_used_bottom_out_hint[match(dat$StuID,kirk$StuID)]
dat <- filter(dat,!is.na(post.total_math_score))


cat(sum(is.na(dat$nbo)&dat$rdm_condition=='ASSISTments'),' students with NA for # bottom out hints; setting to 0\n')
dat$nbo[is.na(dat$nbo)&dat$rdm_condition=='ASSISTments'] <- 0

plot(table(dat$nbo[dat$rdm_condition=='ASSISTments']),xlim=c(0,30))
med <- median(dat$nbo[dat$rdm_condition=='ASSISTments'],na.rm=TRUE)

dat$S <- dat$nbo>med

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

impDat <- psdat%>%
  select(virtual, Gender,raceEth,Performance.Level5,EIP,Scale.Score5:PercentInAttendance6,starts_with("pre",ignore.case=FALSE))%>%
  select(where(~mean(is.na(.))<0.2))%>%
  select(where(~n_distinct(.,na.rm=TRUE)>1))%>%
  map_dfc(~if(is.character(.)|n_distinct(.,na.rm=TRUE)<4) as.factor(.) else .)%>%
  as.data.frame()

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

imp <- missForest(impDat,variablewise = TRUE)

save(imp,file='data/imputations.RData')

mse=which(names(imp$OOBerror)=='MSE')
names(imp$OOBerror) <- names(impDat)

r2imp=1-imp$OOBerror[mse]/sapply(impDat[,mse],var,na.rm=TRUE)

#names(imp$ximp) <- paste0(names(imp$ximp),'Imp')

naFlags <- map_dfc(impDat,is.na)%>%select(where(~sum(.)>1))%>%as.data.frame()
V <- matrix(nrow=ncol(naFlags)-1,ncol=ncol(naFlags)-1,dimnames=list(names(naFlags)[-1],names(naFlags)[-ncol(naFlags)]))
for(i in 2:(ncol(naFlags))) for(j in 1:(i-1)) V[i-1,j] <- cramerV(naFlags[,i],naFlags[,j])

naFlags <- naFlags[,!names(naFlags)%in%rownames(which(V>0.8,arr.ind=TRUE))]
names(naFlags) <- paste0(names(naFlags),'NA')

psdat <- cbind(psdat[,c('StuID','SchIDPre','TeaIDPre_within_school','ClaIDPre','trt','S','Y')],
               imp$ximp,
               naFlags)

psdat <- psdat%>%
  rename(TeaIDPre=TeaIDPre_within_school)%>%
  mutate(across(ends_with('IDPre'),as.factor))

psdat$TeaIDPre <- relevel(psdat$TeaIDPre,ref='78') ## biggest in ASSISTments group

save(psdat,file='data/psdat.RData')
