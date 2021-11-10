dat <- read_csv('assistmentsData/principal_stratification_dataset.csv')

### I THINK: students randomized within problem to different "assigned_tsid"s, some of which contain videos and others don't
## we want problems in which half of the available tsids contain videos

dat%>%group_by(problem_id)%>%summarize(nCond=n_distinct(assigned_tsid))%>%pull(nCond)%>%table()

dat%>%group_by(problem_id,assigned_tsid)%>%summarize(vid=n_distinct(contains_video))%>%pull(vid)%>%table()

table(dat$contains_video)

mean(dat$contains_video)


### how many conditions in each problem have videos?

vidConds <- dat%>%
    group_by(problem_id,assigned_tsid)%>%
    summarize(vid=contains_video[1],nCond=n())%>%
    group_by(problem_id)%>%
    summarize(num_cond=n(),nTot=sum(nCond),perCondVid=mean(vid),perStudVid=sum(nCond[vid])/nTot)

intTS <- function(ts1,ts2){
    rng1 <- range(ts1)
    rng2 <- range(ts2)
    int <- as.numeric(min(c(rng1[2], rng2[2]))-max(c(rng1[1],rng2[1])))
    if(int<0) return(0)
    int
}

overlap <- function(ts1,ts2){
    rng1 <- range(ts1)
    rng2 <- range(ts2)
    int <- as.numeric(min(c(rng1[2], rng2[2]))-max(c(rng1[1],rng2[1])))
    if(int<=0) return(0)
    diff <- as.numeric(sum(abs(rng1-rng2)))

    1-diff/int
}

#### next:
## for each problem id, look for the max intTS between vid and non-vid tsid's
## keep the pair with the max
## delete ovservations outside the int period

tsidInt <- function(ts0,ts1){
    mat <- matrix(0,nrow=length(ts0),ncol=length(ts1))
    for(i in 1:length(ts0))
        for(j in 1:length(ts1))
            mat[i,j] <- intTS(ts0[[i]],ts1[[j]])
    if(all(mat==0)) return(NULL)
    mmm <- which(mat==max(mat),arr.ind=TRUE)
    c(names(ts0)[mmm[1,'row']],names(ts1)[mmm[1,'col']])
}

psProc <- function(psDat){

    ts <- split(psDat,psDat$contains_video)%>%
        map(~split(.$timestamp,.$assigned_tsid))

    tInd <- tsidInt(ts[['FALSE']],ts[['TRUE']])
    if(!length(tInd)) return(NULL)
    psDat%>%
        filter(assigned_tsid%in%tInd)%>%
        mutate(
            MIN=max(min(timestamp[contains_video]),min(timestamp[!contains_video])),
            MAX=min(max(timestamp[contains_video]),max(timestamp[!contains_video]))
        )%>%
        filter(between(timestamp,MIN[1],MAX[1]))
}


datOvr <- dat%>%
    group_by(problem_id)%>%
    mutate(vv=var(contains_video))%>%
    filter(vv>0)%>%
    group_split()%>%
    map(psProc)%>%
    reduce(bind_rows)

mean(datOvr$contains_video)

datOvr <- datOvr%>%
    group_by(problem_id)%>%
    mutate(nCond=n_distinct(assigned_tsid),
           nVid=n_distinct(assigned_tsid[contains_video]),
           nStud=n(),
           perVid=mean(contains_video))%>%
    filter(nCond>1,nStud>100)%>%
    ungroup()


byProb <- datOvr%>%
    group_by(problem_id)%>%
    summarize(nCond=n_distinct(assigned_tsid),
              nVid=n_distinct(assigned_tsid[contains_video]),
              nStud=n(),
              perVid=mean(contains_video),
              Sna=mean(is.na(time_on_task)),
              SnaVid=mean(is.na(time_on_task[contains_video])),
              diffNA=SnaVid-mean(is.na(time_on_task[!contains_video])))

xtabs(~nCond+nVid,data=byProb)
plot(byProb$nStud,byProb$perVid)


### missing time_on_task -> didn't finish problem -> didn't do next problem
xtabs(~is.na(time_on_task)+is.na(next_problem_correct),data=datOvr)


## 8-schools style analysis
library(rstan)
options(mc.cores = parallel::detectCores())

effsByProb <- datOvr%>%
    group_by(problem_id)%>%
    summarize(
        pVid1=mean(next_problem_correct[contains_video],na.rm=TRUE),
        pTxt1=mean(next_problem_correct[!contains_video],na.rm=TRUE),
        eff1=pVid1-pTxt1,
        nVid1=sum(!is.na(next_problem_correct[contains_video])),
        nTxt1=sum(!is.na(next_problem_correct[!contains_video])),
        n1=nVid1+nTxt1,
        n2=n(),
        nVid2=sum(contains_video),
        nTxt2=n2-nVid2,
        pVid2=sum(next_problem_correct[contains_video],na.rm=TRUE)/nVid2,
        pTxt2=sum(next_problem_correct[!contains_video],na.rm=TRUE)/nTxt2,
        eff2=pVid2-pTxt2,
        vVid1=pVid1*(1-pVid1)/nVid1,
        vTxt1=pTxt1*(1-pTxt1)/nTxt1,
        se1=sqrt(vVid1+vTxt1),
        vVid2=pVid2*(1-pVid2)/nVid2,
        vTxt2=pTxt2*(1-pTxt2)/nTxt2,
        se2=sqrt(vVid2+vTxt2)
    )

with(effsByProb,plot(n,eff2))

sdat <- with(filter(effsByProb,vVid2*vTxt2>0),
             list(
                 J=length(eff2),
                 y=eff2,
                 sigma=se2
             )
             )

mod <- stan('code/assistments/schools.stan',data=sdat)

sss <- summary(mod)$summ
theta <- sss[paste0('theta[',1:891,']'),]

hist(sss[,'Rhat'])
hist(sss[,'n_eff'])

theta <- theta[order(theta[,'mean']),]
Hmisc::errbar(1:891,theta[,'mean'],theta[,'97.5%'],theta[,'2.5%'])

max(theta[,'2.5%'])
min(theta[,'97.5%'])

round(sss['mu',],4)
round(sss['tau',],4)

save(mod,sss,sdat,file='assistmentsResults/rubin8schoolsTypeModel.RData')

