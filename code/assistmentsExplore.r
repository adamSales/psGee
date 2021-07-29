library(tidyverse)
library(mosaic)
library(estimatr)
library(geex)
library(rstan)

dat <- read_csv('assistmentsData/principal_stratification_dataset.csv')

### I THINK: students randomized within problem to different "assigned_tsid"s, some of which contain videos and others don't
## we want problems in which half of the available tsids contain videos

psHalfVid <- dat%>%
    group_by(problem_id,assigned_tsid)%>%
    summarize(vid=mean(contains_video))%>%
    group_by(problem_id)%>%
    summarize(ntsid=n_distinct(assigned_tsid),pVid=mean(vid))%>%
    filter(pVid==0.5)%>%
    pull(problem_id)

psThirdVid <- dat%>%
    group_by(problem_id,assigned_tsid)%>%
    summarize(vid=mean(contains_video))%>%
    group_by(problem_id)%>%
    summarize(ntsid=n_distinct(assigned_tsid),pVid=mean(vid))%>%
    filter(pVid==1/3)%>%
    pull(problem_id)

sum(dat$problem_id%in%psHalfVid)
sum(dat$problem_id%in%psThirdVid)


### go with Pr(vid)=1/3

dat <- dat%>%
    filter(problem_id%in%psThirdVid)


### attrition
prop.test(is.na(next_problem_correct)~contains_video,data=dat)
tally(is.na(next_problem_correct)~contains_video,format='proportion',data=dat)

## differential attrition
dat%>%group_by(problem_id)%>%summarize(pAtt=mean(is.na(next_problem_correct)))%>%pull(pAtt)%>%quantile(c(0,0.001,seq(0.01,0.04,0.01),seq(0.05,0.95,by=0.05),seq(0.96,0.99,0.01),0.999,1),na.rm=TRUE)

### drop problem ids with >50% attrition
dat <- dat%>%group_by(problem_id)%>%mutate(pAtt=mean(is.na(next_problem_correct)))%>%filter(pAtt<0.5)%>%ungroup()
prop.test(is.na(next_problem_correct)~contains_video,data=dat)
tally(is.na(next_problem_correct)~contains_video,format='proportion',data=dat)
## still differential attrition... maybe set NAs to 0? (if they didn't do the next problem, they may as well have gotten it wrong)
## figure it out later


### look at video length (only defined if contains_video?)

xtabs(~contains_video+is.na(video_length),dat)
## yup

### define S: compare time-on-task to video length

## time on task has outliers
hist(dat$time_on_task)
quantile(dat$time_on_task[dat$contains_video],c(0,0.001,seq(0.01,0.04,0.01),seq(0.05,0.95,by=0.05),seq(0.96,0.99,0.01),0.999,1),na.rm=TRUE)
mean(dat$time_on_task[dat$contains_video]<1000,na.rm=T)


xtabs(~contains_video+is.na(time_on_task),dat)
tally(is.na(time_on_task)~contains_video,format='proportion',data=dat)
prop.test(is.na(time_on_task)~contains_video,data=dat)

dat%>%filter(time_on_task<1000)%>%group_by(problem_id)%>%summarize(time=mean(time_on_task[contains_video],na.rm=TRUE),video_length=mean(video_length,na.rm=TRUE))%>%gf_point(time~video_length,data=.)+geom_smooth()

dat%>%group_by(problem_id)%>%summarize(time=mean(log(time_on_task[contains_video]),na.rm=TRUE),video_length=mean(log(video_length),na.rm=TRUE))%>%gf_point(time~video_length,data=.)+geom_smooth()

dat%>%filter(time_on_task<1000,contains_video)%>%ggplot(aes(video_length,time_on_task))+stat_summary_bin(geom='point',fun="mean",bins=500)+geom_smooth(method="lm")

dat%>%filter(contains_video)%>%ggplot(aes(log(video_length),log(time_on_task)))+stat_summary_bin(geom='point',fun="mean",bins=500)+geom_smooth(method="lm")

dat%>%filter(time_on_task<1000)%>%group_by(problem_id)%>%summarize(time=mean(time_on_task[!contains_video],na.rm=TRUE),video_length=mean(video_length,na.rm=TRUE))%>%    gf_point(time~video_length,data=.)+geom_smooth()

boxplot(time_on_task~contains_video,data=filter(dat,time_on_task<1000))
wilcox.test(time_on_task~contains_video,data=dat)
median(time_on_task~contains_video,data=dat,na.rm=TRUE)

dat%>%group_by(problem_id)%>%summarize(timeDiff=mean(log(time_on_task[contains_video]),na.rm=TRUE)-mean(log(time_on_task[!contains_video]),na.rm=TRUE),video_length=mean(log(video_length),na.rm=TRUE))%>%    gf_point(timeDiff~video_length,data=.)+geom_smooth()

### doesn't seem to be much of a relationship between how much MORE time vid students are spending vs text students, and length of vid

gf_histogram(~log(time_on_task/video_length),bins=100,data=dat)

mean(dat$time_on_task<dat$video_length,na.rm=TRUE)

## whatever
dat$S <- dat$time_on_task<dat$video_length

dat <- subset(dat,!is.na(time_on_task))

### try out ols modeling
form <- as.formula(paste('next_problem_correct~',paste(grep('prior',names(dat),value=TRUE),collapse='+')))
mod1 <- lm(update(form,.~.+S),data=dat,subset=contains_video)
mean(fitted(mod1)<0)
mean(fitted(mod1)>1)
arm::binnedplot(fitted(mod1),resid(mod1))
##oops
par(mfrow=c(3,5))
mm <- model.matrix(mod1)
walk(colnames(mm)[-c(1,ncol(mm))],
     function(x) arm::binnedplot(mm[,x],resid(mod1),main=x))

walk(colnames(mm)[-c(1,ncol(mm))],
     function(x) hist(mm[,x],main=x))
par(mfrow=c(1,1))


form2 <- next_problem_correct~log(student_prior_started_problem_count)+log(student_prior_average_attempt_count)+
    student_prior_completed_problem_fraction+student_prior_answer_first_fraction+log(student_prior_median_first_response_time)+
    log(student_prior_median_time_on_task)+splines::ns(student_prior_average_correctness,10)+log(problem_prior_started_problem_count)+
    log(problem_prior_average_attempt_count)+
    problem_prior_completed_problem_fraction+problem_prior_answer_first_fraction+log(problem_prior_median_first_response_time)+
    splines::ns(problem_prior_average_correctness,10)+S

mod2 <- lm(form2,data=dat,subset=contains_video&(student_prior_average_attempt_count>0)&(problem_prior_average_attempt_count>0))
mean(fitted(mod2)<0)
mean(fitted(mod2)>1)
arm::binnedplot(fitted(mod2),resid(mod2))

par(mfrow=c(3,5))
mm <- model.matrix(mod2)
walk(grep('spline',colnames(mm)[-c(1,ncol(mm))],value=TRUE,invert=TRUE),
     function(x) arm::binnedplot(mm[,x],resid(mod2),main=x))
arm::binnedplot(dat$student_prior_average_correctness[as.numeric(rownames(mm))],resid(mod2))
arm::binnedplot(dat$problem_prior_average_correctness[as.numeric(rownames(mm))],resid(mod2))
par(mfrow=c(1,1))

### overall ATE
prop.test(next_problem_correct~contains_video,data=subset(dat,!is.na(next_problem_correct)&!is.na(contains_video)))
tally(next_problem_correct~contains_video,data=subset(dat,!is.na(next_problem_correct)&!is.na(contains_video)),format="percent")

### with regression
lm_robust(update(form2,.~.-S+contains_video),data=dat,subset=(student_prior_average_attempt_count>0)&(problem_prior_average_attempt_count>0))

### S model
Smod1 <- glm(update(form2,S~.-S),subset=contains_video&(student_prior_average_attempt_count>0)&(problem_prior_average_attempt_count>0),family=binomial,data=dat)

pred <- function(mod) predict(mod,type='response')
res <- function(mod) mod$y-pred(mod)

arm::binnedplot(pred(Smod1),res(Smod1))

par(mfrow=c(3,5))
mm <- model.matrix(Smod1)
walk(grep('spline',colnames(mm)[-c(1,ncol(mm))],value=TRUE,invert=TRUE),
     function(x) arm::binnedplot(mm[,x],res(Smod1),main=x))
arm::binnedplot(dat$student_prior_average_correctness[as.numeric(rownames(mm))],resid(mod2))
arm::binnedplot(dat$problem_prior_average_correctness[as.numeric(rownames(mm))],resid(mod2))
par(mfrow=c(1,1))

form3 <- update(form2,.~.-log(problem_prior_median_time_on_task)-log(student_prior_median_time_on_task)+splines::ns(log(problem_prior_median_time_on_task),10)+splines::ns(log(student_prior_median_time_on_task),10))

Smod2 <- glm(update(form3,S~.-S),subset=contains_video&(student_prior_average_attempt_count>0)&(problem_prior_average_attempt_count>0),family=binomial,data=dat)

arm::binnedplot(pred(Smod2),res(Smod2))

## let's just use that MM

dat <- droplevels(dat)%>%
    filter((student_prior_average_attempt_count>0),(problem_prior_average_attempt_count>0))

X <- model.matrix(update(form3,.~.-S),data=dat)[,-1]
X <- scale(X)

dat1 <- dat[as.numeric(rownames(X)),]

newDat <- data.frame(cbind(Y=ifelse(dat1$next_problem_correct,1,0),
                           S=ifelse(dat1$S,1,0),
                           Z=ifelse(dat1$contains_video,1,0),
                           X))

estFun <- function(data){
    XX <- as.matrix(select(data,-Y,-S,-Z))

    YY <- data$Y

    ZZ <- data$Z
    SS <- data$S

    p <- ncol(XX)
    function(theta){
        ## YY model (ZZ=1): YY=a10+a11SS+a2XX
        ##         (ZZ=0): YY={a00 or a01}+a2XX
        ## SS model (ZZ=1): logit(SS)=b0+b1XX
        ## effects: a10-a00 and a10+a11-a01

        a10 <- theta[1]
        a11 <- theta[2]
        a00 <- theta[3]
        a01 <- theta[4]
        eff0 <- theta[5] ## a10-a00
        eff1 <- theta[6] ## a10+a11-a01
        effDiff <- theta[7] ## a11+a00-a01
        a2 <- theta[8:(p+7)]
        b0 <- theta[p+8]
        b1 <- theta[(p+9):(2*p+8)]

        ## ols z=1
        ## print(length(SS))
        ## print(dim(XX))
        xb1 <- (cbind(1,SS,XX)%*%c(a10,a11,a2))[,1]

        ## xb for z=0
        xb2 <- (XX%*%a2)[,1] ## (a2 is same in trt groups)

        ## logit z=1
        xb3 <- (cbind(1,XX)%*%c(b0,b1))[,1]
        ps <- plogis(xb3)

        c(
### regression for treatment group
            sum(ifelse(ZZ==1,YY-xb1,0)),
            sum(ifelse(ZZ==1,SS*(YY-xb1),0)),
            as.vector(apply(XX,2,function(x)
                sum(ifelse(ZZ==1,x*(YY-xb1),x*(YY-xb2-a01*ps-a00*(1-ps)))))),
### logistic regression
            sum(ifelse(ZZ==1,SS-ps,0)),
            as.vector(apply(XX,2,function(x) sum(ifelse(ZZ==1,x*(SS-ps),0)))),
### mixture model in control group
            sum(ifelse(ZZ==0,a01*ps+a00*(1-ps)-YY+xb2,0)),
            sum(ifelse(ZZ==0,a01*ps^2+a00*(ps-ps^2)-(YY-xb2)*ps,0)),
            eff0-(a10-a00),
            eff1-(a11+a10-a01),
            effDiff-(eff1-eff0)
        )
    }
}

mod1 <- lm(Y~.-Z,data=newDat,subset=Z==1)
mod2 <- glm(S~.-Z-Y,data=newDat,subset=Z==1,family=binomial)

time <- system.time(
    result <- m_estimate(estFun,newDat,root_control = setup_root_control(
                                           start = c(coef(mod1)[1],
                                                     coef(mod1)['S'],
                                                     rep(.1,5), #a0*, effects, effDiff
                                                     coef(mod1)[-c(1,which(names(coef(mod1))=='S'))],
                                                     coef(mod2))
                                       )
                         )
)
save(result,file='mestAssist.RData')

