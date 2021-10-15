library(tidyverse)
library(estimatr)
library(geex)
library(rstan)
library(rootSolve)
library(numDeriv)


dat <- read_csv('assistmentsData/principal_stratification_dataset.csv')

### I THINK: students randomized within problem to different "assigned_tsid"s, some of which contain videos and others don't
## we want problems in which half of the available tsids contain videos

psThirdVid <- dat%>%
    group_by(problem_id,assigned_tsid)%>%
    summarize(vid=mean(contains_video))%>%
    group_by(problem_id)%>%
    summarize(ntsid=n_distinct(assigned_tsid),pVid=mean(vid))%>%
    filter(pVid==1/3)%>%
    pull(problem_id)

sum(dat$problem_id%in%psThirdVid)


### go with Pr(vid)=1/3

dat <- dat%>%
    filter(problem_id%in%psThirdVid)




## differential attrition

### drop problem ids with >50% attrition
dat <- dat%>%group_by(problem_id)%>%mutate(pAtt=mean(is.na(next_problem_correct)))%>%filter(pAtt<0.5)%>%ungroup()

## whatever
dat$S <- dat$time_on_task<dat$video_length

dat <- subset(dat,!is.na(time_on_task))

### try out ols modeling
form2 <- next_problem_correct~log(student_prior_started_problem_count)+log(student_prior_average_attempt_count)+
    student_prior_completed_problem_fraction+student_prior_answer_first_fraction+log(student_prior_median_first_response_time)+
    log(student_prior_median_time_on_task)+splines::ns(student_prior_average_correctness,10)+log(problem_prior_started_problem_count)+
    log(problem_prior_average_attempt_count)+
    problem_prior_completed_problem_fraction+problem_prior_answer_first_fraction+log(problem_prior_median_first_response_time)+
    splines::ns(problem_prior_average_correctness,10)+S

### S model

form3 <- update(form2,.~.-log(problem_prior_median_time_on_task)-log(student_prior_median_time_on_task)+splines::ns(log(problem_prior_median_time_on_task),10)+splines::ns(log(student_prior_median_time_on_task),10))

dat <- droplevels(dat)%>%
    filter((student_prior_average_attempt_count>0),(problem_prior_average_attempt_count>0))

X <- model.matrix(update(form3,.~.-S),data=dat)[,-1]
X <- scale(X)

dat1 <- dat[as.numeric(rownames(X)),]

newDat <- data.frame(cbind(Y=ifelse(dat1$next_problem_correct,1,0),
                           S=ifelse(dat1$S,1,0),
                           Z=ifelse(dat1$contains_video,1,0),
                           X))

newDat$S[newDat$Z==0] <- 0

### overall ATE
ate1 <- with(newDat,
             prop.test(c(sum(Y[Z==1]),sum(Y[Z==0])),c(sum(Z),sum(Z==0))))

### with regression
ate2 <- tidy(lm_robust(Y~.-S,data=newDat))%>%filter(term=='Z')

ate3 <- tidy(lm_lin(Y~Z,
                    covariates=as.formula(paste('~',paste(names(newDat)[-c(1:3)],collapse='+'))),
                    data=newDat)
             )%>%filter(term=='Z')


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
            sum(ZZ*(YY-xb1)),
            sum(ZZ*SS*(YY-xb1)),
            (t(XX)%*%(YY-(ZZ*xb1+(1-ZZ)*(xb2+a01*ps+a00*(1-ps)))))[,1],
### logistic regression
            sum(ZZ*(SS-ps)),
            (t(XX)%*%(ZZ*(SS-ps)))[,1],
### mixture model in control group
            sum((1-ZZ)*(a01*ps+a00*(1-ps)-YY+xb2)),
            sum((1-ZZ)*(a01*ps^2+a00*(ps-ps^2)-(YY-xb2)*ps)),
            length(ZZ)*(eff0-(a10-a00)),
            length(ZZ)*(eff1-(a11+a10-a01)),
            length(ZZ)*(effDiff-(eff1-eff0))
        )
    }
}

psi <- estFun(newDat)

theta <- multiroot(psi,rep(.1,108))

A <- -jacobian(psi,theta$root)
A <- A/nrow(newDat)

ee <- psi(theta$root)
B <- ee%*%t(ee)/nrow(newDat)

Ainv <- solve(A)

vcv <- Ainv%*%B%*%t(Ainv)

