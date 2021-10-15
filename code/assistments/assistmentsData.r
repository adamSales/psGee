
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



save(X,dat1,newDat,ate1,ate2,ate3,file='assistmentsData/assistmentsAnalysisData.RData')
