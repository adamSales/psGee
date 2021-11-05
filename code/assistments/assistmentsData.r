
dat <- read_csv('assistmentsData/principal_stratification_dataset.csv')

### I THINK: students randomized within problem to different "assigned_tsid"s, some of which contain videos and others don't
## we want problems in which half of the available tsids contain videos

dat%>%group_by(problem_id)%>%summarize(nCond=n_distinct(assigned_tsid))%>%pull(nCond)%>%table()

dat%>%group_by(problem_id,assigned_tsid)%>%summarize(vid=n_distinct(contains_video))%>%pull(vid)%>%table()

table(dat$contains_video)

mean(dat$contains_video)

dat <-
    dat%>%
    group_by(problem_id)%>%
    mutate(nExp=n())%>%
    group_by(assigned_tsid)%>%
    mutate(nCond=n(),perCond=nCond/nExp)%>%
    ungroup()%>%
    filter(nCond>19,perCond>0.2,nExp>100)%>%
    group_by(problem_id)%>%
    mutate(
        nstud=n(),
        nstudVid=sum(contains_video),
        vidCond=n_distinct(assigned_tsid[contains_video]),
        totCond=n_distinct(assigned_tsid),
        perVid=vidCond/totCond,
        vidAss=nstudVid/nstud)%>%
## group_by(perVid)%>%
##      summarize(
##          nStud=sum(nstud),
##          nProb=n(),
##          probAss=weighted.mean(vidAss,nstud)
##          )
    filter(perVid==0.5)%>%
filter(btw(nstudVid,qbinom(c(0.025,0.975),nstud,0.5)))%>%
ungroup()


## differential attrition

dat%>%group_by(contains_video)%>%summarize(mean(is.na(next_problem_correct)))

### drop problem ids with >50% attrition
dat <- dat%>%group_by(problem_id)%>%mutate(pAtt=mean(is.na(next_problem_correct)))%>%filter(pAtt<0.5)%>%ungroup()

dat%>%group_by(contains_video)%>%summarize(mean(is.na(next_problem_correct)))

## whatever

dat$S <- dat$time_on_task>=dat$video_length

dat <- subset(dat,!is.na(time_on_task))

dat%>%group_by(contains_video)%>%summarize(mean(is.na(next_problem_correct)))
with(dat,prop.test(table(contains_video,is.na(next_problem_correct))))

prop.test(table(dat$contains_video))

dat1 <- dat%>%
    filter(!is.na(student_prior_median_time_on_task),student_prior_median_time_on_task>0)%>%
    mutate(Y=ifelse(next_problem_correct,1,0), ### name Y, Z, S, make numeric
           Z=ifelse(contains_video,1,0),
           S=ifelse(S&(Z==1),1,0))%>%
    filter(!is.na(Y),!is.na(Z))


RItools::xBalance(contains_video~student_prior_started_problem_count+student_prior_average_attempt_count+student_prior_completed_problem_fraction+student_prior_answer_first_fraction+student_prior_average_correctness+student_prior_median_first_response_time+student_prior_median_time_on_task,strata=list(prob=~problem_id),report='all',data=dat1)

save(dat1,file='assistmentsData/assistmentsAnalysisData.RData')

