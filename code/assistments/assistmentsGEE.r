if(!all(c('dat1')%in%ls())) load('assistmentsData/assistmentsAnalysisData.RData')

covForm <- ~ splines::ns(log(student_prior_median_time_on_task),
    knots = c(quantile(log(student_prior_median_time_on_task),
         probs = c((1:10)/11)), 6))

estimate <- est(dat1,covForm=covForm,int=FALSE)
effsFromFit(estimate)

binnedplot(estimate$psMod$fitted,residuals(estimate$psMod,type='response'))


est2 <- est(dat1,covForm,int=TRUE)
effsFromFit(est2)

save(estimate,file=paste0('assistmentsResults/Mestimate',Sys.Date(),'.RData')

try(system('drive push --quiet assistmentsResult.RData'))


## try random forest ps model
rf <- randomForest(factor(S)~student_prior_started_problem_count + student_prior_average_attempt_count +
    student_prior_completed_problem_fraction + student_prior_answer_first_fraction +
    student_prior_median_first_response_time + student_prior_average_correctness
     + problem_prior_started_problem_count + problem_prior_average_attempt_count +
    problem_prior_completed_problem_fraction + problem_prior_answer_first_fraction +
    problem_prior_median_first_response_time + problem_prior_average_correctness
     + problem_prior_median_time_on_task
    + student_prior_median_time_on_task,data=subset(dat1,Z==1))

ps1 <- predict(rf,type='prob')[,'1']

arm::binnedplot(ps1,(rf$y=='1')-ps1)

plot(ps1,estimate$psMod$fitted.values)

ps <- predict(rf,dat1,type='prob')
Sp <- ifelse(dat1$Z==1,dat1$S,ps)

outMod2 <- lm(formula(estimate$outMod),data=within(dat1,Sp <- Sp))

#effsFromFit(list(outMod=outMod2,vcv=sandwich(outMod2)))arm::binnedplot(predict

outlier <- function(x,thresh=1.5){
    qq <- quantile(x,c(0.25,0.75),na.rm=TRUE)
    iqr <- qq[2]-qq[1]
    (x>qq[2]+thresh*iqr)|(x<qq[1]-thresh*iqr)
}


psMod3 <- glm(S~
                log(student_prior_started_problem_count)+
                log(student_prior_median_first_response_time)+
                log(problem_prior_average_attempt_count)+
                problem_prior_completed_problem_fraction+
                log(problem_prior_median_first_response_time)+
                outlier(log(student_prior_started_problem_count))+
                outlier(log(student_prior_median_first_response_time))+
                outlier(log(problem_prior_average_attempt_count))+
                outlier(problem_prior_completed_problem_fraction)+
                outlier(log(problem_prior_median_first_response_time)),
            dat=subset(dat1,Z==1),family=binomial(cloglog))
estimate3 <- est(dat1,psMod=psMod3)
effs3 <- effsFromFit(estimate3)

with(estimate3,arm::binnedplot(fitted(outMod),resid(outMod)))
with(estimate3,arm::binnedplot(psMod$fitted,psMod$y-psMod$fitted))


binnedLink <- function(x,y,linkfun,RESID=TRUE,...){
    br <- arm::binned.resids(x,y)[[1]]

    newX <- br[,'xbar']
    newY <- linkfun(br[,'ybar'])
    if(RESID) newY <- newY-newX

    data.frame(x=newX,y=newY)
}


