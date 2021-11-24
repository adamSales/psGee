### start with deleting missing S

covForm <- ~ splines::ns(log(student_prior_median_time_on_task),
    knots = c(quantile(log(student_prior_median_time_on_task),
         probs = c((1:20)/21)), 6,8,10,12))

## check usage model
modU1 <- glm(update(covForm,S~.),data=dat,family=binomial,
             subset=contains_video&!is.na(S))

bpMod <- function(mod) 
  arm::binnedplot(predict(mod,type='response'),resid(mod,type='response'))

bpMod(modU1)

with(subset(dat,!is.na(S)),
     arm::binnedplot(log(student_prior_median_time_on_task),S))

lines(seq(2,14,length=100),predict(modU1,data.frame(student_prior_median_time_on_task=exp(seq(2,14,length=100))),type='response'),col='red',lwd=2)

with(subset(dat,!is.na(S)),
     arm::binnedplot(log(student_prior_median_time_on_task),
                     resid(modU1,type='response')))

with(subset(dat,!is.na(S)),
     arm::binnedplot(log(student_prior_median_time_on_task),
                     resid(modU1,type='response'),xlim=c(2,6)))

### test out outcome models
dat1 <- datNames(dat,trt = 'contains_video',out='npc',use='S',block='alts')

AUCmod(modU1)

ps <- predict(modU1,dat1,type='response')
dat1 <- within(dat1,Sp <- ifelse(is.na(S),ps,S))
modY <- lm(Yadj~Zadj*Sp+splines::ns(log(student_prior_median_time_on_task),
                                    knots = c(quantile(log(student_prior_median_time_on_task),
                                                       probs = c((1:20)/21)), 6,8,10,12)),data=dat1)

dat2 <- dat1%>%
  filter(!is.na(student_prior_average_correctness))%>%
  ungroup()%>%
  mutate(corFac=cut(student_prior_average_correctness,2))

modY2 <- feols(Y~1#Z*Sp+#ns(log(student_prior_median_time_on_task),
                          #           knots = c(seq(2,14,2)))+
#corFac                 # I(student_prior_average_correctness==0)#+
                 #bs(log(student_prior_started_problem_count),df=5)
                 |block
,data=dat2,vcov='hetero')

bp(fitted(modY2),resid(modY2),na.rm=TRUE)


bpMod(modY)

bp <- arm::binnedplot
with(dat,bp(log(student_prior_median_time_on_task),Yadj))


lines(seq(2,14,length=100),predict(modY,data.frame(student_prior_median_time_on_task=exp(seq(2,14,length=100)),Zadj=0.5,Sp=1),type='response'),col='red',lwd=2)

lines(seq(2,14,length=100),predict(modY,data.frame(student_prior_median_time_on_task=exp(seq(2,14,length=100)),Zadj=-0.5,Sp=1),type='response'),col='blue',lwd=2)

lines(seq(2,14,length=100),predict(modY,data.frame(student_prior_median_time_on_task=exp(seq(2,14,length=100)),Zadj=-0.5,Sp=0),type='response'),col='green',lwd=2)

lines(seq(2,14,length=100),predict(modY,data.frame(student_prior_median_time_on_task=exp(seq(2,14,length=100)),Zadj=0.5,Sp=0),type='response'),col='purple',lwd=2)


### test indepdendece in Z==1
modY11 <- lm(Yadj~S+splines::ns(log(student_prior_median_time_on_task),
                                    knots = c(quantile(log(student_prior_median_time_on_task),
                                                       probs = c((1:20)/21)), 6,8,10,12)),data=subset(dat1,!is.na(S)))

modY12 <- lm(Yadj~S*splines::ns(log(student_prior_median_time_on_task),
                                knots = c(quantile(log(student_prior_median_time_on_task),
                                                   probs = c((1:20)/21)), 6,8,10,12)),data=subset(dat1,!is.na(S)))

anova(modY11,modY12)



pe <- pointEst(dat1,covForm)

bpMod(pe$outMod)
library(splines)
pe2 <- pointEst(dat1,covForm,covFormY=~ns(log(student_prior_median_time_on_task),knots=seq(2,14,2))+ns(student_prior_average_correctness,knots=c(seq(.2,.8,.2),.9))+I(student_prior_average_correctness==0))

bpMod(pe2$outMod)


with(subset(dat,!is.na(student_prior_average_correctness)),
     arm::binnedplot(student_prior_average_correctness,npc))

pe3 <- pointEst(dat,covForm,covFormY=update(covForm,.~.+ns(student_prior_average_correctness,10)+I(student_prior_average_correctness==0)),trt = 'Zadj',out='Yadj',strat='S')

bpMod(pe3$outMod)

estimate <- est(dat,covFormU=covForm,
                covFormY=update(covForm,.~.+ns(student_prior_average_correctness,10)+I(student_prior_average_correctness==0)),trt = 'Zadj',out='Yadj',strat='S')

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


