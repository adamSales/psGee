library(splines)

datBig <- filter(dat,alts==byExp$alts[which.max(byExp$n)])

### randomization prob
mean(datBig$contains_video)

### does this make sense?
pbinom(sum(datBig$contains_video),nrow(datBig),1/3)


### missing S
mean(is.na(datBig$time_on_task))
datBig%>%group_by(contains_video)%>%summarize(mean(is.na(time_on_task)))
feols(is.na(time_on_task)~contains_video,data=datBig,vcov='hetero')%>%summary()

### is there an overall effect?
feols(npc~contains_video,data=datBig,vcov='hetero')%>%summary()

### what if we delete with missing time_on_task?
feols(npc~contains_video,data=subset(datBig,!is.na(time_on_task)),vcov='hetero')%>%summary()

datBig <- filter(datBig,!is.na(time_on_task))

datBig$ltot <- log(datBig$student_prior_median_time_on_task)

bp <- arm::binnedplot

with(filter(datBig,!is.na(S)),bp(ltot,S))

modU1 <- glm(S~bs(ltot,knots=c(2.5,3.3,4,4.5)),data=datBig,family=binomial)

bp(predict(modU1,type='response'),resid(modU1,type='response'))

modU2 <- glm(S~bs(ltot,df=10),data=datBig,family=binomial)
bp(predict(modU2,type='response'),resid(modU2,type='response'))
AUCmod(modU2)


modU3 <- update(modU2,.~bs(ltot,df=15))
bp(predict(modU3,type='response'),resid(modU3,type='response'))
AUCmod(modU3)
anova(modU2,modU3,test = 'Chisq')


### test out outcome models
ps <- predict(modU2,datBig,type='response')
datBig <- within(datBig,Sp <- ifelse(is.na(S),ps,S))

modY <- lm(npc~contains_video*Sp+bs(ltot,df=10),data=datBig)
bp(fitted(modY),resid(modY))

with(subset(datBig,!is.na(student_prior_average_correctness)),bp(ltot,npc))

### what about adding in prior correctness?

with(subset(datBig,!is.na(student_prior_average_correctness)),bp(student_prior_average_correctness,npc))

modY2 <- update(modY,.~.#-bs(ltot,df=10)+ns(ltot,df=5)+
                  +ns(student_prior_average_correctness,df=10)+
                  I(student_prior_average_correctness==0)+
                  I(student_prior_average_correctness==1))
bp(fitted(modY2),resid(modY2))

summary(modY)
summary(modY2)

with(subset(datBig,!is.na(student_prior_average_correctness)),bp(student_prior_average_correctness,npc))
lines(seq(0,1,length=100),
      predict(modY2,
              data.frame(
                Sp=0,
                contains_video=FALSE,
                ltot=mean(datBig$ltot),
                student_prior_average_correctness=seq(0,1,length=100)
              )))

with(subset(datBig,!is.na(student_prior_average_correctness)),bp(log(student_prior_started_problem_count),npc))


## OK time to estimate stuff!
datBig <- subset(datBig,!is.na(student_prior_average_correctness))
estimate1 <- 
  est(data = datBig,
     covFormU=~bs(ltot,df=10),
     covFormY=~bs(ltot,df=10)+
       ns(student_prior_average_correctness,df=10)+
       I(student_prior_average_correctness==0)+
       I(student_prior_average_correctness==1),
     trt="contains_video",
     out="npc",
     use="S")

effsFromFit(estimate1)    

effs(data = datBig,
covFormU=~bs(ltot,df=10),
covFormY=~bs(ltot,df=10)+
ns(student_prior_average_correctness,df=10)+
I(student_prior_average_correctness==0)+
I(student_prior_average_correctness==1),
trt="contains_video",
out="npc",
use="S")

### bootstrap SEs
BS <- replicate(1000,
                effs(data = datBig[sample(1:nrow(datBig),nrow(datBig),
                                          replace=TRUE),],
                     covFormU=~bs(ltot,df=10),
                     covFormY=~bs(ltot,df=10)+
                       ns(student_prior_average_correctness,df=10)+
                       I(student_prior_average_correctness==0)+
                       I(student_prior_average_correctness==1),
                     trt="contains_video",
                     out="npc",
                     use="S"))

