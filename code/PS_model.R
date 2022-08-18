#### Principal Startification Analysis ####



####Installing & Loading Packages###
#create list of packages
packages = c(
  "tidyverse",
  "plyr",
  "ggExtra",
  "xts",
  "lubridate",
  "readxl",
  "data.table",
  "RSQLite",
  "DBI",
  "mice"
  ,"psych",
  "stringr",
  "sjmisc",
  "lme4",
  "campfin",
  "pROC",
  "cutpointr",
  "missForest",
  "corrplot",
  "tidyverse",
  "arm"
)
#load install
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  } 
) 
rm(package.check, packages)

ifnull <- function(x,y) ifelse(is.na(x), y, x)
ifnull0 <- function(x) ifelse(is.na(x), 0, x)

options(scipen = 100)
check <- function(X) { list(
  class = class(X),
  distint_count = length(unique(X)),
  (if  (class(X) == "numeric")
  {
    descibe = psych::describe(X)
  } else {
    table = head(table(X), 20)
  }),
  
  is_na=table(is.na(X)))
}
options(digits.secs = 3) 

### Connect to SQLite DB ####
# set db path
ies_research_con<- dbConnect(RSQLite::SQLite(), "/Users/kirkvanacore/Documents/WPI Analyses/MAPLE_IES_DB_Creation/ies_research schema/maple_ies_research.db")

#### Load Data ####

dat <-dbGetQuery(ies_research_con, "
   select cws.StuID,
        cws.assistments_user_id,
        cws.condition_assignment,
        sr.SchIDEnd,
        ifnull(sr.TeaIDEnd, sr.TeaIDPre) as TeaIDEnd,
        sr.virtual,
        sr.inperson,
        sd.raceEthnicityFed,
        sd.hispanicEthnicity,
        sd.IEP,
        sd.EIP,
        sd.ESOL,
     /* sd.ESOL_FORMER, to much missing data to be useful */
        sd.GIFTED,
     /* att.PresentDays6,
        att.AbsentDays6,
        att.Absent_Days7,
        att.Present_Days7, to much missing data to be useful -- check with Jieun */
        sf.fidelity_complete_percent,
        sf.fidelity_started_sum,
        sf.fidelity_complete_sum,
        ast.post_total_math_score,
        ast.Scale_Score5,
        ast.Performance_Level5,
        ast.pre_total_math_score,
        ast.pre_percentage_math_score,
        ast.pre_sub_P_score,
        ast.pre_sub_C_score,
        ast.pre_sub_F_score,
        ast.pre_math_completed_num,
        ast.pre_math_completed_percent,
        ast.pre_total_time_on_tasks,
        ast.pre_avg_time_on_tasks,
        ast.pre_MA_total_score,
        ast.pre_MA_avg_score,
        ast.pre_MSE_total_score,
        ast.pre_MSE_avg_score,
        ast.pre_PS_tasks_total_score,
        ast.pre_PS_part1_score,
        ast.pre_PS_part2E_score,
        ast.pre_PS_part2NE_score,
        ast.pre_PS_completed_num,
        ast.pre_PS_completed_percent,
        ast.pre_PS_total_RT_sec,
        ast.pre_PS_total_RT_min,
        ast.pre_PS_avg_RT_sec,
        ast.Scale_Score7,
        ast.Performance_Level7,
        ast.pre_total_math_score,
        assist.*
        
 from crosswalk_student cws
     left join student_demo  sd on sd.StuID = cws.StuID
     left join student_attendance att on att.StuID = cws.STUDID
     left join student_fidelity sf on sf.StuID = cws.StuID
     left join assess_student ast on ast.StuID = cws.StuID
     left join assist_student assist on assist.StuID = cws.StuID
     inner join student_roster sr on sr.StuID = cws.StuID 
            and sr.DROPSCH1 = 0 
            and sr.DROPSCH2 = 0 
            and sr.rdm_condition in ('BAU', 'ASSISTments')
 where cws.assistments_user_id is not null;

                             ")
colnames(dat)

write.csv(dat, "02_data/data_uptated_8_18.csv")
### Clean Data ####

# drop dup columns (from my lazy query)
unique_names <- unique(colnames(dat))

# Keep Only Unique Column Names
dat<-dat[unique_names]

table(dat$condition_assignment)

# DROP ALL POST TEST SCORES 
dat <- dat %>%
  filter(!is.na(dat$post_total_math_score))

# no cross contamination
table(dat$condition_assignment, 
      dat$num_problem_parts_used_bottom_out_hint > 0)

table(dat$condition_assignment, 
      dat$total_hints_accessed > 0)

# non of the students with missing assignments data did more than 3 assignments (most likely the assessments)
table(dat$fidelity_started_sum , is.na(dat$num_assignments_started))

# if NA for ASSISTments problem data than zero
colnames(dat)
dat[44:67] <- lapply(dat[44:67], ifnull0)
table(is.na(dat$num_assignments_started))

# # group race/ethnicity
dat <- dat %>%
  mutate(raceEthnicityFed=raceEthnicityFed %>%
                  factor()%>%
                  fct_lump_min(100),
                  fct_recode(`Hispanic/Latino`="1",Asian="3",White="6")) 
      
table(dat$raceEthnicityFed)



#### Missing Data ######

colSums(is.na(dat))
  # there are 58 students who have no ASSITments data --> Non compliers
colnames(dat)
lapply(dat, class)
dat[(1:13)] <- lapply(lapply(dat[1:13], as.character), as.factor)
lapply(dat[1:15], class)
dat[19] <- lapply(lapply(dat[19], as.character), as.factor)
dat[(14:18)] <- lapply(dat[14:18], as.numeric)
dat[(20:42)] <- lapply(dat[20:42], as.numeric)
dat[43] <- lapply(lapply(dat[43], as.character), as.factor)
dat[(44:67)] <- lapply(dat[44:67], as.numeric)
lapply(dat, n_distinct)

# add missing values count flag
dat$count_na <- rowSums(is.na(dat))
hist(dat$count_na)


imp <- dat%>%
  select(-(1:3), 
       #  -post_total_math_score, # should we be imputing the post test scores
         -(43:68)) %>%
  missForest( 
    variablewise = T
    )

# what cont variable has less than 5 levels?

length(imp$ximp)
names(imp$ximp)
imp_eval_metrics<- as.data.frame(cbind(lapply(imp$ximp, class), names(imp$OOBerror), imp$OOBerror))
imp_eval_metrics


dat_imp <-cbind(dat[1:3],  
                imp$ximp, 
                dat[43:68])

colSums(is.na(dat_imp))
write.csv(dat_imp, "02_data/data_uptated_8_18_IMPUTSED.csv")

#### Instant Data ####

table(dat$condition_assignment)
table(dat_imp$condition_assignment)

ps_Instant <-dat_imp %>%
  filter(condition_assignment == 'Instant')



#### Bottom-Outer Cut Point ####

table( ps_Instant$fidelity_started_sum < ps_Instant$num_assignments_started) # this shouldnt happen --> weird

# 


# cut point based on raw numbers
quantile(ps_Instant$num_problem_parts_used_bottom_out_hint)

quantile(ps_Instant[ps_Instant$num_problems_started > 0, ]$num_problem_parts_used_bottom_out_hint, na.rm = T)

ps_Instant %>%
  ggplot(aes(x=num_problem_parts_used_bottom_out_hint)) +
  geom_histogram() +
  geom_vline(xintercept =11, color = "red") +
  geom_vline(xintercept =28.5, color = "blue") +
  theme_minimal()


# MEDIAN SPLIT Raw definition of bottom-outer: students use used bottom out hints for > 8% of the problem parts completed
ps_Instant$bottom_outer_raw <- ifelse(ps_Instant$num_problem_parts_used_bottom_out_hint > 11, 1, 0)
table(ps_Instant$bottom_outer_raw)
table(ps_Instant[ps_Instant$num_problems_started > 0, ]$bottom_outer_raw)

# cut point based on % of problems where bottom out hint was used
ps_Instant$per_problems_parts_used_bottom_out_hint <- ifelse(ps_Instant$num_graded_problem_parts_attempted == 0, 0,
                                                             ps_Instant$num_problem_parts_used_bottom_out_hint/ps_Instant$num_graded_problem_parts_attempted)

hist(ps_Instant$per_problems_parts_used_bottom_out_hint, breaks = 50)
cor.test(ps_Instant$num_graded_problem_parts_attempted, ps_Instant$per_problems_parts_used_bottom_out_hint) # not correlated with numer of problem parts attempted
cor.test(ps_Instant$num_graded_problem_parts_attempted, ps_Instant$num_problem_parts_used_bottom_out_hint) # not correlated with numer of problem parts attempted
psych::describe(ps_Instant$per_problems_parts_used_bottom_out_hint)
quantile(ps_Instant$per_problems_parts_used_bottom_out_hint)
table(is.na(ps_Instant$per_graded_problem_parts_attempted))
ps_Instant %>%
  ggplot(aes(x=per_problems_parts_used_bottom_out_hint)) +
  geom_histogram() +
  geom_vline(xintercept = .08, color = "red") +
  geom_vline(xintercept = .17, color = "blue") +
  theme_minimal()

# MEDIAN SPLIT % definition of bottom-outer: students use used bottom out hints for > 8% of the problem parts completed
ps_Instant$bottom_outer_pre <- ifelse(ps_Instant$per_problems_parts_used_bottom_out_hint > .08, 1, 0)

### compare 2 definitions
ps_Instant %>%
  ggplot(aes(x=num_problem_parts_used_bottom_out_hint, y =num_graded_problem_parts_attempted,
             color = as.factor(bottom_outer_pre),
             shape = as.factor(bottom_outer_raw))) +
  geom_point() +
  geom_vline(xintercept = 11, color = "red") +
  theme_minimal()

table(ps_Instant$bottom_outer_pre, ps_Instant$bottom_outer_raw)


ps_Instant %>%
  select(
    starts_with("pre")
  ) %>%
  cor %>%
  corrplot::corrplot( method = 'number',
                      type = 'lower'
                      ) 
### teacher level variance
ps_Instant %>%
  dplyr::group_by(TeaIDEnd) %>%
  dplyr::summarise(
    n = n(),
    n_bottom = sum(bottom_outer_raw),
    p_bottom = round(n_bottom/ n, 2)
  )  %>%
  ggplot(aes(x= TeaIDEnd, y = p_bottom, fill = n)) +
  geom_bar(stat="identity")

### school level variance
ps_Instant %>%
  dplyr::group_by(SchIDEnd) %>%
  dplyr::summarise(
    n = n(),
    n_bottom = sum(bottom_outer_raw),
    p_bottom = round(n_bottom/ n, 2)
  )  %>%
  ggplot(aes(x= SchIDEnd, y = p_bottom, fill = n)) +
  geom_bar(stat="identity")

#### PS Models based on DEMO and PRETEST ####

#### Testing and Train data sets #####
set.seed(50513)
train <- ps_Instant %>%
  sample_frac(size = .75)


test <- ps_Instant %>%
  filter(!StuID  %in% train$StuID ) 

# model with all pretest data
  ps_m_all <- glmer(
    bottom_outer_raw ~
      virtual +
      raceEthnicityFed +
      EIP +
      ESOL +
      GIFTED +
      pre_total_math_score +
      pre_percentage_math_score +
      pre_sub_P_score +
      pre_sub_C_score +
      pre_sub_F_score +
      pre_math_completed_percent +
      scale(pre_total_time_on_tasks) +
      pre_MA_total_score +
      pre_MSE_total_score +
      pre_PS_tasks_total_score +
      pre_PS_part1_score +
      pre_PS_part2E_score +
      pre_PS_part2NE_score +
      pre_PS_completed_num +
      scale(pre_PS_total_RT_min) +
      (1 | TeaIDEnd) + (1|SchIDEnd)
    ,
    data = train,
    family = "binomial"
  )
summary(ps_m_all)

summary(ps_m_all)
pROC::auc(test$bottom_outer_raw, predict(ps_m_all, test, type = "response",  allow.new.levels = T))

test$bottom_outer_ALL <-  predict(ps_m_all, test, type = "response",  allow.new.levels = T)
cpALL<-cutpointr(test,
                 bottom_outer_ALL, bottom_outer_raw,
               method = maximize_metric, metric = sum_sens_spec
)
summary(cpALL)

#### P of Bottom_outer Models ######
  
ps_m1 <- glmer(
                 bottom_outer_raw ~
                   virtual +
                 pre_percentage_math_score +
                 scale(pre_total_time_on_tasks) +
                   Scale_Score5*(pre_MA_total_score +
                 pre_MSE_total_score) +
                 pre_PS_tasks_total_score +
                 scale(pre_PS_total_RT_min) +
                 Scale_Score5 +
                 (1 | TeaIDEnd),
  data = train,
  family = binomial()
)
summary(ps_m1)
pROC::auc(test$bottom_outer_raw, predict(ps_m1, test, type = "response",  allow.new.levels = T))

test$bottom_outer_p1 <-  predict(ps_m1, test, type = "response",  allow.new.levels = T)
cp1<-cutpointr(test,
   bottom_outer_p1, bottom_outer_pre,
  method = maximize_metric, metric = sum_sens_spec
)
summary(cp1)
plot(cp1)

binnedplot(predict(ps_m1, train, type = "response",  allow.new.levels = T), 
          resid(ps_m1))

binnedplot(train$Scale_Score5, 
           resid(ps_m1),
           xlab = "Scale Score G5")

binnedplot(train$num_assignments_started, 
           resid(ps_m1),
           xlab = "Number of Assignment Started")


binnedplot(train$fidelity_started_sum, 
           resid(ps_m1),
           xlab = "Number of Assignment Started (Including Assesments)",
           nclass = 25)

par(mfrow=c(2,1))
binnedplot(train$count_na, 
           resid(ps_m1),
           xlab = "Number of Imputed Variables",
           xlim = c(0, 24),
           nclass = 25
           )
boxplot(train$count_na, horizontal = T, ylim = c(0, 24),
        xlab = "Number of Imputed Variables")

### models with in program data which are exists for both instant and delayed students
colnames(ps_Instant)

table(dat_imp$condition_assignment, dat_imp$fidelity_started_sum)
t.test(dat_imp$fidelity_started_sum ~ dat_imp$condition_assignment )

describeBy(dat$avg_accuracy_first_attempt_before_hint, dat_imp$condition_assignment)
t.test(dat$avg_accuracy_first_attempt_before_hint ~ dat_imp$condition_assignment) # instant students did better on next problems correctness

describeBy(dat$num_graded_problem_parts_attempted, dat_imp$condition_assignment)
t.test(dat$num_graded_problem_parts_attempted ~ dat_imp$condition_assignment)# instant students did fewer problems




ps_m3 <- glmer(
    bottom_outer_raw ~
    virtual +
    pre_percentage_math_score +
    scale(pre_total_time_on_tasks) +
    Scale_Score5*(pre_MA_total_score +
                    pre_MSE_total_score) +
    pre_PS_tasks_total_score +
    scale(pre_PS_total_RT_min) +
    fidelity_started_sum +
    
    (1 | TeaIDEnd),
  data = train,
  family = binomial()
)
summary(ps_m3)
pROC::auc(ps_Instant$bottom_outer_raw, predict(ps_m3, ps_Instant, type = "response",  allow.new.levels = T))

ps_Instant$bottom_outer_p3 <-  predict(ps_m3, ps_Instant, type = "response",  allow.new.levels = T)
cp3<-cutpointr(ps_Instant,
               bottom_outer_p3, bottom_outer_raw,
               method = maximize_metric, metric = sum_sens_spec,
               na.rm = T
)
summary(cp3)
plot(cp3)


binnedplot(predict(ps_m3, train, type = "response",  allow.new.levels = T), 
           resid(ps_m3))

binnedplot(train$Scale_Score7, 
           resid(ps_m3),
           xlab = "Scale Score G5")


binnedplot(train$fidelity_started_sum, 
           resid(ps_m3),
           xlab = "Number of Assignment Started (Including Assesments)",
           nclass = 25)

par(mfrow=c(2,1))
binnedplot(train$count_na, 
           resid(ps_m3),
           xlab = "Number of Imputed Variables",
           xlim = c(0, 24),
           nclass = 25
)
boxplot(train$count_na, horizontal = T, ylim = c(0, 24),
        xlab = "Number of Imputed Variables")
par(mfrow=c(1,1))




#### compare models with and without fidelity metric

par(mfrow=c(2,1))
binnedplot(train$Scale_Score7, 
           resid(ps_m1,
                 nclass = 25),
           xlab = "Scale Score G5")
binnedplot(train$Scale_Score7, 
           resid(ps_m3,
                 nclass = 25),
           xlab = "Scale Score G5")
par(mfrow=c(1,1))

par(mfrow=c(2,1))
binnedplot(train$fidelity_started_sum, 
           resid(ps_m1),
           xlab = "Scale Score G5",
           nclass = 25)
binnedplot(train$fidelity_started_sum, 
           resid(ps_m3),
           xlab = "Scale Score G5",
           nclass = 25)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
binnedplot(train$count_na, 
           resid(ps_m1),
           xlab = "Scale Score G5",
           nclass = 25)
binnedplot(train$count_na, 
           resid(ps_m3),
           xlab = "Scale Score G5",
           nclass = 25)
par(mfrow=c(1,1))
