## R code for LCMM and multinomial logistic analysis

library(lcmm)
library(reshape)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stargazer)
library(stargazer)

# PTSD SYMPTOM TRAJECTORIES VIA LCMM ####
### 1) ED DISCHARGE --- iSTAR ####
# "iSTAR" IS DATA FROM PARENT STUDY OF HOSPITALIZED PARTICIPANTS
# FILTER TO RETAIN ONLY BLACK PARTICIPANTS
set.seed(87)
iSTAR_PCL_timepoints<-iSTAR %>% filter(Racial_Category == 4) %>% select(record_id, PCL5_Total.2WkDay2a, PCL5_Total.3Mo, PCL5_Total.6Mo)

# REMOVE SUBS MISSING 2 OR MORE (OF 3) TIME POINTS
iSTARmissingPCL<-rowSums(is.na(iSTAR_PCL_timepoints))
iSTAR_PCL_timepoints<-iSTAR_PCL_timepoints[-which(iSTARmissingPCL>=2),]

# RENAME COLUMNS TO SOMETHING READABLE
colnames(iSTAR_PCL_timepoints)<-c("ID", "T1", "T2", "T3")
iSTAR_PCL_timepoints<-as.data.frame(iSTAR_PCL_timepoints)
# RESHAPE DATA TO LONG FORMAT (by PCL at the 3 time points)
iSTAR_reshape_PCL_timepoints<- melt(data=iSTAR_PCL_timepoints, id=c("ID"), measures=c(colnames(iSTAR_PCL_timepoints[,2:ncol(iSTAR_PCL_timepoints)])))

# RENAME COLUMNS TO SOMETHING READABLE
colnames(iSTAR_reshape_PCL_timepoints)<-c("ID", "TIME", "PCL")
iSTAR_reshape_PCL_timepoints$TIME<-as.numeric(iSTAR_reshape_PCL_timepoints$TIME)

# CALCULATE LCMM MODELS FOR N=1-6 CLASS SOLUTIONS
# LINEAR TERM ONLY
d1_istar3<-lcmm(PCL~TIME, random=~TIME, subject="ID", ng=1, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear")
d2_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=2, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d3_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=3, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d4_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=4, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d5_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=5, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d6_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=6, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))

# COMPARE MODEL FIT - LINEAR TERM ONLY
summarytable(d1_istar3, d2_istar3, d3_istar3, d4_istar3, d5_istar3, d6_istar3, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

# TRY QUADRATIC TERM
# CALCULATE QUADRATIC TERM
iSTAR_reshape_PCL_timepoints$TIME2 <- as.numeric(iSTAR_reshape_PCL_timepoints$TIME)^2
# GRID SEARCH EMPLOYED FOR ng=2-6, STARTING FROM INITIAL VALUES FROM ng=1 SOLUTION ABOVE
# 50 REPETITIONS (DEPARTURES FROM RANDOM INTIAL VALUES)
# 100 ITERATIONS TO OPTIMIZE ALGORITHM
d1_istar3_2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", ng=1, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear")
d2_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=2, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d3_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=3, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d4_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=4, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d5_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=5, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d6_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=6, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))

# COMPARE MODEL FIT - QUADRATIC TERM
#summarytable(d1_istar3_2, d2_istar3_2, d3_istar3_2, d4_istar3_2, d5_istar3_2, d6_istar3_2, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

# "iSTAR" DATASET IS DATA FROM PARENT STUDY OF HOSPITALIZED PARTICIPANTS
# FILTER TO RETAIN ONLY BLACK PARTICIPANTS 
iSTAR_PCL_timepoints<-iSTAR %>% filter(Racial_Category == 4) %>% select(record_id, PCL5_Total.2WkDay2a, PCL5_Total.3Mo, PCL5_Total.6Mo) 

# REMOVE SUBS MISSING 2 OR MORE (OF 3) TIME POINTS
iSTARmissingPCL<-rowSums(is.na(iSTAR_PCL_timepoints))
iSTAR_PCL_timepoints<-iSTAR_PCL_timepoints[-which(iSTARmissingPCL>=2),]

# RENAME COLUMNS TO SOMETHING READABLE
colnames(iSTAR_PCL_timepoints)<-c("ID", "T1", "T2", "T3")
iSTAR_PCL_timepoints<-as.data.frame(iSTAR_PCL_timepoints)
# RESHAPE DATA TO LONG FORMAT (by PCL at the 3 time points)
iSTAR_reshape_PCL_timepoints<- melt(data=iSTAR_PCL_timepoints, id=c("ID"), measures=c(colnames(iSTAR_PCL_timepoints[,2:ncol(iSTAR_PCL_timepoints)])))

# RENAME COLUMNS TO SOMETHING READABLE
colnames(iSTAR_reshape_PCL_timepoints)<-c("ID", "TIME", "PCL")
iSTAR_reshape_PCL_timepoints$TIME<-as.numeric(iSTAR_reshape_PCL_timepoints$TIME)

# CALCULATE LCMM MODELS FOR N=1-6 CLASS SOLUTIONS
# LINEAR TERM ONLY
d1_istar3<-lcmm(PCL~TIME, random=~TIME, subject="ID", ng=1, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear")
d2_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=2, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d3_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=3, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d4_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=4, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d5_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=5, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d6_istar3<-gridsearch(minit=d1_istar3, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=6, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))


# COMPARE MODEL FIT - LINEAR TERM ONLY
summarytable(d1_istar3, d2_istar3, d3_istar3, d4_istar3, d5_istar3, d6_istar3, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

# TRY QUADRATIC TERM
# CALCULATE QUADRATIC TERM
iSTAR_reshape_PCL_timepoints$TIME2 <- as.numeric(iSTAR_reshape_PCL_timepoints$TIME)^2
# GRID SEARCH EMPLOYED FOR ng=2-6, STARTING FROM INITIAL VALUES FROM ng=1 SOLUTION ABOVE
# 50 REPETITIONS (DEPARTURES FROM RANDOM INTIAL VALUES)
# 100 ITERATIONS TO OPTIMIZE ALGORITHM
d1_istar3_2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", ng=1, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear")
d2_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=2, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d3_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=3, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d4_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=4, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d5_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=5, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))
d6_istar3_2<-gridsearch(minit=d1_istar3_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=6, idiag=FALSE, data=iSTAR_reshape_PCL_timepoints, link="linear"))

# COMPARE MODEL FIT - QUADRATIC TERM
summarytable(d1_istar3_2, d2_istar3_2, d3_istar3_2, d4_istar3_2, d5_istar3_2, d6_istar3_2, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))


# FOR BEST FITTING MODEL GRAB POSTERIOR PROBABILITIES
people4_istar3 <- as.data.frame(d4_istar3$pprob[,1:2])
# PUT CLASS ASSIGNMENTS FROM POSTERIOR PROBABILITIES INTO OUR DATAFRAME WITH THE PCL SCORES SO WE CAN PLOT CLASSES
iSTAR_reshape_PCL_timepoints$group4 <- factor(people4_istar3$class[sapply(iSTAR_reshape_PCL_timepoints$ID, function(x) which(people4_istar3$ID==x))])

iSTAR_reshape_PCL_timepoints$group4<-factor(iSTAR_reshape_PCL_timepoints$group4, labels = c("Non-remitting","Moderate",  "Remitting", "Resilient"), levels = c(3,1,4,2))


### 2) HOSPITALIZED --- STAR 2.0 ####
# "STAR2" DATASET IS DATA FROM PARENT STUDY OF HOSPITALIZED PARTICIPANTS
# FILTER TO RETAIN ONLY BLACK PARTICIPANTS 
set.seed(87)
STAR2_PCL_timepoints_black<-STAR2 %>% filter(Racial_Category.bl==3) %>% select(record_id, PCL5_Total.bl, PCL5_Total.1mo, PCL5_Total.6mo) 
STAR2_PCL_timepoints_black<-as.data.frame(STAR2_PCL_timepoints_black)
# REMOVE SUBS MISSING 2 OR MORE (OF 3) TIME POINTS
STAR2missingPCL_black<-rowSums(is.na(STAR2_PCL_timepoints_black))
STAR2_PCL_timepoints_black<-STAR2_PCL_timepoints_black[-which(STAR2missingPCL_black>=2),]

# RENAME COLUMNS TO SOMETHING READABLE
colnames(STAR2_PCL_timepoints_black)<-c("ID", "T1", "T2", "T3")

# RESHAPE DATA TO LONG FORMAT (by PCL at the 3 time points)
STAR2_reshape_PCL_timepoints_black<- melt(data=STAR2_PCL_timepoints_black, id=c("ID"), measures=c(colnames(STAR2_PCL_timepoints_black[,2:ncol(STAR2_PCL_timepoints_black)])))

# RENAME COLUMNS TO SOMETHING READABLE
colnames(STAR2_reshape_PCL_timepoints_black)<-c("ID", "TIME", "PCL")
STAR2_reshape_PCL_timepoints_black$TIME<-as.numeric(STAR2_reshape_PCL_timepoints_black$TIME)

# CALCULATE LCMM MODELS FOR N=1-6 CLASS SOLUTIONS
# LINEAR TERM ONLY
d1_star2_black<-lcmm(PCL~TIME, random=~TIME, subject="ID", ng=1, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear")
d2_star2_black<-gridsearch(minit=d1_star2_black, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=2, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear"))
d3_star2_black<-gridsearch(minit=d1_star2_black, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=3, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear"))
d4_star2_black<-gridsearch(minit=d1_star2_black, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=4, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear"))
d5_star2_black<-gridsearch(minit=d1_star2_black, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=5, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear"))
d6_star2_black<-gridsearch(minit=d1_star2_black, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=6, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear"))

# COMPARE MODEL FIT - LINEAR TERM ONLY
summarytable(d1_star2_black, d2_star2_black, d3_star2_black, d4_star2_black, d5_star2_black, d6_star2_black, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

# TRY QUADRATIC TERM
# CALCULATE QUADRATIC TERM
STAR2_reshape_PCL_timepoints_black$TIME2 <- as.numeric(STAR2_reshape_PCL_timepoints_black$TIME)^2
# GRID SEARCH EMPLOYED FOR ng=2-6, STARTING FROM INITIAL VALUES FROM ng=1 SOLUTION ABOVE
# 50 REPETITIONS (DEPARTURES FROM RANDOM INTIAL VALUES)
# 100 ITERATIONS TO OPTIMIZE ALGORITHM
d1_star2_2_black<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", ng=1, idiag=FALSE, data=STAR2_reshape_PCL_timepoints_black, link="linear")
d2_star2_2<-gridsearch(minit=d1_star2_2_black, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=2, idiag=FALSE, data=STAR2_reshape_PCL_timepoints, link="linear"))
d3_star2_2<-gridsearch(minit=d1_star2_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=3, idiag=FALSE, data=STAR2_reshape_PCL_timepoints, link="linear"))
d4_star2_2<-gridsearch(minit=d1_star2_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=4, idiag=FALSE, data=STAR2_reshape_PCL_timepoints, link="linear"))
d5_star2_2<-gridsearch(minit=d1_star2_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=5, idiag=FALSE, data=STAR2_reshape_PCL_timepoints, link="linear"))
d6_star2_2<-gridsearch(minit=d1_star2_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=6, idiag=FALSE, data=STAR2_reshape_PCL_timepoints, link="linear"))

# COMPARE MODEL FIT - QUADRATIC TERM
summarytable(d1_istar3_2, d2_istar3_2, d3_istar3_2, d4_istar3_2, d5_istar3_2, d6_istar3_2, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))


# FOR EACH MODEL GRAB POSTERIOR PROBABILITIES
people3_star2_black <- as.data.frame(d3_star2_black$pprob[,1:2])

# PUT CLASS ASSIGNMENTS FROM POSTERIOR PROBABILITIES INTO OUR DATAFRAME WITH THE PCL SCORES SO WE CAN PLOT CLASSES
people3_star2_black$ID<-as.character(people3_star2_black$ID)

STAR2_reshape_PCL_timepoints_black<-
  STAR2_reshape_PCL_timepoints_black %>% 
  left_join(people3_star2_black, by = c("ID" = "ID")) %>% 
  rename("group3"="class") %>% 
  recode(STAR2_reshape_PCL_timepoints_black$group3, `1`="Delayed", `2`="Resilient", `3`="Non-remitting")


#
#
#
# FIGURE 1 --- PLOT LCMM SOLUTIONS
#
#
#
plot_istar_trajectory<- 
  ggplot(iSTAR_reshape_PCL_timepoints, aes(TIME, PCL, group=ID, colour=group4)) + 
  geom_line(data=na.omit(iSTAR_reshape_PCL_timepoints)) + 
  geom_smooth(aes(group=group4), method="loess", size=2, se=T) + 
  labs(title="", x="Time",y="PCL-5 Total Severity",colour="Trajectory") + 
  ggtitle("ED Discharged (n=111)") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=22,face="bold"), 
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        axis.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20))  + 
  geom_hline(yintercept=32, linetype="dashed", size=0.5) + 
  scale_x_continuous(breaks=c(1,2,3)) + 
  scale_x_continuous(breaks=c(1,2,3), labels=c("T1 (2-3 weeks)", "T2 (3 MO)", "T3 (6 MO)")) + 
  scale_color_manual(values=c("Non-remitting"="#F8766D", "Moderate"="#D39200", "Remitting"="#AC88FF",  "Resilient"="#00BFC4"))

plot_star2_trajectory<- 
  ggplot(STAR2_reshape_PCL_timepoints_black, aes(TIME, PCL, group=ID, colour=group3)) + 
  geom_line(data=na.omit(STAR2_reshape_PCL_timepoints_black)) + 
  geom_smooth(aes(group=group3), method="loess", size=2, se=T) + 
  labs(title="", x="Time",y="PCL-5 Total Severity",colour="Trajectory") + 
  ggtitle("Hospitalized (n=81)") + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size=22,face="bold"), 
        legend.title=element_text(size=20, face="bold"), 
        legend.text=element_text(size=20), 
        axis.title = element_text(size=20, face="bold"), 
        axis.text=element_text(size=20))  + 
  geom_hline(yintercept=32, linetype="dashed", size=0.5) + 
  scale_x_continuous(breaks=c(1,2,3)) + 
  scale_x_continuous(breaks=c(1,2,3), labels=c("T1 (hospital)", "T2 (3 MO)", "T3 (6 MO)")) + 
  scale_color_manual(values=c("Non-remitting"="#F8766D", "Delayed"="#7CAE00", "Resilient"="#00BFC4"))

# FIGURE 1 IN THE MANUSCRIPT
plot_trajectory_grid<-grid.arrange(plot_istar_trajectory, plot_star2_trajectory, nrow=1)


#
#
#
#
# DISCRIMINATION PREDICTS PTSD TRAJECTORY CLASS ####
#
#
#
#
#
#
### 1) ED DISCHARGE --- iSTAR ####
#
# FOLD TRAJECTORY ASSIGNMENTS BACK INTO FULL DATASET FOR REGRESSION ANALYSIS
iSTAR_traj<-
  iSTAR_reshape_PCL_timepoints %>% 
  left_join(iSTAR, by = c("ID" = "record_id"))

# RE-LEVEL TRAJECTORY FACTOR SO RESILIENT IS THE REFERENCE GROUP
iSTAR_traj$traj_relevel<-as.factor(iSTAR_lpa_traj$group4)
iSTAR_traj$traj_relevel<-relevel(iSTAR_lpa_traj$traj_relevel, ref="Resilient")

# MULTINOMIAL LOGISTIC REGRESSION
# NO COVARS
istar_multinom_nocovars<-multinom(traj_relevel~PEDQ_TOTAL, data=iSTAR_traj)
# TABLE OF RESULTS
stargazer(istar_multinom_nocovars, type='text', coef=list(exp(coef(istar_multinom_nocovars))), p.auto=FALSE)

# ALL COVARS
istar_multinom<-multinom(traj_relevel~PEDQ_TOTAL + age_at_time_of_data_collection.2WkDay1 + Gender + MOI +  LEC_Weighted_Total_2Wk + Education_relevel + time_since_injury.2WkDay2a + iss, data=iSTAR_traj)
# TABLE OF RESULTS
stargazer(istar_multinom, type='text', coef=list(exp(coef(istar_multinom))), p.auto=FALSE)


#
#
### 2) HOSPITALIZED --- STAR 2.0 ####
#
#

STAR2_traj<-
  STAR2_reshape_PCL_timepoints %>% 
  left_join(STAR2, by = c("ID" = "record_id"))
# RE-LEVEL TRAJECTORY FACTOR SO RESILIENT IS THE REFERENCE GROUP
STAR2_traj$traj_relevel<-as.factor(STAR2_lpa_traj$group3)
STAR2_traj$traj_relevel<-relevel(STAR2_lpa_traj$traj_relevel, ref="Resilient")


# MULTINOMIAL LOGISTIC REGRESSION
# NO COVARS
star2_multinom_nocovars<-multinom(traj_relevel~PEDQ_TOTAL, data=STAR2_traj)
# TABLE OF RESULTS
stargazer(star2_multinom_nocovars, type='text', coef=list(exp(coef(star2_multinom_nocovars))), p.auto=FALSE)

# ALL COVARS
star2_multinom<-multinom(traj_relevel~PEDQ_TOTAL + age.bl + gender.bl + MOI.bl + LEC_Weighted_Total.bl + education.bl + time_since_injury.bl + iss.bl, data=STAR2_traj)
# TABLE OF RESULTS
stargazer(star2_multinom, type='text', coef=list(exp(coef(star2_multinom))), p.auto=FALSE)


#
#
#
# FIGURE 2 --- PLOT PROBABILITY OF PTSD TRAJECTORY FROM DISCRIMINATION
#
#
#
# CALCULATE POSTERIOR PROBABILITIES FROM MULTINOM MODEL
pprob_istar_multinom<-ggeffects::ggeffect(istar_multinom, terms = c("PEDQ_TOTAL"))

plot_multinom_istar_covars<-
  ggplot(data=pprob_istar_multinom, aes(x=x, y=predicted, color=response.level, group=response.level)) + 
  geom_line(size=2) + 
  geom_point(size=3) + 
  theme_classic() + 
  labs(x="BPEDQ Total", 
       y="Probability", 
       title="ED Discharged (N=111)")  + 
  scale_color_manual(name="PTSD Trajectory", 
                     breaks=c("Nonremitting", "Moderate", "Remitting", "Resilient"), 
                     labels=c("Non-remitting", "Moderate", "Remitting", "Resilient"), 
                     values=c("#F8766D","#D39200","#AC88FF", "#00BFC4")) + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=14), 
        plot.title = element_text(hjust=.5, face="bold", size = 18)) 


# CALCULATE POSTERIOR PROBABILITIES FROM MULTINOM MODEL
pprob_star2_multinom<-ggeffects::ggeffect(star2_multinom, terms = "PEDQ_TOTAL")

plot_multinom_star2_covars<-
  ggplot(data=pprob_star2_multinom, aes(x=x, y=predicted, color=response.level, group=response.level)) + 
  geom_line(size=2) + 
  geom_point(size=3) + 
  theme_classic() + 
  labs(x="BPEDQ Total", y="Probability", title="Hospitalized (N=81)") + 
  scale_color_manual(name="PTSD Trajectory", 
                     breaks=c("Nonremitting", "Delayed", "Resilient"), 
                     labels=c("Non-remitting", "Delayed", "Resilient"), 
                     values=c("#F8766D", "#7CAE00", "#00BFC4")) + 
  theme(axis.text = element_text(size=14), 
        axis.title = element_text(size=16), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=14), 
        plot.title = element_text(hjust=.5, face="bold", size = 18)) 


# FIGURE 2 IN THE MANUSCRIPT
plot_pred_grid_covars<-grid.arrange(plot_multinom_istar_covars, plot_multinom_star2_covars, nrow=1)