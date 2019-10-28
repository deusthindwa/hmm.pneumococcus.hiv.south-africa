#Written by Deus Thindwa
#Are HIV+ adult PNC carriers more likely to transmit pneumococci than HIV- adult PNC carriers?
#Continuous-time time-inhomogeneous hidden Markov modelling study, PhD chapter 1.
#13/08/2019

#---------------load required packages into memory
phirst.packages <-c("tidyverse","plyr","msm","timetk","gridExtra","dplyr","minqa","lubridate", "magrittr","data.table","parallel")
lapply(phirst.packages, library, character.only=TRUE)

#---------------load a simulated phirst dataset
phirst <- cav[1:2846,]
phirst$PTNUM <- as.character(phirst$PTNUM)
phirst <- subset(phirst, select=c(PTNUM,years,state,age,pdiag,cumrej,firstobs))
phirst$age <- if_else(phirst$age<40 & phirst$firstobs==1,0L,if_else(phirst$age>=40 & phirst$firstobs==1,1L,NULL))
phirst <- phirst%>%group_by(PTNUM)%>%fill(starts_with("age"), .direction="down")
colnames(phirst) <- c("iid","wks","state","age","hiv","serotype","firstobs")
phirst$hiv <- if_else(phirst$hiv=="IHD",1L,0L)
phirst$state <- if_else(phirst$state==3,1L,if_else(phirst$state==4,4L,if_else(phirst$state==1,1L,2L)))
phirst$serotype <- if_else(phirst$state==2,sample(c("1","3","4","5","6A","6B","7F","9V","14","18C","19A","19F","23F"),size=2846,replace=TRUE),NULL) #serotype

#---------------show transition frequency in state table
statetable.msm(state,iid,data=phirst)

#---------------initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.1), 
                  c(0.1,0.0))
rownames(matrix.Q) <- c("Clear","Carrying")
colnames(matrix.Q) <- c("Clear","Carrying")

#---------------initiate emission matrix E
matrix.E <- rbind(c(0.0,0.0), 
                  c(0.1,0.0))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carrying")

#---------------fit a simpler model without misclassification
p.model1<-msm(state~wks, subject=iid, data=phirst, qmatrix=matrix.Q, censor=4, gen.inits=TRUE, covariates=~age+hiv, control=list(trace=1,REPORT=1), cl=0.95)
printnew.msm(p.model1)

#---------------assessment of models with and without covariates
#(a) observed vs expected
plot.prevalence.msm(p.model1, mintime=0, maxtime=15)
phirst.prev <-tk_tbl(prevalence.msm(p.model1), preserve_index = TRUE, rename_index = "time") 
phirst.prev <- subset(phirst.prev, select=c(time,Observed.percentages.State.1,Observed.percentages.State.2,Expected.percentages.Clear,Expected.percentages.Carrying))
colnames(phirst.prev) <- c("time","obs.clear","obs.carry","exp.clear","exp.carry")


A<-ggplot() +
  geom_line(data=phirst.prev, aes(x=time, y=obs.clear), colour="darkgreen") +
  geom_line(data=phirst.prev, aes(x=time, y=exp.clear), colour="darkgreen", lty="dashed") + 
  scale_y_discrete(limits = c(0,20,40,60,80,100)) + 
  scale_x_discrete(limits = c(0,5,10,15,20)) +
  coord_cartesian(ylim=c(0,100)) + 
  labs(title="Clear", x="time (weeks)", y="prevalence (%)") + 
  theme_classic()

B<-ggplot() +
  geom_line(data=phirst.prev, aes(x=time, y=obs.carry), colour="darkred") +
  geom_line(data=phirst.prev, aes(x=time, y=exp.carry), colour="darkred", lty="dashed") + 
  scale_y_discrete(limits = c(0,20,40,60,80,100)) + 
  scale_x_discrete(limits = c(0,5,10,15,20)) +
  coord_cartesian(ylim=c(0,100)) + 
  labs(title="Carrying", x="time (weeks)", y="prevalence (%)") + 
  theme_classic()

grid.arrange(grobs=list(A, B), ncol=2,nrow=2)

#(b) pearson-type good of fit
options(digits=2)
pearson.msm(p.model1, timegroups=2)

#---------------fit Markov model1 with covariates seperately
qmatrix.msm(p.model1, covariates=list(hiv=0))
qmatrix.msm(p.model1, covariates=list(hiv=1))
qmatrix.msm(p.model1, covariates="mean")

qmatrix.msm(p.model1, covariates=list(age=0))
qmatrix.msm(p.model1, covariates=list(age=1))
qmatrix.msm(p.model1, covariates="mean")

#---------------fit Markov model1 with covariates in combination
qmatrix.msm(p.model1, covariates=list(hiv=1, age=0))
qmatrix.msm(p.model1, covariates=list(hiv=1, age=1))
qmatrix.msm(p.model1, covariates=list(hiv=0, age=1))
qmatrix.msm(p.model1, covariates=list(hiv=0, age=0))



#create an outcome distribution of each observed state (swabs)
model1.HMM <- list(hmmBinom(1,0.5), hmmBinom(2,0.5))
model1.HMM <- list(hmmBinom(1,0.5), hmmBinom(1,0.5))

#fit a 2-state hidden markov model accounting for false negatives
model1.msm <- msm(state~sample.p, subject=iid, data=phirst,
                  qmatrix=model1.Q,
                  ematrix=model1.E,
                  hmodel=model1.HMM,
                  covariates= ~agecat+hiv,
                  pci=4,
                  est.initprobs=T,
                  opt.method="bobyqa")

#diagnosdtics plots
plot.survfit.msm(model3.msm, main="model 3", mark.time=FALSE)
plot.prevalence.msm(model3.msm)









