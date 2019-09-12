#Written by Deus Thindwa
#Does living in a household with HIV-infected individuals increase child risk of pneumococcal carriage acquisition?
#An analysis using hidden markov transition models, PhD chapter 1.
#13/08/2019

#===============load required packages into memory===============
phirst.packages <-c("tidyverse","plyr","msm","dplyr","minqa","lubridate", "magrittr","data.table")
lapply(phirst.packages, library, character.only=TRUE)

#===============simulate a phirst dataset===============
set.seed(12345)
sample.size=45994

#8 individuals per household (320 total obs, 8 HH members, each with 40 obs)
phirst1 <- data.frame(hid=rep(letters[1:26], each=320), pid=rep.int(1:8, 40),prefix="a")
phirst1$hhid <- apply(phirst1, 1 , function(x) paste0(toString(x[1]), toString(x[3])))
phirst1$iid <- apply(phirst1, 1 , function(x) paste0(toString(x[4]), toString(x[2])))
phirst1 <- subset(phirst1, select=c(hhid, iid))
phirst1 <- ddply(phirst1, .(iid), mutate, sample.p = seq(0,by=4,length.out=40), state=as.integer(sample(c("1","2"), size=40, replace=TRUE)), age=rep(abs(rnorm(n=sample.size,mean=20,sd=9.5)),each=40,length.out=40), agecat=ifelse(age<5,1,2))
phirst1$hiv <- if_else(phirst1$age<1,"1",rep(sample(c("1","2"),size=8320, replace=TRUE), each=40, length.out=8320))

#4 individuals per household (240 total obs, 4 HH members, each with 60 obs)
phirst2 <- data.frame(hid=rep(letters[1:26], each=240), pid=rep.int(1:4, 60),prefix="b")
phirst2$hhid <- apply(phirst2, 1 , function(x) paste0(toString(x[1]), toString(x[3])))
phirst2$iid <- apply(phirst2, 1 , function(x) paste0(toString(x[4]), toString(x[2])))
phirst2 <- subset(phirst2, select=c(hhid, iid))
phirst2 <- ddply(phirst2, .(iid), mutate, sample.p = seq(0,by=4,length.out=60), state=as.integer(sample(c("1","2"), size=60, replace=TRUE)), age=rep(abs(rnorm(n=sample.size,mean=20,sd=9.5)),each=60,length.out=60), agecat=ifelse(age<5,1,2))
phirst2$hiv <- if_else(phirst2$age<1,"1",rep(sample(c("1","2"),size=6240, replace=TRUE), each=60, length.out=6240))

#12 individuals per household (384 total obs, 12 HH members, each with 32 obs)
phirst3 <- data.frame(hid=rep(letters[1:26], each=384), pid=rep.int(1:12, 32),prefix="d")
phirst3$hhid <- apply(phirst3, 1 , function(x) paste0(toString(x[1]), toString(x[3])))
phirst3$iid <- apply(phirst3, 1 , function(x) paste0(toString(x[4]), toString(x[2])))
phirst3 <- subset(phirst3, select=c(hhid, iid))
phirst3 <- ddply(phirst3, .(iid), mutate, sample.p=seq(0,by=4,length.out=32), state=as.integer(sample(c("1","2"),size=32,replace=TRUE)), age=rep(abs(rnorm(n=sample.size,mean=20,sd=9.5)),each=32,length.out=32), agecat=ifelse(age<5,1,2))
phirst3$hiv <- if_else(phirst3$age<1,"1",rep(sample(c("1","2"),size=9984, replace=TRUE), each=32, length.out=9984))

#15 individuals per household (825 total obs, 15 HH members, each with 55 obs)
phirst4 <- data.frame(hid=rep(letters[1:26], each=825), pid=rep.int(1:15, 55),prefix="e")
phirst4$hhid <- apply(phirst4, 1 , function(x) paste0(toString(x[1]), toString(x[3])))
phirst4$iid <- apply(phirst4, 1 , function(x) paste0(toString(x[4]), toString(x[2])))
phirst4 <- subset(phirst4, select=c(hhid, iid))
phirst4 <- ddply(phirst4, .(iid), mutate, sample.p=seq(0,by=4,length.out=55), state=as.integer(sample(c("1","2"),size=55,replace=TRUE)), age=rep(abs(rnorm(n=sample.size,mean=20,sd=9.5)),each=55,length.out=55), agecat=ifelse(age<5,1,2))
phirst4$hiv <- if_else(phirst4$age<1,"1",rep(sample(c("1","2"),size=21450, replace=TRUE), each=55, length.out=21450))

#combine all datasets
phirst <- rbind(phirst1, phirst2, phirst3,phirst4)
remove(phirst1,phirst2,phirst3,phirst4)

#add missing state with 2.2% missing swabs
phirst$state[sample(1:length(phirst$state),2000)]<-NA
phirst$state_m <-phirst$state
phirst$state_m[is.na(phirst$state)] <-3

#add new serotype state for subset of available PCV13 serotypes
phirst$serotype <- if_else(phirst$state==2,sample(c("1","3","4","5","6A","6B","7F","9V","14","18C","19A","19F","23F"),size=sample.size,replace=TRUE),NULL) #serotype
phirst$serotype[sample(1:length(phirst$serotype),2000)]<-NA
phirst$ST <- phirst$serotype
phirst <- phirst[order(phirst$iid,phirst$sample.p),]

equal_first <- function(x) {x %>% equals(x[1]) %>% not %>% as.numeric}
phirst %>% group_by(iid) %>% mutate_each(funs(equal_first),starts_with("serotype"))


#===============hidden markov models using "msm" R package===============

#---------------MODEL1 (FALSE NEGATIVES)---------------

#show the state table
statetable.msm(state,iid,data=phirst)
                   
#define initial values of a transition intensity matrix Q
model1.Q <- rbind(c(0.0,0.1), 
                  c(0.1,0.0))
rownames(model1.Q) <- c("Clear","Carrying")
colnames(model1.Q) <- c("Clear","Carrying")

#define initial values of an emission matrix E
model1.E <- rbind(c(0.0,0.0), 
                  c(0.1,0.0))
colnames(model1.E) <- c("SwabNeg","SwabPos")
rownames(model1.E) <- c("Clear","Carrying")

#create an outcome distribution of each observed state (swabs)
model1.HMM <- list(hmmBinom(1,0.2), hmmBinom(2,0.8))

#fit a 2-state hidden markov model accounting for false negatives
model1.msm <- msm(state~sample.p, subject=iid, data=phirst,
                  qmatrix=model1.Q,
                  ematrix=model1.E,
                  hmodel=model1.HMM,
                  covariates= ~agecat,
                  est.initprobs=T,
                  opt.method="bobyqa")

#maximum likelihood estimates of intensity matrix and outcome distributions
printold.msm(model1.msm)

#transition rates and 95%CI of HMM
sojourn.msm(model1.msm)

#most likely true series of states underlying the data
viterbi.msm(model1.msm)

#diagnosdtics plots
plot.survfit.msm(model1.msm, main="model 1", mark.time=FALSE)
plot.prevalence.msm(model1.msm)

#---------------MODEL3 (FALSE NEGATIVES + MISSING SWABS)---------------

#define initial values of a transition intensity matrix Q
model3.Q <- rbind(c(0.0,0.1,0.0),
                  c(0.1,0.0,0.0),
                  c(0.0,0.0,0.0))
rownames(model3.Q) <- c("Clear","Carrying","Missing")
colnames(model3.Q) <- c("Clear","Carrying","Missing")

#define initial values of an emission matrix E
model3.E <- rbind(c(0.0,0.0,0.1), 
                  c(0.1,0.0,0.1), 
                  c(0.0,0.0,0.0))
rownames(model3.E) <- c("Clear","Carrying","Missing")
colnames(model3.E) <- c("SwabNeg","SwabPos","Missing")

#create an outcome distribution of each observed state (swabs)
model3.HMM <- list(hmmCat(0.2,1),hmmCat(0.7,1),hmmIdent())

#fit a 2-state hidden markov model accounting for false negatives
model3.msm <- msm(state_m~sample.p, subject=iid, data=phirst,
                  qmatrix=model3.Q,
                  ematrix=model3.E,
                  hmodel=model3.HMM,
                  est.initprobs=T,
                  opt.method="bobyqa")

#maximum likelihood estimates of intensity matrix and outcome distributions
printold.msm(model3.msm)

#transition rates and 95%CI of HMM
sojourn.msm(model3.msm)

#most likely true series of states underlying the data
viterbi.msm(model3.msm)

#diagnosdtics plots
plot.survfit.msm(model3.msm, main="model 3", mark.time=FALSE)
plot.prevalence.msm(model3.msm)

#---------------MODEL2 (FALSE NEGATIVES + MULTIPLE ACQUISITIONS)---------------







