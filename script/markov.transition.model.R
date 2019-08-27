#Written by Deus Thindwa
#Does living in a howusehold with HIV-infected individuals increase child risk of pneumococcal carriage acquisition?
#An analysis using markov transition model, PhD chapter 1.
#13/08/2019

#===============load required packages into memory===============
phirst.packages <-c("tidyverse","plyr","msm","dplyr")
lapply(phirst.packages, library, character.only=TRUE)

#===============simulate a phirst dataset===============
set.seed(12345)
sample.size=4862

#8 individuals per household (32=total FU-frequency, 8=no. of HH members, 4=individual FU-frequency)
phirst1 <- data.frame(hid=rep(letters[1:26], each=32), pid=rep.int(1:8, 4))
phirst1$sample.date <- seq(as.Date('2018/01/01'), as.Date('2020/04/11'), by="day")
phirst1$prefix <- "a"
phirst1$hhid <- apply(phirst1, 1 , function(x) paste0(toString(x[4]), toString(x[1])))
phirst1$iid <- apply(phirst1, 1 , function(x) paste0(toString(x[5]), toString(x[2])))
phirst1 <- subset(phirst1, select=c(hhid, iid, sample.date))
phirst1 <- ddply(phirst1, .(iid), mutate, age=rep(abs(rnorm(n=sample.size, mean=20, sd=9.5)), each=4, length.out=4))
phirst1$hiv <- if_else(phirst1$age<1,rep(sample(c("0","0"),size=832, replace=TRUE), each=4, length.out=832),
               if_else(phirst1$age>=1,rep(sample(c("0","1"),size=832, replace=TRUE), each=4, length.out=832),NULL))

#4 individuals per household
phirst2 <- data.frame(hid=rep(letters[1:26], each=20), pid=rep.int(1:4, 5))
phirst2$sample.date <- seq(as.Date('2018/01/01'), as.Date('2019/06/04'), by="day")
phirst2$prefix <- "b"
phirst2$hhid <- apply(phirst2, 1 , function(x) paste0(toString(x[4]), toString(x[1])))
phirst2$iid <- apply(phirst2, 1 , function(x) paste0(toString(x[5]), toString(x[2])))
phirst2 <- subset(phirst2, select=c(hhid, iid, sample.date))
phirst2 <- ddply(phirst2, .(iid), mutate, age=rep(abs(rnorm(n=sample.size, mean=20, sd=9.5)), each=5, length.out=5))
phirst2$hiv <-0

#5 individuals per household
phirst3 <- data.frame(hid=rep(letters[1:26], each=30), pid=rep.int(1:5, 6))
phirst3$sample.date <- seq(as.Date('2018/01/01'), as.Date('2020/02/19'), by="day")
phirst3$prefix <- "c"
phirst3$hhid <- apply(phirst3, 1 , function(x) paste0(toString(x[4]), toString(x[1])))
phirst3$iid <- apply(phirst3, 1 , function(x) paste0(toString(x[5]), toString(x[2])))
phirst3 <- subset(phirst3, select=c(hhid, iid, sample.date))
phirst3 <- ddply(phirst3, .(iid), mutate, age=rep(abs(rnorm(n=sample.size, mean=20, sd=9.5)), each=6, length.out=6))
phirst3$hiv <- if_else(phirst3$age<1,rep(sample(c("0","0"),size=780, replace=TRUE), each=6, length.out=780),
               if_else(phirst3$age>=1,rep(sample(c("0","1"),size=780, replace=TRUE), each=6, length.out=780),NULL))

#15 individuals per household
phirst4 <- data.frame(hid=rep(letters[1:26], each=105), pid=rep.int(1:15, 7))
phirst4$sample.date <- seq(as.Date('2018/01/01'), as.Date('2025/06/22'), by="day")
phirst4$prefix <- "d"
phirst4$hhid <- apply(phirst4, 1 , function(x) paste0(toString(x[4]), toString(x[1])))
phirst4$iid <- apply(phirst4, 1 , function(x) paste0(toString(x[5]), toString(x[2])))
phirst4 <- subset(phirst4, select=c(hhid, iid, sample.date))
phirst4 <- ddply(phirst4, .(iid), mutate, age=rep(abs(rnorm(n=sample.size, mean=20, sd=9.5)), each=7, length.out=7))
phirst4$hiv <-0

#combine all datasets
phirst <- rbind(phirst1, phirst2, phirst3, phirst4)
remove(phirst1,phirst2,phirst3,phirst4)

#create sampling points for each individual
phirst <- ddply(phirst, .(iid), mutate, sample.point = seq_along(sample.date))

#define carriage state for each individuals
phirst$state <- as.integer(sample(c("1","2"), size=sample.size, replace=TRUE))

#define missing state with 15% missing swabs
phirst$state[sample(1:length(phirst$state),730)] <- NA
phirst$missing <- if_else(is.na(phirst$state),1,0)

#define serotype carried by each individuals
phirst$serotype <- if_else(phirst$state==2,sample(c("1","3","5","6A","6B","14","19A","19F","23F"),size=sample.size,replace=TRUE),
                   if_else(phirst$state==1, NULL,NULL)) #serotype
                   phirst$serotype[sample(1:length(phirst$serotype),1460)] <- NA


#phirst$Tobs <- c(1L, phirst$uid[-1]!=phirst$uid[-nrow(phirst)]) 


#===============hidden markov models using "msm" R package===============

#---------------MODEL1---------------
#show the state table
statetable.msm(state,iid,data=phirst)
                   
#define initial values of a transition intensity matrix Q
model1.Q <- rbind(c(0.5,0.5), c(0.5,0.5))
rownames(model1.Q) <- c("Clear","Carrying")
colnames(model1.Q) <- c("Clear","Carrying")

#define initial values of an emission matrix E
model1.E <- rbind(c(0.5,0.5), c(0,0.5))
rownames(model1.E) <- c("SwabNeg","SwabPos")
colnames(model1.E) <- c("Clear","Carrying")

#create an outcome distribution of each observed state (swabs)
model1.HMM <- list(hmmBinom(1,0.5), hmmBinom(1,0.5))

#fit a 2-state hidden markov model accounting for false negatives
model1.msm <- msm(state~sample.point, subject=iid, data=phirst,
                  qmatrix=model1.Q,
                  ematrix=model1.E,
                  hmodel=model1.HMM,
                  method="BFGS")

#maximum likelihood estimates of intensity matrix and outcome distributions
printold.msm(model1.msm)

#transition rates and 95%CI of HMM
sojourn.msm(model1.msm)

#most likely true series of states underlying the data
viterbi.msm(model1.msm)

#---------------MODEL2---------------
#define initial values of a transition intensity matrix Q
model2.Q <- rbind(c(0,0.3,0), c(0.05,0,0),c(0,0,0))
rownames(model2.Q) <- c("Clear","Carrying","Missing")
colnames(model2.Q) <- c("Clear","Carrying","Missing")

#define initial values of an emission matrix E
model2.E <- rbind(c(0,0.15,0), c(0,0,0), c(0.1,0.2,0))
rownames(model2.E) <- c("SwabNeg","SwabPos","Missing")
colnames(model2.E) <- c("Clear","Carrying","Missing")

#create an outcome distribution of each observed state (swabs)
model2.HMM <- list(hmmBinom(1,0.5), hmmBinom(1,0.5), hmmBinom(1,0.15))

#fit a 2-state hidden markov model accounting for false negatives
model2.msm <- msm(carrige~sample.point, subject=uid, data=phirst,
                  qmatrix=model2.Q,
                  ematrix=model2.E,
                  hmodel=model2.HMM,
                  method="BFGS")

#maximum likelihood estimates of intensity matrix and outcome distributions
printold.msm(model2.msm)

#transition rates and 95%CI of HMM
sojourn.msm(model2.msm)

#most likely true series of states underlying the data
viterbi.msm(model2.msm)



