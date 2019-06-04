#Written by Deus Thindwa
#Does household HIV-exposure increase infant risk of pneumococcal carriage acquisition?
#An analysis using markov transition model
#PhD chapter 2
#14/05/2019

#load required packages into memory
phirst.packages <-c("tidyverse","plyr", "msm")
lapply(phirst.packages, library, character.only=TRUE)

#simulate a phirst dataset
set.seed(12345)
sample.size=520
phirst <- data.frame(hhid=rep(letters[1:26], each=20), pid=rep.int(1:4, 5)) #household id, individual id
phirst$uid <- apply(phirst, 1 , function(x) paste0(toString(x[1]), toString(x[2]))) #household individual id
phirst$sample.date <- seq(as.Date('2018/01/01'), as.Date('2019/06/04'), by="day") #sample date
phirst <- ddply(phirst, .(uid), mutate, sample.point = seq_along(sample.date)) #sample point
phirst <- ddply(phirst, .(uid), mutate, age=rep(rnorm(n=sample.size, mean=20, sd=9.5), each=5, length.out=5)) #age
phirst$carrige <- sample(c("negative","positive"), size=sample.size, replace=TRUE) #carriage
phirst$serotype <- if_else(phirst$carrige=="positive", sample(c("1","3","5","6A","6B","14","19A","19F","23F"), size=sample.size, replace=TRUE),
                   if_else(phirst$carrige=="negative", NULL,NULL)) #serotype
phirst$hivstatus <- if_else(phirst$age<1, rep(sample(c("unexposed","exposed"), size=sample.size, replace=TRUE), each=5, length.out=520),
                           if_else(phirst$age>=1, rep(sample(c("negative","positive"), size=sample.size, replace=TRUE), each=5, length.out=520),NULL)) #hiv status

#analysis using "msm" R package
#load coronary allograft vasculopathy dataset
cav <- msm::cav

#exclude pdiag missing obs
cav <- cav[!is.na(cav$pdiag),]

#summarise number of transitions between states
statetable.msm(cav$state, cav$PTNUM)

#construct a transition matrix with initial guess values
transitionM <- rbind(c(0, 0.25, 0, 0.25), c(0.166, 0, 0.166, 0.166), c(0, 0.25, 0, 0.25), c(0, 0, 0, 0))
rownames(transitionM) <- colnames(transitionM) <- c("Well", "Mild", "Severe", "Death")

#maximum likelihood of the transition matrix
transitionMLE <- msm(state ~ years, subject=PTNUM, data=cav, qmatrix=transitionM, death=4)

#display the fitted transition probabilities at specified time interval
pmatrix.msm(transitionMLE, t=1, ci="normal")

#fit transition probability matrix with covariates
cav$ihd <- as.numeric(cav[,"pdiag"]=="IHD")
transitionMLE.cov <- msm(state ~ years, subject=PTNUM, data=cav, covariates=~dage+ihd, qmatrix=transitionM,death=4,
                         method="BFGS", control=list(fnscale=4000, maxit=10000))

#hazard ratios (cov value comparison) of transitioning between states
hazard.msm(transitionMLE.cov)

#calculate transition matrix of specified cov value
qmatrix.msm(transitionMLE.cov, covariates=list(dage=50,ihd=1))

#model comparisons
lrtest.msm(transitionMLE,transitionMLE.cov)










