#Written by Deus Thindwa
#Does household HIV-exposure increase infant risk of pneumococcal carriage acquisition?
#An analysis using markov transition model
#PhD chapter 2
#14/05/2019

#load required packages into memory
phirst.packages <-c("tidyverse","plyr")
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

#descriptive analysis






