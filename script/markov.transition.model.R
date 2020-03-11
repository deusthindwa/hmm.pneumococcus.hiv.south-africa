#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#1/10/2019 - 24/1/2020


#===============load required packages into memory

phirst.packages <- c("tidyverse","dplyr","plyr","msm","timetk","gridExtra","curl","minqa","table1",
                    "lubridate","magrittr","data.table","parallel","foreign","readstata13","ggpubr",
                    "wakefield","zoo","janitor","rethinking","doParallel","scales","msmtools","nlme")
lapply(phirst.packages, library, character.only=TRUE)

#---------------load all phirst datasets (household, master and follow-up)
phirst.hh <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_household.dta",generate.factors=T)
phirst.ms <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_master2.dta",generate.factors=T)
phirst.fu <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_follow-up1.dta",generate.factors=T)
phirst.fu.abx <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_follow-up2.dta",generate.factors=T)

#---------------subset to get required variables in household, and master datasets
phirst.hh <- subset(phirst.hh, is.na(reason_hh_not_inc))
phirst.hh <- subset(phirst.hh, select=c(hh_id,hh_mems_11swabs))

phirst.ms <- subset(phirst.ms, ind_elig_inc=="Yes")
phirst.ms <- subset(phirst.ms, select=c(hh_id,ind_id,year,site,sex,age_at_consent,hiv_status,arv_current_self,arv_current_vl,cd4_count_1,
                                        cd4_count_2,cd4_count_3,cd4_count_4,cd4_count_5,cd4_count_6,cd4_count_7,
                                        cd4_count_8,pcv6wks,pcv14wks,pcv9mth,alcohol,anysmokenow))
phirst.ms <- merge(phirst.ms, phirst.hh, by="hh_id")


#===============prepare master (ms) dataset for baseline demographic characteristics

#---------------studysite
phirst.ms$site <- as.factor(if_else(phirst.ms$site=="Agincourt","Agincourt",
                          if_else(phirst.ms$site=="Klerksdorp","Klerksdorp",NULL)))

#---------------sex category
phirst.ms$sex <- as.factor(if_else(phirst.ms$sex=="Male","Male",
                         if_else(phirst.ms$sex=="Female","Female",NULL)))

#---------------age category
phirst.ms$agecont <- phirst.ms$age_at_consent
phirst.ms$age <- as.factor(if_else(phirst.ms$age_at_consent<5,"Child",
                         if_else(phirst.ms$age_at_consent>=5,"Adult",NULL)))

#---------------hiv status
phirst.ms$hiv <- as.factor(if_else(phirst.ms$hiv_status=="Negative","Negative",
                         if_else(phirst.ms$hiv_status=="Positive","Positive",NULL)))

#---------------ART status
phirst.ms$art1 <- as.factor(if_else(phirst.ms$arv_current_vl=="Yes","Yes",
                         if_else(phirst.ms$arv_current_vl=="No","No",NULL)))

phirst.ms$art2 <- as.factor(if_else(phirst.ms$arv_current_self=="Yes","Yes",
                                    if_else(phirst.ms$arv_current_self=="No","No",NULL)))

#---------------CD4+ cell count
phirst.ms$cd4 <- rowMeans(cbind(phirst.ms$cd4_count_1,phirst.ms$cd4_count_2,phirst.ms$cd4_count_3,phirst.ms$cd4_count_4,
                            phirst.ms$cd4_count_5,phirst.ms$cd4_count_6,phirst.ms$cd4_count_7,phirst.ms$cd4_count_8), na.rm=TRUE)
phirst.ms$cd4cont <- if_else(phirst.ms$hiv=="Positive",as.integer(phirst.ms$cd4),NULL)

phirst.ms$cd4 <- if_else(phirst.ms$cd4<=350 & phirst.ms$age=="Adult","Low",
                         if_else(phirst.ms$cd4>350 & phirst.ms$age=="Adult","High",
                                 if_else(phirst.ms$cd4<=750 & phirst.ms$age=="Child","Low",
                                         if_else(phirst.ms$cd4>750 & phirst.ms$age=="Child","High", NULL))))

phirst.ms$cd4 <- as.factor(if_else(phirst.ms$hiv=="Positive",phirst.ms$cd4,NULL))

#---------------pcv status
phirst.ms$pcv6w <- if_else(phirst.ms$pcv6wks=="Yes","Yes",if_else(phirst.ms$pcv6wks=="No","No",NULL))
phirst.ms$pcv6w <- as.factor(if_else(phirst.ms$age=="Child",phirst.ms$pcv6w,NULL))

phirst.ms$pcv14w <- if_else(phirst.ms$pcv14wks=="Yes","Yes",if_else(phirst.ms$pcv14wks=="No","No",NULL))
phirst.ms$pcv14w <- as.factor(if_else(phirst.ms$age=="Child",phirst.ms$pcv14w,NULL))

phirst.ms$pcv9m <- if_else(phirst.ms$pcv9mth=="Yes","Yes",if_else(phirst.ms$pcv9mth=="No","No",NULL))
phirst.ms$pcv9m <- as.factor(if_else(phirst.ms$age=="Child",phirst.ms$pcv9m,NULL))

#---------------alcohol consumption
phirst.ms$alcohol <- if_else(phirst.ms$alcohol==1,"Yes",if_else(phirst.ms$alcohol==0,"No",NULL))
phirst.ms$alcohol <- as.factor(if_else(phirst.ms$age=="Adult",phirst.ms$alcohol,NULL))

#---------------smoking status
phirst.ms$smoke <- if_else(phirst.ms$anysmokenow==1,"Yes",if_else(phirst.ms$anysmokenow==0,"No",NULL))
phirst.ms$smoke <- as.factor(if_else(phirst.ms$age=="Adult",phirst.ms$smoke,NULL))

#---------------household size
phirst.ms$hhsize <- as.integer(phirst.ms$hh_mems_11swabs)

#---------------number of HIV+ adults in the household
phirst.hhhiv <- subset(phirst.ms,select=c(hh_id,age,hiv))
phirst.hhhiv$hhiv <- if_else(phirst.hhhiv$age=="Adult" & phirst.hhhiv$hiv=="Positive",1L,
                             if_else(phirst.hhhiv$age=="Adult" & phirst.hhhiv$hiv=="Negative",0L,
                                     if_else(phirst.hhhiv$age=="Child" & phirst.hhhiv$hiv=="Negative",0L,
                                             if_else(phirst.hhhiv$age=="Child" & phirst.hhhiv$hiv=="Positive",0L,NULL))))
phirst.hhhiv$age <- phirst.hhhiv$hiv <- NULL
phirst.hhhiv <- setDT(phirst.hhhiv)[,list(ahiv=sum(hhiv)), by=.(hh_id)]
phirst.hhhiv$ahivc <- as.factor(if_else(phirst.hhhiv$ahiv>=1,"Yes",if_else(phirst.hhhiv$ahiv==0,"No",NULL)))

#---------------final master dataset
phirst.ms <- subset(phirst.ms, select=c(hh_id,ind_id,year,hhsize,age,agecont,site,sex,hiv,art1,art2,cd4,cd4cont,pcv6w,pcv14w,pcv9m,alcohol,smoke))
phirst.ms <- merge(x=phirst.ms, y=phirst.hhhiv, by="hh_id", all.x=TRUE)
remove(phirst.hhhiv)

#---------------baseline demographic characteristics (Table 1)
phirst.render <- function(x,name,...){
  if(!is.numeric(x)) return(render.categorical.default(na.omit(x)))
  }
table1(~site+sex+hiv+art1+art2+cd4+ahivc+pcv6w+pcv14w+pcv9m+smoke+alcohol|age, data=phirst.ms,
       topclass="Rtable1-grid Rtable1-shade Rtable1-times",render=phirst.render)

#age count
mean(phirst.ms[phirst.ms$age=="Child","agecont"], na.rm=TRUE)
sd(phirst.ms[phirst.ms$age=="Child","agecont"], na.rm=TRUE)
mean(phirst.ms[phirst.ms$age=="Adult","agecont"], na.rm=TRUE)
sd(phirst.ms[phirst.ms$age=="Adult","agecont"], na.rm=TRUE)
mean(phirst.ms$agecont, na.rm=TRUE)
sd(phirst.ms$agecont, na.rm=TRUE)   

#cd4 count
mean(phirst.ms[phirst.ms$age=="Child","cd4cont"], na.rm=TRUE)
sd(phirst.ms[phirst.ms$age=="Child","cd4cont"], na.rm=TRUE)
mean(phirst.ms[phirst.ms$age=="Adult","cd4cont"], na.rm=TRUE)
sd(phirst.ms[phirst.ms$age=="Adult","cd4cont"], na.rm=TRUE)
mean(phirst.ms$cd4cont, na.rm=TRUE)
sd(phirst.ms$cd4cont, na.rm=TRUE)         


#===============prepare follow-up dataset for modelling

#---------------create follow-up continuous time in days
phirst.fu.abx <- subset(phirst.fu.abx, select=c(finalid,antibiotic))
phirst.fu.abx <- rename(phirst.fu.abx, c("finalid"="visit_id","antibiotic"="abx"))
phirst.fu <- subset(phirst.fu, select=c(ind_id,visit_id,visit_date,visit,npspne,npspneload))
phirst.fu <- arrange(phirst.fu,visit_id)
phirst.fu$visit_date <- ymd(phirst.fu$visit_date)
phirst.fu$startdate <- if_else(phirst.fu$visit==1L,phirst.fu$visit_date,NULL)
phirst.fu$startdate <- na.locf(phirst.fu$startdate)
phirst.fu$dys <- difftime(phirst.fu$visit_date,phirst.fu$startdate,units="days")

#---------------merge follow-up dataset to master dataset
phirst.ms <- arrange(phirst.ms, ind_id)
phirst.fu <- merge(phirst.fu, phirst.ms, by="ind_id")
phirst.fu <- merge(x=phirst.fu, y=phirst.fu.abx, by="visit_id", all.x=TRUE)

#---------------number of PNC+ adults in the household per visit
phirst.hhpnc <- subset(phirst.fu,select=c(visit_id,age,visit,npspne))
phirst.hhpnc$hh_id <- substr(phirst.hhpnc$visit_id,0,4)
phirst.hhpnc$hpnc <- if_else(phirst.hhpnc$age=="Adult" & phirst.hhpnc$npspne==1,1L,
                             if_else(phirst.hhpnc$age=="Adult" & phirst.hhpnc$npspne==0,0L,
                                     if_else(phirst.hhpnc$age=="Child" & phirst.hhpnc$npspne==0,0L,
                                             if_else(phirst.hhpnc$age=="Child" & phirst.hhpnc$npspne==1,0L,NULL))))
phirst.hhpnc$age <- phirst.hhpnc$npspne <- NULL
phirst.hhpnc <- arrange(phirst.hhpnc,hh_id,visit_id)
phirst.hhpnc <- setDT(phirst.hhpnc)[,list(apnc=sum(hpnc)),by=.(hh_id,visit)]
phirst.fu <- merge(x=phirst.fu, y=phirst.hhpnc, by=c("hh_id","visit"), all.x=TRUE)
remove(phirst.hhpnc)

#---------------subset only required variables in follow-up dataset
phirst.fu <- subset(phirst.fu, select=c(visit_id,ind_id,dys,npspne,npspneload,age,sex,hhsize,hiv,art1,ahiv,ahivc,apnc,abx))

#---------------integerise variable categories for modelling
phirst.fu$dys <- as.integer(phirst.fu$dys)

phirst.fu <- rename(phirst.fu, c("npspne" = "state"))
phirst.fu$state <- recode_factor(phirst.fu$state,`1`="Carry",`0`="Clear")
phirst.fu$state <- recode_factor(phirst.fu$state,`Carry`=2,`Clear`=1)
phirst.fu$state <- if_else(phirst.fu$state==2,2L,if_else(phirst.fu$state==1,1L,NULL))
phirst.fu$statem <- phirst.fu$state
phirst.fu$statem[is.na(phirst.fu$statem)] <- 999L
phirst.fu$obst <- if_else(phirst.fu$statem==999L,1L,0L)

phirst.fu$age <- recode_factor(phirst.fu$age,`Adult`=1L,`Child`=0L)
phirst.fu$age <- if_else(phirst.fu$age==1,1L,if_else(phirst.fu$age==0,0L,NULL))

phirst.fu$sex <- recode_factor(phirst.fu$sex,`Male`=1L,`Female`=0L)
phirst.fu$sex <- if_else(phirst.fu$sex==1,1L,if_else(phirst.fu$sex==0,0L,NULL))

phirst.fu$hiv <- recode_factor(phirst.fu$hiv,`Positive`=1L,`Negative`=0L)
phirst.fu$hiv <- if_else(phirst.fu$hiv==1,1L,if_else(phirst.fu$hiv==0,0L,NULL))

phirst.fu$art1 <- recode_factor(phirst.fu$art1,`Yes`=1L,`No`=0L)
phirst.fu$art1 <- if_else(phirst.fu$art1==1,1L,if_else(phirst.fu$art1==0,0L,NULL))
phirst.fu$art1[is.na(phirst.fu$art1) & phirst.fu$hiv==0] <- 0L

phirst.fu$ahiv <- as.integer(phirst.fu$ahiv)
phirst.fu$ahivc <- if_else(phirst.fu$ahiv>=1,1L,if_else(phirst.fu$ahiv==0,0L,NULL))

phirst.fu$apnc <- as.integer(phirst.fu$apnc)
phirst.fu$apncc <- if_else(phirst.fu$apnc>=1,1L,if_else(phirst.fu$apnc==0,0L,NULL))

phirst.fu$abx <- as.integer(phirst.fu$abx)
phirst.fu$abx[is.na(phirst.fu$abx)] <- 0L
phirst.fu$hh_id <- substr(phirst.fu$visit_id,1,4)

phirst.trans <- subset(phirst.fu, select=c(hh_id,dys,state))
phirst.trans <- subset(phirst.trans, dys==0)
phirst.trans <- subset(phirst.trans, select=c(hh_id,state))
phirst.trans <- subset(phirst.trans, !is.na(phirst.trans$state))
phirst.trans$seqx <- seq(from=1, to=nrow(phirst.trans))
phirst.trans <- dcast(phirst.trans, hh_id ~ state, value.var="seqx")
phirst.trans <- rename(phirst.trans, c("1"="clear", "2"="carry"))
phirst.trans <- subset(phirst.trans,carry==0)
phirst.trans$trans <- 1L
phirst.trans$clear <- phirst.trans$carry <- NULL
phirst.fu <- merge(x=phirst.fu, y=phirst.trans, all.x=TRUE)
phirst.fu$trans[is.na(phirst.fu$trans)] <- 2L

#---------------final variables in follow-up dataset
phirst.fu <- subset(phirst.fu, select=c(visit_id,ind_id,dys,state,statem,obst,npspneload,hhsize,age,sex,hiv,art1,ahiv,ahivc,apnc,apncc,abx,trans))
phirst.fu <- subset(phirst.fu, !is.na(hiv))

#---------------follow-up characteristics of participants (figure 2)
source('~/Rproject/Markov.Model/script/fig2.R')


#===============hidden Markov modelling of carriage dynamics within houshold and from community

#---------------show transition frequency
phirst.fu <- arrange(phirst.fu,visit_id)
statetable.msm(state,ind_id,data=phirst.fu)

#---------------initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.9),
                  c(0.9,0.0))
rownames(matrix.Q) <- c("Clear","Carry")
colnames(matrix.Q) <- c("Clear","Carry")

#---------------initiate emission matrix E
matrix.E <- rbind(c(1.0,0.0), 
                  c(0.1,0.9))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carry")

#---------------fitting hidden Markov models with misclassifications
phirst.fu <- arrange(phirst.fu,visit_id)

p.model1 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans,
                censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model2 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans+apncc,
                censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model3 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans+art1,
                censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model4 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans+abx,
                censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model5 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans+apncc+art1,
                censor=999,censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model6 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans+apncc+abx,
                censor=999,censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model7 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+trans+art1+abx,
                censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model8 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=~age+hiv+ahivc+apncc+art1+abx,
                censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=200000))

#---------------comparing the AIC and BIC of the fitted models 
AIC(p.model1,p.model2,p.model3,p.model4,p.model5,p.model6,p.model7)
AIC(p.model1,k=log(length(phirst.fu)))
AIC(p.model2,k=log(length(phirst.fu)))
AIC(p.model3,k=log(length(phirst.fu)))
AIC(p.model4,k=log(length(phirst.fu)))
AIC(p.model5,k=log(length(phirst.fu)))
AIC(p.model6,k=log(length(phirst.fu)))
AIC(p.model7,k=log(length(phirst.fu)))

#print out the baseline intensities, and emission probability of the selected model
printnew.msm(p.model5)

#---------------convergence of the selected model
phirst.fu <- arrange(phirst.fu,visit_id)
j=0.05
k=2.00
for(i in 1:5){
matrix.Qc <- rbind(c(0.0,j), c(k,0.0))
rownames(matrix.Qc) <- c("Clear","Carry")
colnames(matrix.Qc) <- c("Clear","Carry")

matrix.Ec <- rbind(c(1.0,0.0), c(0.15,0.85))
colnames(matrix.Ec) <- c("SwabNeg","SwabPos")
rownames(matrix.Ec) <- c("Clear","Carry")

sink("~/Rproject/Markov.Model.Resources/data/convergence.txt", append=TRUE)
m.converge <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                   qmatrix=matrix.Qc, ematrix=matrix.Ec,
                   covariates=~age+hiv+ahivc+apncc+art1,
                   censor=999, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                   control=list(maxit=100000,trace=1,REPORT=1))
sink()
j=j+0.1
k=k-0.4
}

#---------------onbserved versus predicted values of selected model
p.ObsExp <-msm(statem~dys, subject=ind_id, data=phirst.fu,
              qmatrix=matrix.Q,
              covariates=~age+hiv+ahivc+trans+apncc+art1,
              censor=999, 
              censor.states=c(1,2), 
              est.initprobs=T,
              opt.method="bobyqa", control=list(maxfun=100000))

m.OEDS <- tk_tbl(prevalence.msm(p.ObsExp,times=c(14,28,42,56,70,84,98,112,126,140,154,168,182,196,210,224,238,252,266,280)),preserve_index=TRUE,rename_index="Time")
m.OEDS <- subset(m.OEDS, select=c(Time, Observed.State.1, Observed.State.2, Observed.percentages.State.1, Observed.percentages.State.2,
                                  Expected.Clear,Expected.Carry,Expected.percentages.Clear,Expected.percentages.Carry))
m.OEDS <- rename(m.OEDS, c("Observed.State.1"="obs.clear", "Observed.State.2"="obs.carry","Observed.percentages.State.1"="obs.p.clear","Observed.percentages.State.2"="obs.p.carry",
                           "Expected.Clear"="exp.clear","Expected.Carry"="exp.carry","Expected.percentages.Clear"="exp.p.clear","Expected.percentages.Carry"="exp.p.carry"))

m.OEDS$lci.clear=m.OEDS$exp.p.clear/100-(1.96*sqrt(m.OEDS$exp.p.clear/100*(1-m.OEDS$exp.p.clear/100)/m.OEDS$exp.clear)) 
m.OEDS$uci.clear=m.OEDS$exp.p.clear/100+(1.96*sqrt(m.OEDS$exp.p.clear/100*(1-m.OEDS$exp.p.clear/100)/m.OEDS$exp.clear))
m.OEDS$lci.carry=m.OEDS$exp.p.carry/100-(1.96*sqrt(m.OEDS$exp.p.carry/100*(1-m.OEDS$exp.p.carry/100)/m.OEDS$exp.carry)) 
m.OEDS$uci.carry=m.OEDS$exp.p.carry/100+(1.96*sqrt(m.OEDS$exp.p.carry/100*(1-m.OEDS$exp.p.carry/100)/m.OEDS$exp.carry))

#---------------plot of model parameter convergence, and observed and predictions (figure 3)
dev.off()
source('~/Rproject/Markov.Model/script/fig3.R')

#---------------viterbi algorithm
phirst.vi <- viterbi.msm(p.model5)
phirst.vi$time <- as.integer(phirst.vi$time)
phirst.vi$observed <- as.integer(phirst.vi$observed)
phirst.vi$fitted <- as.integer(phirst.vi$fitted)
pstate <- as.data.frame(phirst.vi$pstate)
phirst.vi$probhs1 <- pstate$V1
phirst.vi$probhs2 <- pstate$V2
remove(pstate); phirst.vi$pstate <- NULL
phirst.vi$observed2 <- if_else(phirst.vi$observed==1L,"Clear","Carry")
phirst.vi$fitted2 <- if_else(phirst.vi$fitted==1L,"Clear","Carry")

#---------------plot of observed states, viterbi states and viterbi probability (figure 4)
source('~/Rproject/Markov.Model/script/fig4.R')

#---------------plot of within household acquisition and clearance rates and probabilities from selected model (figure 5)
p.model5a <- msm(statem~dys, subject=ind_id, data=subset(phirst.fu,trans==2),
                 qmatrix=matrix.Q, ematrix=matrix.E,
                 covariates=~age+hiv+ahivc+apncc+art1,
                 censor=999,censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                 opt.method="bobyqa", control=list(maxfun=100000))
source('~/Rproject/Markov.Model/script/fig5.R')

#---------------plot of community acquisition rates and probabilities from selected model (figure 6)
p.model5b <- msm(statem~dys, subject=ind_id, data=subset(phirst.fu, trans==1),
                 qmatrix=matrix.Q, ematrix=matrix.E,
                 covariates=~age+hiv+ahivc+apncc+art1,
                 censor=999,censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                 opt.method="bobyqa", control=list(maxfun=100000))
source('~/Rproject/Markov.Model/script/fig6.R')

#---------------plot of other important epidemiological characterisations of carriage (figure 7)
source('~/Rproject/Markov.Model/script/fig7.R')
