#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#1/10/2019 - 24/1/2020

#===============load required packages into memory
phirst.packages <-c("tidyverse","plyr","msm","timetk","gridExtra","curl","dplyr","minqa","lubridate","magrittr","data.table","parallel","foreign","readstata13","wakefield","zoo","janitor","rethinking","doParallel","scales","msmtools")
lapply(phirst.packages, library, character.only=TRUE)

#---------------load all phirst datasets (household, master and follow-up)
phirst.hh <- read.dta13("~/Rproject/Markov.Model.Resources/phirst_household.dta",generate.factors=T)
phirst.ms <- read.dta13("~/Rproject/Markov.Model.Resources/phirst_master.dta",generate.factors=T)
phirst.fu <- read.dta13("~/Rproject/Markov.Model.Resources/phirst_follow-up.dta",generate.factors=T)

#---------------subset to get required variables in household, and master datasets
phirst.hh <- subset(phirst.hh, is.na(reason_hh_not_inc))
phirst.hh <- subset(phirst.hh, select=c(hh_id,hh_mems_11swabs))
phirst.ms <- subset(phirst.ms, ind_elig_inc=="Yes")
phirst.ms <- subset(phirst.ms, select=c(hh_id,ind_id,year,site,sex,age_at_consent,hiv_status,arv_current,cd4_count_1,cd4_count_2,cd4_count_3,cd4_count_4,cd4_count_5,cd4_count_6,cd4_count_7,cd4_count_8,vir_ld_rna_cop_3,vir_ld_rna_cop_4,vir_ld_rna_cop_5,vir_ld_rna_cop_6,vir_ld_rna_cop_7,vir_ld_rna_cop_8,pcv6wks,pcv14wks,pcv9mth,alcohol,anysmokenow))
phirst.ms <- merge(phirst.ms,phirst.hh)
rm(phirst.hh)

#===============prepare master (ms) dataset for baseline demographic characteristics

#---------------studysite
phirst.ms$site <- if_else(phirst.ms$site=="Agincourt","Agincourt",
                          if_else(phirst.ms$site=="Klerksdorp","Klerksdorp",NULL))

#---------------sex category
phirst.ms$sex <- if_else(phirst.ms$sex=="Male","Male",
                         if_else(phirst.ms$sex=="Female","Female",NULL))

#---------------age category
phirst.ms$age <- if_else(phirst.ms$age_at_consent<5,"Child",
                         if_else(phirst.ms$age_at_consent>=5,"Adult",NULL))

#---------------hiv status
phirst.ms$hiv <- if_else(phirst.ms$hiv_status=="Negative","Negative",
                         if_else(phirst.ms$hiv_status=="Positive","Positive",NULL))

#---------------ART status
phirst.ms$art <- if_else(phirst.ms$arv_current=="Yes","Yes",
                         if_else(phirst.ms$arv_current=="No","No","No"))

#---------------CD4+ cell count
phirst.ms$cd4 <- if_else(!is.na(phirst.ms$cd4_count_8),phirst.ms$cd4_count_8,
                         if_else(!is.na(phirst.ms$cd4_count_7),phirst.ms$cd4_count_7,
                                 if_else(!is.na(phirst.ms$cd4_count_6),phirst.ms$cd4_count_6,
                                         if_else(!is.na(phirst.ms$cd4_count_5),phirst.ms$cd4_count_5,
                                                 if_else(!is.na(phirst.ms$cd4_count_4),phirst.ms$cd4_count_4,
                                                         if_else(!is.na(phirst.ms$cd4_count_3),phirst.ms$cd4_count_3,
                                                                 if_else(!is.na(phirst.ms$cd4_count_2),phirst.ms$cd4_count_2,
                                                                         if_else(!is.na(phirst.ms$cd4_count_1),phirst.ms$cd4_count_1, NULL))))))))
phirst.ms$cd4 <- if_else(phirst.ms$cd4<=350 & phirst.ms$age=="Adult","Low",
                         if_else(phirst.ms$cd4>350 & phirst.ms$age=="Adult","High",
                                 if_else(phirst.ms$cd4<=750 & phirst.ms$age=="Child","Low",
                                         if_else(phirst.ms$cd4>750 & phirst.ms$age=="Child","High", NULL))))
phirst.ms$cd4 <- if_else(phirst.ms$hiv=="Positive",phirst.ms$cd4,NULL)

#---------------viral load
phirst.ms$vl <- if_else(!is.na(phirst.ms$vir_ld_rna_cop_8),phirst.ms$vir_ld_rna_cop_8,
                        if_else(!is.na(phirst.ms$vir_ld_rna_cop_7),phirst.ms$vir_ld_rna_cop_7,
                                if_else(!is.na(phirst.ms$vir_ld_rna_cop_6),phirst.ms$vir_ld_rna_cop_6,
                                        if_else(!is.na(as.integer(phirst.ms$vir_ld_rna_cop_5)), as.integer(phirst.ms$vir_ld_rna_cop_5),
                                                if_else(!is.na(phirst.ms$vir_ld_rna_cop_4),phirst.ms$vir_ld_rna_cop_4,
                                                        if_else(!is.na(phirst.ms$vir_ld_rna_cop_3),phirst.ms$vir_ld_rna_cop_3, NULL))))))
phirst.ms$vl <- if_else(phirst.ms$vl<=1000,"Low",
                        if_else(phirst.ms$vl>1000,"High", NULL))

#---------------PCV status
phirst.ms$pcv6w <- if_else(phirst.ms$pcv6wks=="Yes","Yes",
                           if_else(phirst.ms$pcv6wks=="No","No",NULL))
phirst.ms$pcv14w <- if_else(phirst.ms$pcv14wks=="Yes","Yes",
                            if_else(phirst.ms$pcv14wks=="No","No",NULL))
phirst.ms$pcv9m <- if_else(phirst.ms$pcv9mth=="Yes","Yes",
                           if_else(phirst.ms$pcv9mth=="No","No",NULL))

#---------------alcohol consumption
phirst.ms$alcohol <- if_else(phirst.ms$alcohol==1,"Yes",
                             if_else(phirst.ms$alcohol==0,"No",NULL))

#---------------smoking status
phirst.ms$smoke <- if_else(phirst.ms$anysmokenow==1,"Yes",
                           if_else(phirst.ms$anysmokenow==0,"No",NULL))

#---------------household size
phirst.ms$hhsize <- phirst.ms$hh_mems_11swabs

#---------------number of HIV+ adults in the household
phirst.hhhiv <- subset(phirst.ms,select=c(hh_id,age,hiv))
phirst.hhhiv$hhiv <- if_else(phirst.hhhiv$age=="Adult" & phirst.hhhiv$hiv=="Positive",1L,
                             if_else(phirst.hhhiv$age=="Adult" & phirst.hhhiv$hiv=="Negative",0L,
                                     if_else(phirst.hhhiv$age=="Child" & phirst.hhhiv$hiv=="Negative",0L,
                                             if_else(phirst.hhhiv$age=="Child" & phirst.hhhiv$hiv=="Positive",0L,NULL))))
phirst.hhhiv$age <- phirst.hhhiv$hiv <- NULL
phirst.hhhiv <- setDT(phirst.hhhiv)[,list(ahiv=sum(hhiv)), by=.(hh_id)]

#---------------final master dataset
phirst.ms <- subset(phirst.ms, select=c(hh_id,ind_id,year,hhsize,age,site,sex,hiv,art,cd4,vl,pcv6w,pcv14w,pcv9m,alcohol,smoke))
phirst.ms <- merge(phirst.ms,phirst.hhhiv)
remove(phirst.hhhiv)

#---------------baseline demographic characteristics (Table 1)
phirst.ms %>% tabyl(age, show_na=FALSE) %>% adorn_pct_formatting(digits=1)

for(i in colnames(phirst.ms[c(5:17)])){
  print(tabyl(phirst.ms[[i]],show_na=FALSE) %>% adorn_pct_formatting(digits=1))
};remove(i)

for(i in colnames(phirst.ms[c(6:17)])){
  print(tabyl(phirst.ms[[i]][phirst.ms$age=="Child"],show_na=FALSE) %>% adorn_pct_formatting(digits=1))
};remove(i)

for(i in colnames(phirst.ms[c(6:17)])){
  print(tabyl(phirst.ms[[i]][phirst.ms$age=="Adult"],show_na=FALSE) %>% adorn_pct_formatting(digits=1))
};remove(i)

#===============prepare follow-up (fu) dataset for modelling

#---------------merge follow-up to master dataset
phirst.fu <- subset(phirst.fu, select=c(ind_id,visit_id,visit_date,visit,npspne,npspneload))
phirst.fu$visit_date <- ymd(phirst.fu$visit_date)
phirst.fu$startdate <- if_else(phirst.fu$visit==1L,phirst.fu$visit_date,NULL)
phirst.fu$startdate <- na.locf(phirst.fu$startdate)
phirst.fu$dys <- difftime(phirst.fu$visit_date,phirst.fu$startdate,units="days")
phirst.fu$pndensity <- if_else(phirst.fu$npspneload>=1000,"High",
                               if_else(phirst.fu$npspneload<1000,"Low",NULL))
phirst.fu <- merge(phirst.fu,phirst.ms)

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
phirst.fu <- merge(phirst.fu,phirst.hhpnc,by=c("hh_id","visit"))
remove(phirst.hhpnc)

#---------------gather only required variables in follow-up dataset
phirst.fu <- subset(phirst.fu, select=c(hh_id,ind_id,visit_id,hhsize,dys,npspne,pndensity,age,sex,hiv,art,cd4,vl,ahiv,apnc))

#---------------integerise variable categories for modelling
phirst.fu$dys <- as.integer(phirst.fu$dys)

phirst.fu <- rename(phirst.fu, c("npspne" = "state"))
phirst.fu$state <- recode_factor(phirst.fu$state,`1`="Carry",`0`="Clear")
phirst.fu$state <- recode_factor(phirst.fu$state,`Carry`=2,`Clear`=1)
phirst.fu$state <- if_else(phirst.fu$state==2,2L,if_else(phirst.fu$state==1,1L,NULL))
phirst.fu$statem <- phirst.fu$state
phirst.fu$statem[is.na(phirst.fu$statem)] <- 999L
phirst.fu$obst <- if_else(phirst.fu$statem==999L,1L,0L)

phirst.fu$pndensity <- recode_factor(phirst.fu$pndensity,`High`=1L,`Low`=0L)
phirst.fu$pndensity <- if_else(phirst.fu$pndensity==1,1L,if_else(phirst.fu$state==0,0L,NULL))

phirst.fu$age <- recode_factor(phirst.fu$age,`Adult`=1L,`Child`=0L)
phirst.fu$age <- if_else(phirst.fu$age==1,1L,if_else(phirst.fu$age==0,0L,NULL))

phirst.fu$sex <- recode_factor(phirst.fu$sex,`Male`=1L,`Female`=0L)
phirst.fu$sex <- if_else(phirst.fu$sex==1,1L,if_else(phirst.fu$sex==0,0L,NULL))

phirst.fu$hiv <- recode_factor(phirst.fu$hiv,`Positive`=1L,`Negative`=0L)
phirst.fu$hiv <- if_else(phirst.fu$hiv==1,1L,if_else(phirst.fu$hiv==0,0L,NULL))

phirst.fu$art <- recode_factor(phirst.fu$art,`Yes`=1L,`No`=0L)
phirst.fu$art <- if_else(phirst.fu$art==1,1L,if_else(phirst.fu$art==0,0L,NULL))

phirst.fu$cd4 <- recode_factor(phirst.fu$cd4,`High`=1L,`Low`=0L)
phirst.fu$cd4 <- if_else(phirst.fu$cd4==1,1L,if_else(phirst.fu$cd4==0,0L,NULL))

phirst.fu$vl <- recode_factor(phirst.fu$vl,`High`=1L,`Low`=0L)
phirst.fu$vl <- if_else(phirst.fu$vl==1,1L,if_else(phirst.fu$vl==0,0L,NULL))

phirst.fu$ahiv <- as.integer(phirst.fu$ahiv)
phirst.fu$ahivc <- if_else(phirst.fu$ahiv>=1,1L,if_else(phirst.fu$ahiv==0,0L,NULL))

phirst.fu$apnc <- as.integer(phirst.fu$apnc)
phirst.fu$apncc <- if_else(phirst.fu$apnc>=1,1L,if_else(phirst.fu$apnc==0,0L,NULL))

#---------------final variables in follow-up dataset
phirst.fu <- subset(phirst.fu, select=c(hh_id,ind_id,visit_id,dys,state,statem,obst,pndensity,hhsize,age,sex,hiv,art,cd4,vl,ahiv,ahivc,apnc,apncc))

#===============markov modelling without transmission assumptions within houshold

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

#---------------fitting time-homogeneous nested models with misclassifications
phirst.fu <- arrange(phirst.fu,visit_id)

p.model1 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model2 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model3 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc+apncc,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model4 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+apncc+art,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model5 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc+art,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model6 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc+apncc+art,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model7 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+apncc+abx,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model8 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc+abx,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model9 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc+apncc+abx,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model10 <- msm(statem~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+ahivc+apncc+art+abx,
                censor=999,censor.states=c(1,2),
                obstrue=obst,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

#---------------comparing the Akaike Information Criterion of the fitted nested models 
AIC(p.model1,p.model2,p.model3,p.model4,p.model5,p.model6)

#print out the baseline intensities, and emission probability of the model with smallest AIC
printnew.msm(p.model3)

#---------------model convergence
q.list <- boot.msm(p.model3,B=3,cores=3,stat=function(p.model3){qmatrix.msm(p.model3)$minus2loglik})

#---------------plot carriage and clearence prevalence of the model with smallest AIC
m.prev.data <- as.data.frame(prevalence.msm(times=c(14,28,42,56,70,84,98,112,126,140,154,168,182,196,210,224,238,252,266,280), p.model3, ci="normal"))
m.prev.perc <- subset(X, select=c(Observed.State.1,Observed.State.2,Expected.estimates.Clear))

dev.off()
par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
plot.prevalence.msm(p.model3, mintime=0,maxtime=288,legend.pos=c(0,100),lwd.obs=2,lty.obs=1,col.obs="blue",lwd.exp=2,lty.exp=1,ci="normal",col.exp="red",cex=0.7,xlab="Time (days)",ylab="% Prevalence")
legend(0,100, legend=c("Observed carriage", "Fitted model"),col=c("blue","red"), lty=1:1, cex=1, lwd=3)

#---------------fitting a model with misclassification
phirst.fu$dyz <- phirst.fu$dys/365.26
p.model2 <- msm(state~dyz, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates=~age+hiv+nahiv1,
                est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=60000))

#---------------pearson-type goodness of fit
options(digits=2)
pearson.msm(p.model2, timegroups=2)

#---------------expected vs observed carriage
dev.off()
par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
plot.prevalence.msm(p.model2, mintime=0,maxtime=288,legend.pos=c(0,100),lwd.obs=2.5,lwd.exp=2.5,cex=0.7,xlab="Time (days)",ylab="% Prevalence")
legend(0, 100, legend=c("Observed carriage", "Fitted model"),col=c("blue","red"), lty=1:3, cex=1, lwd=3)

#---------------Likelihood ratio test
lrtest.msm(p.model1,p.model2)
AIC(p.model1,p.model2)
BIC(p.model1,p.model2)

#---------------transition intensity matrix for a misclassification model with covariates
#acquisition and clearance rates in HIV+ and HIV- children and adults from all HH
set.seed(1988)
qmatrix.msm(p.model2, covariates=list(hiv=1,age=0), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=0,age=0), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=1,age=1), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=0,age=1), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance rates in HIV+ or HIV- children and adults from HH without HIV+ adult(s)
set.seed(1988)
qmatrix.msm(p.model2, covariates=list(hiv=1,age=0,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=0,age=0,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=1,age=1,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=0,age=1,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance rates in HIV+ or HIV- children and adults from HH with HIV+ adult(s)
set.seed(1988)
qmatrix.msm(p.model2, covariates=list(hiv=1,age=0,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=0,age=0,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=1,age=1,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)
qmatrix.msm(p.model2, covariates=list(hiv=0,age=1,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance rates in all children from HH without HIV+ adults
set.seed(1988)
qmatrix.msm(p.model2, covariates=list(age=0,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance rates in all children from HH with HIV+ adults
set.seed(1988)
qmatrix.msm(p.model2, covariates=list(age=0,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)


#---------------transition probability matrix for model with misclassification with covariates (Table 3)
#acquisition and clearance probability in HIV+ and HIV- children and adults from all HH
set.seed(1988)
pmatrix.msm(p.model2, covariates=list(hiv=1,age=0), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0,age=0), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=1,age=1), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0,age=1), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance probability in HIV+ or HIV- children and adults from HH without HIV+ adult(s)
set.seed(1988)
pmatrix.msm(p.model2, covariates=list(hiv=1,age=0,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0,age=0,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=1,age=1,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0,age=1,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance probability in HIV+ or HIV- children and adults from HH with HIV+ adult(s)
set.seed(1988)
pmatrix.msm(p.model2, covariates=list(hiv=1,age=0,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0,age=0,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=1,age=1,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0,age=1,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance probability in all children from HH without HIV+ adult(s)
set.seed(1988)
pmatrix.msm(p.model2, covariates=list(age=0,nahiv1=0), ci="normal", cl=0.95, B=10, cores=3)

#acquisition and clearance probability in all children from HH with HIV+ adult(s)
set.seed(1988)
pmatrix.msm(p.model2, covariates=list(age=0,nahiv1=1), ci="normal", cl=0.95, B=10, cores=3)


#---------------average period in days in a single stay in a state (sojourn times)
#carriage duration for HIV+ and HIV- children and adults from all HH
sojourn.msm(p.model2,covariates=list(hiv=1,age=0),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=0,age=0),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=1,age=1),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=0,age=1),ci="normal",cl=0.95,B=10,cores=3)

#carriage duration for HIV+ and HIV- children or adults from HH without HIV+ adult(s)
sojourn.msm(p.model2,covariates=list(hiv=1,age=0,nahiv1=0),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=0,age=0,nahiv1=0),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=1,age=1,nahiv1=0),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=0,age=1,nahiv1=0),ci="normal",cl=0.95,B=10,cores=3)

#carriage duration for HIV+ and HIV- children or adults from HH with HIV+ adult(s)
sojourn.msm(p.model2,covariates=list(hiv=1,age=0,nahiv1=1),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=0,age=0,nahiv1=1),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=1,age=1,nahiv1=1),ci="normal",cl=0.95,B=10,cores=3)
sojourn.msm(p.model2,covariates=list(hiv=0,age=1,nahiv1=1),ci="normal",cl=0.95,B=10,cores=3)

#---------------forecasted total length of time spent in each trasient state
totlos.msm(p.model2,fromt=0,tot=289,covariates=list(hiv=1,age=0),ci="normal",B=5,cl=0.95,cores=3)
totlos.msm(p.model2,fromt=0,tot=289,covariates=list(hiv=0,age=0),ci="normal",B=5,cl=0.95,cores=3)
totlos.msm(p.model2,fromt=0,tot=120,covariates=list(hiv=1,age=1),ci="normal",B=5,cl=0.95,cores=3)
totlos.msm(p.model2,fromt=0,tot=289,covariates=list(hiv=0,age=1),ci="normal",B=5,cl=0.95,cores=3)

#---------------expected time until Markov process first enters a carrying state (hitting time)
efpt.msm(p.model2,tostate=2,covariates=list(hiv=1,age=0),ci="normal",B=5,cl=0.95,cores=3)
efpt.msm(p.model2,tostate=2,covariates=list(hiv=0,age=0),ci="normal",B=5,cl=0.95,cores=3)
efpt.msm(p.model2,tostate=2,covariates=list(hiv=1,age=0),ci="normal",B=5,cl=0.95,cores=3)
efpt.msm(p.model2,tostate=2,covariates=list(hiv=0,age=0),ci="normal",B=5,cl=0.95,cores=3)

#---------------expected number of visits to a state
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=1, age=0),ci="normal",B=5,cl=0.95,cores=3)
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=0, age=0),ci="normal",B=5,cl=0.95,cores=3)
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=1, age=1),ci="normal",B=5,cl=0.95,cores=3)
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=0, age=1),ci="normal",B=5,cl=0.95,cores=3)

#---------------viterbi algorithm
phirst.vi <- viterbi.msm(p.model2)
phirst.vi$time <- as.integer(phirst.vi$time)
phirst.vi$observed <- as.integer(phirst.vi$observed)
phirst.vi$fitted <- as.integer(phirst.vi$fitted)
pstate <- as.data.frame(phirst.vi$pstate)
phirst.vi$probhs1 <- pstate$V1
phirst.vi$probhs2 <- pstate$V2
remove(pstate); phirst.vi$pstate <- NULL
phirst.vi$observed2 <- if_else(phirst.vi$observed==1L,"Clear","Carry")
phirst.vi$fitted2 <- if_else(phirst.vi$fitted==1L,"Clear","Carry")

dev.off()
hmm.plot0 <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_point(aes(time,observed2), color='red', size=2) + 
  theme_bw() +
  ylab("Observed") +
  xlab("") +
  theme(strip.text.x=element_text(size=12, face="bold", color="black")) +
  theme(axis.text.x=element_text(face="bold", size=12), axis.text.y=element_text(face="bold", size=12),axis.title=element_text(size=14)) 

hmm.plot1 <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_point(aes(time,fitted2), color='red', size=2) + 
  theme_bw() +
  ylab("Fitted (Viterbi)") +
  xlab("") +
  theme(strip.text.x = element_text(size=12, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=12), axis.text.y = element_text(face="bold", size=12),axis.title=element_text(size=14)) 

hmm.plot2 <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_line(aes(time,probhs2), color='blue', size=1.2) +
  theme_bw() +
  ylab("Probability (For/Back)") +
  scale_y_continuous(labels=percent) +
  xlab("Time (days)") +
  theme(strip.text.x = element_text(size=12, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=12), axis.text.y = element_text(face="bold", size=12),axis.title=element_text(size=14)) 

grid.arrange(grobs=list(hmm.plot0,hmm.plot1,hmm.plot2),ncol=1,nrow=3)

#---------------plot the transition intensities of a fitted HM model2
hmm.plot3 <-read.csv(curl("https://raw.githubusercontent.com/deusthindwa/markov.chain.model.pneumococcus.hiv.rsa/master/data/hmm_acq_plots.csv"))
dev.off()
hmm.plot3$Age = factor(hmm.plot3$Age,levels(hmm.plot3$Age)[c(2,4,1,3)])
ggplot(hmm.plot3, aes(x=Age, y=Intensity, color=Age)) + 
  geom_point(size=2.5,position=position_dodge(width=0.3),stat="identity") +
  geom_errorbar(aes(ymin=Lintensity, ymax=Uintensity), width=0.2,size=1,position=position_dodge(width=0.3),stat="identity") +
  facet_grid(.~HHIVstatus, scales="free_y") +
  ylim(c(10,25)) +
  theme_bw() +
  ylab("Acquisition per year") +
  xlab("") +
  theme(strip.text.x=element_text(size=12,face="bold",color="black")) +
  theme(axis.text.x=element_text(face="bold",size=0), axis.text.y=element_text(face="bold",size=12)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
  guides(color=guide_legend(title="")) +
  theme(legend.text=element_text(size=12), legend.title=element_text(face="bold",size=12)) 

#---------------plot the transition intensities of a fitted HM model2
hmm.plot4 <-read.csv(curl("https://raw.githubusercontent.com/deusthindwa/markov.chain.model.pneumococcus.hiv.rsa/master/data/hmm_dur_plots.csv"))
dev.off()
hmm.plot4$Age = factor(hmm.plot4$Age,levels(hmm.plot4$Age)[c(2,4,1,3)])
ggplot(hmm.plot4, aes(x=Age, y=Intensity, color=Age)) + 
  geom_point(size=2.5,position=position_dodge(width=0.3),stat="identity") +
  geom_errorbar(aes(ymin=Lintensity, ymax=Uintensity), width=0.2,size=1,position=position_dodge(width=0.3),stat="identity") +
  facet_grid(.~HHIVstatus, scales="free_y") +
  ylim(c(0,100)) +
  theme_bw() +
  ylab("Duration of carriage (days)") +
  xlab("") +
  theme(strip.text.x=element_text(size=12,face="bold",color="black")) +
  theme(axis.text.x=element_text(face="bold",size=0), axis.text.y=element_text(face="bold",size=12)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=15)) +
  guides(color=guide_legend(title="")) +
  theme(legend.text=element_text(size=12), legend.title=element_text(face="bold",size=12)) 


#===============hidden Markov modelling with transmission assumptions

#---------------define possible sequence of observed carriage data
phirst.tx <- phirst.fu
phirst.tx$visitno <- as.integer(substr(phirst.tx$visit_id,10,12))
phirst.tx$hstate <- if_else(phirst.tx$state==1 & phirst.tx$age==0 & phirst.tx$hiv==0,1,
                            if_else(phirst.tx$state==2 & phirst.tx$age==0 & phirst.tx$hiv==0,2,
                                    if_else(phirst.tx$state==1 & phirst.tx$age==1 & phirst.tx$hiv==0,3,
                                            if_else(phirst.tx$state==2 & phirst.tx$age==1 & phirst.tx$hiv==0,4,
                                                    if_else(phirst.tx$state==1 & phirst.tx$age==0 & phirst.tx$hiv==1,5,
                                                            if_else(phirst.tx$state==2 & phirst.tx$age==0 & phirst.tx$hiv==1,6,
                                                                    if_else(phirst.tx$state==1 & phirst.tx$age==1 & phirst.tx$hiv==1,7,
                                                                            if_else(phirst.tx$state==2 & phirst.tx$age==1 & phirst.tx$hiv==1,8,NULL))))))))

phirst.tx <- subset(phirst.tx, select=c(hh_id,visitno,hstate))
phirst.tx <- reshape(phirst.tx, idvar=c("hh_id","visitno"), timevar="hstate",v.names="hstate", direction="wide")
phirst.tx <- subset(phirst.tx, select=c(hh_id,visitno,hstate.1,hstate.2,hstate.3,hstate.4,hstate.5,hstate.6,hstate.7,hstate.8))

#---------------generate single household state per unit time based on observed carriage sequence and transmission assumption
phirst.tx <- phirst.tx %>%
  mutate(phirst.tx, state=if_else(!is.na(hstate.1) & is.na(hstate.2) & !is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),1,
                                  if_else(!is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),2,
                                          if_else(!is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & !is.na(hstate.8),3,
                                                  if_else(is.na(hstate.1) & !is.na(hstate.2) & !is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),4,
                                                          if_else(is.na(hstate.1) & !is.na(hstate.2) & is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),5,
                                                                  if_else(is.na(hstate.1) & !is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & !is.na(hstate.8),6,
                                                                          if_else(is.na(hstate.1) & is.na(hstate.2) & !is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),7,
                                                                                  if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),8,
                                                                                          if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & !is.na(hstate.8),9,
                                                                                                  
                                                                                                  if_else(!is.na(hstate.1) & !is.na(hstate.2) & !is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),5,
                                                                                                          if_else(!is.na(hstate.1) & !is.na(hstate.2) & !is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),4,
                                                                                                                  if_else(!is.na(hstate.1) & is.na(hstate.2) & !is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),2,
                                                                                                                          if_else(!is.na(hstate.1) & !is.na(hstate.2) & is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),5,
                                                                                                                                  if_else(is.na(hstate.1) & !is.na(hstate.2) & !is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),5,
                                                                                                                                          
                                                                                                                                          if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & !is.na(hstate.5) & !is.na(hstate.6) & !is.na(hstate.7) & !is.na(hstate.8),9,
                                                                                                                                                  if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & !is.na(hstate.5) & !is.na(hstate.6) & !is.na(hstate.7) & is.na(hstate.8),7,
                                                                                                                                                          if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & !is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & !is.na(hstate.8),9,
                                                                                                                                                                  if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & !is.na(hstate.6) & !is.na(hstate.7) & !is.na(hstate.8),9,
                                                                                                                                                                          
                                                                                                                                                                          if_else(!is.na(hstate.1) & !is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & !is.na(hstate.7) & !is.na(hstate.8),6,
                                                                                                                                                                                  if_else(!is.na(hstate.1) & !is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & is.na(hstate.7) & !is.na(hstate.8),6,
                                                                                                                                                                                          if_else(!is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & !is.na(hstate.7) & !is.na(hstate.8),3,
                                                                                                                                                                                                  if_else(is.na(hstate.1) & !is.na(hstate.2) & is.na(hstate.3) & is.na(hstate.4) & is.na(hstate.5) & is.na(hstate.6) & !is.na(hstate.7) & !is.na(hstate.8),6,
                                                                                                                                                                                                          
                                                                                                                                                                                                          if_else(is.na(hstate.1) & is.na(hstate.2) & !is.na(hstate.3) & !is.na(hstate.4) & !is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),8,
                                                                                                                                                                                                                  if_else(is.na(hstate.1) & is.na(hstate.2) & !is.na(hstate.3) & !is.na(hstate.4) & is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),8,
                                                                                                                                                                                                                          if_else(is.na(hstate.1) & is.na(hstate.2) & is.na(hstate.3) & !is.na(hstate.4) & !is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),8,
                                                                                                                                                                                                                                  if_else(is.na(hstate.1) & is.na(hstate.2) & !is.na(hstate.3) & is.na(hstate.4) & !is.na(hstate.5) & !is.na(hstate.6) & is.na(hstate.7) & is.na(hstate.8),7,NULL)))))))))))))))))))))))))))


#---------------show transition frequency in state table
statetable.msm(state,hh_id,data=phirst.tx)

#---------------initiate transition intensity matrix Q
matrix.Q.tx <- rbind(c(0.0,0.1,0.1,0.1,0.0,0.0,0.0),
                     c(0.1,0.0,0.0,0.0,0.1,0.0,0.0),
                     c(0.1,0.0,0.0,0.0,0.0,0.1,0.0),
                     c(0.1,0.0,0.0,0.0,0.1,0.1,0.0),
                     c(0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                     c(0.0,0.0,0.0,0.0,0.0,0.0,0.0),
                     c(0.0,0.0,0.0,0.0,0.0,0.0,0.0))
rownames(matrix.Q.tx) <- c("-C-,-A-","-C-,+A-","-C-,+A+","+C-,-A-","+C-,+A-","+C-,+A+","+C+,-A-")
colnames(matrix.Q.tx) <- c("-C-,-A-","-C-,+A-","-C-,+A+","+C-,-A-","+C-,+A-","+C-,+A+","+C+,-A-")

#---------------initiate hidden Markov matrix E
matrix.E.tx <- rbind(c(0.6,0.1,0.1,0.1,0.0,0.0,0.0),
                     c(0.1,0.8,0.0,0.0,0.1,0.0,0.0),
                     c(0.1,0.0,0.8,0.0,0.0,0.1,0.0),
                     c(0.1,0.0,0.0,0.7,0.1,0.1,0.0),
                     c(0.0,0.0,0.0,0.0,1.0,0.0,0.0),
                     c(0.0,0.0,0.0,0.0,0.0,1.0,0.0),
                     c(0.0,0.0,0.0,0.0,0.0,0.0,0.9))
rownames(matrix.E.tx) <- c("-C-,-A-","-C-,+A-","-C-,+A+","+C-,-A-","+C-,+A-","+C-,+A+","+C+,-A-")
colnames(matrix.E.tx) <- c("-C-,-A-","-C-,+A-","-C-,+A+","+C-,-A-","+C-,+A-","+C-,+A+","+C+,-A-")

#---------------fitting a transmission model without misclassification
phirst.tx <- arrange(phirst.tx,hh_id,visitno)
p.model3<-msm(state~visitno, subject=hh_id, data=phirst.tx, 
              qmatrix=matrix.Q.tx,
              ematrix=matrix.E.tx,
              est.initprobs=T,
              opt.method="bobyqa", control=list(maxfun=100000))
printnew.msm(p.model3)
save(p.model3)

#---------------model diagnostics
logLik.msm(p.model3, by.subject=TRUE)
