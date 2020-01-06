#Written by Deus Thindwa
#The effect of parental HIV infection status on pneumococcal acquisition and clearance rates in children
#Continuous-time time-inhomogeneous hidden Markov modelling study, PhD chapter 1.
#13/08/2019

#---------------load required packages into memory
phirst.packages <-c("tidyverse","plyr","msm","timetk","gridExtra","curl","dplyr","minqa","lubridate","magrittr","data.table","parallel","foreign","readstata13","wakefield","zoo","janitor","rethinking","doParallel")
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

#---------------prepare dataset for baseline demographic characteristics
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
                                if_else(phirst.ms$arv_current=="No","No",NULL))

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

#final master dataset
phirst.ms <- subset(phirst.ms, select=c(hh_id,ind_id,year,hhsize,age,site,sex,hiv,art,cd4,vl,pcv6w,pcv14w,pcv9m,alcohol,smoke))

#---------------baseline demographic characteristics (Table 1)
phirst.ms %>% tabyl(age, show_na=FALSE) %>% adorn_pct_formatting(digits=1)

for(i in colnames(phirst.ms[c(5:16)])){
  print(tabyl(phirst.ms[[i]],show_na=FALSE) %>% adorn_pct_formatting(digits=1))
};remove(i)

for(i in colnames(phirst.ms[c(6:16)])){
  print(tabyl(phirst.ms[[i]][phirst.ms$age=="Child"],show_na=FALSE) %>% adorn_pct_formatting(digits=1))
};remove(i)

for(i in colnames(phirst.ms[c(6:16)])){
  print(tabyl(phirst.ms[[i]][phirst.ms$age=="Adult"],show_na=FALSE) %>% adorn_pct_formatting(digits=1))
};remove(i)

#---------------prepare and merge follow-up to master dataset
phirst.fu <- subset(phirst.fu, select=c(ind_id,visit_id,visit_date,visit,npspne,npspneload))
phirst.fu$visit_date <- ymd(phirst.fu$visit_date)
phirst.fu$startdate <- if_else(phirst.fu$visit==1L,phirst.fu$visit_date,NULL)
phirst.fu$startdate <- na.locf(phirst.fu$startdate)
phirst.fu$dys <- difftime(phirst.fu$visit_date,phirst.fu$startdate,units="days")
phirst.fu$pndensity <- if_else(phirst.fu$npspneload>=1000,"High",
                               if_else(phirst.fu$npspneload<1000,"Low",NULL))
phirst.fu <- merge(phirst.fu,phirst.ms)

#---------------final follow-up dataset
phirst.fu <- subset(phirst.fu, select=c(hh_id,ind_id,visit_id,hhsize,dys,npspne,pndensity,age,sex,hiv,art,cd4,vl))

#---------------prepare follow-up dataset for hidden Markov model fitting
phirst.fu$dys <- as.integer(phirst.fu$dys)
phirst.fu <- rename(phirst.fu, c("npspne" = "state"))
phirst.fu$state <- if_else(phirst.fu$state==0,1L,
                          if_else(phirst.fu$state==1,2L,NULL))
phirst.fu$pndensity <- recode_factor(phirst.fu$pndensity,`High`=1L,`Low`=0L)
phirst.fu$age <- recode_factor(phirst.fu$age,`Adult`=1L,`Child`=0L)
phirst.fu$sex <- recode_factor(phirst.fu$sex,`Male`=1L,`Female`=0L)
phirst.fu$hiv <- recode_factor(phirst.fu$hiv,`Positive`=1L,`Negative`=0L)
phirst.fu$art <- recode_factor(phirst.fu$art,`Yes`=1L,`No`=0L)
phirst.fu$cd4 <- recode_factor(phirst.fu$cd4,`High`=1L,`Low`=0L)
phirst.fu$vl <- recode_factor(phirst.fu$vl,`High`=1L,`Low`=0L)

#===============hidden Markov modelling without transmission assumptions

#---------------show transition frequency in state table
statetable.msm(state,ind_id,data=phirst.fu)

#---------------initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.1),
                  c(0.1,0.0))
rownames(matrix.Q) <- c("Clear","Carry")
colnames(matrix.Q) <- c("Clear","Carry")

#---------------fitting a time-homogeneous model without misclassification
phirst.fu <- arrange(phirst.fu,ind_id,dys)
p.model0<-msm(state~dys, subject=ind_id, data=phirst.fu, 
              qmatrix=matrix.Q, 
              covariates=~age+hiv,
              opt.method="bobyqa")

printnew.msm(p.model0)

#---------------fitting a time-inhomogeneous model without misclassification
p.model1<-msm(state~dys, subject=ind_id, data=phirst.fu, 
              qmatrix=matrix.Q, 
              covariates=~age+hiv,
              pci=c(100,150,200,225,275),
              opt.method="bobyqa")

printnew.msm(p.model1)

#---------------Likelihood ratio test
lrtest.msm(p.model0,p.model1)

#---------------pearson-type goodness of fit
options(digits=2)
pearson.msm(p.model0, timegroups=2)
pearson.msm(p.model1, timegroups=2)

#---------------expected vs observed carriage
dev.off()
par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
plot.prevalence.msm(p.model0, mintime=0,maxtime=288,legend.pos=c(0,100),lwd.obs=2,lwd.exp=2,cex=0.7,xlab="Time (days)",ylab="% Prevalence")
legend(0,100, legend=c("Observed carriage", "Fitted model"),col=c("blue","red"), lty=1:3, cex=1, lwd=3)

dev.off()
par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
plot.prevalence.msm(p.model1, mintime=0,maxtime=288,legend.pos=c(0,100),lwd.obs=2,lwd.exp=2,cex=0.7,xlab="Time (days)",ylab="% Prevalence")
legend(0,100, legend=c("Observed carriage", "Fitted model"),col=c("blue","red"), lty=1:3, cex=1, lwd=3)

#---------------fitting a simpler model without misclassification with covariates seperately
qmatrix.msm(p.model1, covariates="mean")

qmatrix.msm(p.model1, covariates=list(hiv=0))
qmatrix.msm(p.model1, covariates=list(hiv=1))

qmatrix.msm(p.model1, covariates=list(age=0))
qmatrix.msm(p.model1, covariates=list(age=1))

#---------------fitting a model without misclassification with covariates in combination
qmatrix.msm(p.model1, covariates=list(hiv=1, age=0))
qmatrix.msm(p.model1, covariates=list(hiv=0, age=0))
qmatrix.msm(p.model1, covariates=list(hiv=1, age=1))
qmatrix.msm(p.model1, covariates=list(hiv=0, age=1))

#---------------initiate emission matrix E
matrix.E <- rbind(c(1.0,0.0), 
                  c(0.1,0.9))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carry")

#---------------fitting a model with misclassification
p.model2 <- msm(state~dys, subject=ind_id, data=phirst.fu,
                qmatrix=matrix.Q,
                ematrix=matrix.E,
                covariates= ~age+hiv,
                pci=c(100,150,200,225,275),
                est.initprobs=T,
                opt.method="bobyqa")

printnew.msm(p.model2)

#---------------transition intensity matrix for a model with misclassification with covariates (Table 2)
qmatrix.msm(p.model2, covariates=list(hiv=1, age=0))
qmatrix.msm(p.model2, covariates=list(hiv=0, age=0))
qmatrix.msm(p.model2, covariates=list(hiv=1, age=1))
qmatrix.msm(p.model2, covariates=list(hiv=0, age=1))

#---------------transition probability matrix for model with misclassification with covariates (Table 3)
pmatrix.msm(p.model2, covariates=list(hiv=1, age=0),t1=0,t=289,ci="bootstrap",B=5,cl=0.95,cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0, age=0),t1=0,t=289,ci="bootstrap",B=5,cl=0.95,cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=1, age=1),t1=0,t=289,ci="bootstrap",B=5,cl=0.95,cores=3)
pmatrix.msm(p.model2, covariates=list(hiv=0, age=1),t1=0,t=289,ci="bootstrap",B=5,cl=0.95,cores=3)

#---------------pearson-type goodness of fit
#options(digits=2)
#pearson.msm(p.model2, timegroups=2)

#---------------expected vs observed carriage
#dev.off()
#par(mgp=c(2,1,0),mar=c(6,4,2,2)+0.1)
#plot.prevalence.msm(p.model2, mintime=0,maxtime=288,legend.pos=c(0,100),lwd.obs=2.5,lwd.exp=2.5,cex=0.7,xlab="Time (days)",ylab="% Prevalence")
#legend(0, 100, legend=c("Observed carriage", "Fitted model"),col=c("blue","red"), lty=1:3, cex=1, lwd=3)

#---------------average period in a single stay in a state (sojourn)
sojourn.msm(p.model2)

#---------------Probability that each state is next
pnext.msm(p.model2)

#---------------forecasted total length of time spent in each trasient state
totlos.msm(p.model2, tot=289)

#---------------expected time until Markov process first enters a carrying state (hitting time)
efpt.msm(p.model2, tostate=2)

#---------------expected number of visits to a state
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=1, age=0),ci="bootstrap",B=5,cl=0.95,cores=3)
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=0, age=0),ci="bootstrap",B=5,cl=0.95,cores=3)
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=1, age=1),ci="bootstrap",B=5,cl=0.95,cores=3)
envisits.msm(p.model2, fromt=0,tot=289,covariates=list(hiv=0, age=1),ci="bootstrap",B=5,cl=0.95,cores=3)

#---------------viterbi algorithm
phirst.vi <- viterbi.msm(p.model1)
phirst.vi$time <- as.integer(phirst.vi$time)
phirst.vi$observed <- as.integer(phirst.vi$observed)
phirst.vi$fitted <- as.integer(phirst.vi$fitted)
pstate <- as.data.frame(phirst.vi$pstate)
phirst.vi$probhs1 <- pstate$V1
phirst.vi$probhs2 <- pstate$V2
remove(pstate); phirst.vi$pstate <- NULL

dev.off()
hmm.plot1 <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_point(aes(time,fitted), color='red') + 
  scale_y_continuous(breaks=c(1L,2L)) +
  theme_bw() +
  ylab("Fitted states") +
  xlab("") +
  theme(strip.text.x = element_text(size=11, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=10)) 

hmm.plot2 <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_line(aes(time,probhs2), color='blue') +
  theme_bw() +
  ylab("Probability of hidden state") +
  xlab("Time (days)") +
  theme(strip.text.x = element_text(size=11, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=10)) 

  grid.arrange(grobs=list(hmm.plot1,hmm.plot2),ncol=1,nrow=2)

#---------------plot the transition intensities of a fitted HM model2
hmm.plot3 <-read.csv(curl("https://raw.githubusercontent.com/deusthindwa/markov.chain.model.pneumococcus.hiv.rsa/master/data/hmm_general_plots.csv"))
dev.off()
ggplot(hmm.plot3, aes(x=Age, y=Intensity*100, color=Age)) + 
  geom_point(size=2,position=position_dodge(width=0.3),stat="identity") +
  geom_errorbar(aes(ymin=Lintensity*100, ymax=Uintensity*100), width=0.2,size=1,position=position_dodge(width=0.3),stat="identity") +
  facet_grid(.~State, scales="free_y") +
  ylim(c(0,10)) +
  theme_bw() +
  ylab("Transition per 100 days") +
  xlab("") +
  theme(strip.text.x = element_text(size = 11, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=0), axis.text.y = element_text(face="bold", size=10)) + 
  guides(color=guide_legend(title="")) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(face="bold", size=11)) 

#---------------plot the transition probabilities of a fitted HM model2
hmm.plot4 <-read.csv(curl("https://raw.githubusercontent.com/deusthindwa/markov.chain.model.pneumococcus.hiv.rsa/master/data/hmm_prob_plots.csv"))
dev.off()
ggplot(hmm.plot4, aes(x=Age, y=Intensity*100, color=Age)) + 
  geom_point(size=2,position=position_dodge(width=0.3),stat="identity") +
  geom_errorbar(aes(ymin=Lintensity*100, ymax=Uintensity*100), width=0.2,size=1,position=position_dodge(width=0.3),stat="identity") +
  facet_grid(.~State, scales="free_y") +
  ylim(c(0,100)) +
  theme_bw() +
  ylab("Probability (%)") +
  xlab("") +
  theme(strip.text.x = element_text(size = 11, face="bold", color="black")) +
  theme(axis.text.x = element_text(face="bold", size=0), axis.text.y = element_text(face="bold", size=10)) + 
  guides(color=guide_legend(title="")) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(face="bold", size=11)) 

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
              opt.method="bobyqa", control=list(maxfun=250000))
printnew.msm(p.model3)
save(p.model3)

#---------------model diagnostics
logLik.msm(p.model3, by.subject=TRUE)
 