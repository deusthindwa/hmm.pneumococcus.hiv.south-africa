#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020


#===============LOAD PACKAGES AND DATASETS IN MEMORY

phirst.packages <- c("tidyverse","plyr","dplyr","msm","timetk","gridExtra","curl","minqa","table1",
                    "lubridate","magrittr","data.table","parallel","foreign","readstata13","ggpubr",
                    "wakefield","zoo","janitor","rethinking","scales","msmtools","nlme")
lapply(phirst.packages, library, character.only=TRUE)

#load all phirst datasets (household-level, individual-level, follow-up & antibiotic use)
phirst.hh <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_household.dta",generate.factors=T)
phirst.ms <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_master2.dta",generate.factors=T)
phirst.fu <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_follow-up1.dta",generate.factors=T)
phirst.ax <- read.dta13("~/Rproject/Markov.Model.Resources/data/phirst_follow-up2.dta",generate.factors=T)

#subset to get required variables in household-level dataset
phirst.hh <- subset(subset(phirst.hh, is.na(reason_hh_not_inc)),select=c(hh_id,hh_mems_11swabs))

#subset to get required variables in individual-level dataset
phirst.ms <- subset(subset(phirst.ms, ind_elig_inc=="Yes"), select=c(hh_id,ind_id,year,site,sex,age_at_consent,
                    hiv_status,arv_current_self,arv_current_vl,cd4_count_1,cd4_count_2,cd4_count_3,cd4_count_4,
                    cd4_count_5,cd4_count_6,cd4_count_7,cd4_count_8,pcv6wks,pcv14wks,pcv9mth,alcohol,anysmokenow))

#merge the household-level to individual-level dataset
phirst.ms <- merge(x=phirst.ms,y=phirst.hh,by="hh_id",all.y=TRUE)


#===============INDIVIDUAL-LEVEL DATASET DESCRIPTION

#rename analysis variables
phirst.ms <- rename(phirst.ms, c("age_at_consent"="age","hiv_status"="hiv","arv_current_vl"="artv","arv_current_self"="artr","anysmokenow"="smoke","hh_mems_11swabs"="hhsize"))

#age category
phirst.ms$agecat <- as.factor(if_else(phirst.ms$age<5,"Child",if_else(phirst.ms$age>=5,"Adult",NULL)))

#study site
phirst.ms$site <- as.factor(phirst.ms$site)

#sex
phirst.ms$sex <- as.factor(phirst.ms$sex)

#hiv status
phirst.ms$hiv <- as.factor(if_else(phirst.ms$hiv=="Negative","Neg",if_else(phirst.ms$hiv=="Positive","Pos",NULL)))

#art status
phirst.ms$artv <- as.factor(if_else(phirst.ms$artv=="Yes","Yes",if_else(phirst.ms$artv=="No","No",NULL)))
phirst.ms$artr <- as.factor(if_else(phirst.ms$artr=="Yes","Yes",if_else(phirst.ms$artr=="No","No",NULL)))

#cd4+ cell count
phirst.ms$cd4 <- rowMeans(cbind(phirst.ms$cd4_count_1,phirst.ms$cd4_count_2,phirst.ms$cd4_count_3,phirst.ms$cd4_count_4,
                 phirst.ms$cd4_count_5,phirst.ms$cd4_count_6,phirst.ms$cd4_count_7,phirst.ms$cd4_count_8),na.rm=TRUE)
phirst.ms$cd4 <- if_else(phirst.ms$hiv=="Pos",as.numeric(phirst.ms$cd4),NULL)
phirst.ms$cd4cat <- if_else(phirst.ms$cd4<=350 & phirst.ms$agecat=="Adult","Low",if_else(phirst.ms$cd4>350 & phirst.ms$agecat=="Adult","High",
                    if_else(phirst.ms$cd4<=750 & phirst.ms$agecat=="Child","Low",if_else(phirst.ms$cd4>750 & phirst.ms$agecat=="Child","High",NULL))))

#pcv status
phirst.ms$pcv6wks <- as.factor(if_else(phirst.ms$pcv6wks=="Yes" & phirst.ms$agecat=="Child","Yes",if_else(phirst.ms$pcv6wks=="No" & phirst.ms$agecat=="Child","No",NULL)))
phirst.ms$pcv14wks <- as.factor(if_else(phirst.ms$pcv14wks=="Yes" & phirst.ms$agecat=="Child","Yes",if_else(phirst.ms$pcv14wks=="No" & phirst.ms$agecat=="Child","No",NULL)))
phirst.ms$pcv9mth <- as.factor(if_else(phirst.ms$pcv9mth=="Yes" & phirst.ms$agecat=="Child","Yes",if_else(phirst.ms$pcv9mth=="No" & phirst.ms$agecat=="Child","No",NULL)))

#alcohol consumption
phirst.ms$alcohol <- as.factor(if_else(phirst.ms$alcohol==1 & phirst.ms$agecat=="Adult","Yes",if_else(phirst.ms$alcohol==0 & phirst.ms$agecat=="Adult","No",NULL)))

#smoking status
phirst.ms$smoke <- as.factor(if_else(phirst.ms$smoke==1 & phirst.ms$agecat=="Adult","Yes",if_else(phirst.ms$smoke==0 & phirst.ms$agecat=="Adult","No",NULL)))

#household size
phirst.ms$hhsize <- as.integer(phirst.ms$hhsize)

#number of hiv+ adults in the household
phirst.hi <- subset(subset(subset(phirst.ms,select=c(hh_id,agecat,hiv)),agecat !="Child"),!is.na(hiv))
phirst.hi$ahiv <- if_else(phirst.hi$age=="Adult" & phirst.hi$hiv=="Pos",1L,if_else(phirst.hi$age=="Adult" & phirst.hi$hiv=="Neg",0L,NULL))
phirst.hi <- setDT(phirst.hi)[,list(ahiv=sum(ahiv)),by=.(hh_id)]
phirst.hi$ahivcat <- as.factor(if_else(phirst.hi$ahiv==0,"No","Yes"))

#final individual-level dataset
phirst.ms <- merge(x=phirst.ms,y=phirst.hi,by="hh_id",all.x=TRUE)
phirst.ms <- subset(phirst.ms,select=c(hh_id,ind_id,hhsize,age,agecat,site,sex,hiv,artv,artr,cd4,cd4cat,pcv6wks,pcv14wks,pcv9mth,alcohol,smoke,ahiv,ahivcat))

#baseline demographic characteristics (Table 1)
table1(~age+cd4|agecat,data=phirst.ms,topclass="Rtable1-grid Rtable1-shade Rtable1-times")
phirst.render <- function(x,name,...){if(!is.numeric(x)) return(render.categorical.default(na.omit(x)))}
table1(~site+sex+hiv+artv+artr+cd4cat+ahivcat+pcv6wks+pcv14wks+pcv9mth+smoke+alcohol|agecat, data=phirst.ms,
       topclass="Rtable1-grid Rtable1-shade Rtable1-times",render=phirst.render)


#===============FOLLOW-UP DATASET DESCRIPTION

#subset individual-level dataset that will merge follow-up dataset
phirst.ms <- arrange(subset(subset(phirst.ms,select=c(ind_id,hhsize,agecat,hiv,artv,artr,ahiv,ahivcat)),!is.na(hiv)),ind_id)

#subset antibiotic dataset that will merge follow-up dataset
phirst.ax <- arrange(rename(subset(phirst.ax,select=c(finalid,antibiotic)),c("finalid"="visit_id","antibiotic"="abx")),visit_id)
phirst.ax$abx <- as.factor(if_else(phirst.ax$abx==1,"Yes","No"))

#subset follow-up dataset that will merge individual-level & antibiotic datasets
phirst.fu <- arrange(subset(phirst.fu, select=c(visit_id,ind_id,visit_date,visit,npspne,npspneload)),visit_id)
phirst.fu$dys <- as.integer(difftime(ymd(phirst.fu$visit_date),na.locf(if_else(phirst.fu$visit==1L,phirst.fu$visit_date,NULL)),units="days"))
phirst.fu <- rename(subset(phirst.fu,select=c(visit_id,ind_id,dys,npspne,npspneload)),c("npspne"="state","npspneload"="npdensity"))

#merge follow-up, antibiotic and individual-level datasets
phirst.fu <- arrange(merge(x=phirst.fu,y=phirst.ax,by="visit_id",all.x=TRUE),visit_id)
phirst.fu <- arrange(merge(x=phirst.fu,y=phirst.ms,by="ind_id",all.y=TRUE),visit_id)

#tidying the follow-up dataset
phirst.fu$state <- if_else(phirst.fu$state==1,2L,if_else(phirst.fu$state==0,1L,NULL));phirst.fu$state[is.na(phirst.fu$state)]<-9L
phirst.fu <- mutate(phirst.fu,hh_id=substr(visit_id,1,4),visit_no=as.integer(substr(visit_id,10,12)))
phirst.fu <- phirst.fu %>% mutate(obst=if_else(phirst.fu$state==9L,1L,0L)) %>% select(visit_id,ind_id,hh_id,visit_no,dys,state,obst,npdensity,abx,hhsize,agecat,hiv,artv,artr,ahiv,ahivcat)

#define community or household acquisition source
phirst.tx <- arrange(subset(subset(phirst.fu,select=c(hh_id,visit_no,state)),state !=9),visit_no)
phirst.tx$state <- if_else(phirst.tx$state==1,0L,1L)
phirst.tx <- setDT(phirst.tx)[,list(tx=sum(state)),by=.(hh_id,visit_no)]
phirst.tx$tx <- as.factor(if_else(phirst.tx$tx>=1,"hhtx","cmtx"))
phirst.fu <- subset(merge(x=phirst.fu, y=phirst.tx, by=c("hh_id","visit_no"),all.y=TRUE),select=c(visit_id,dys,state,obst,npdensity,abx,hhsize,agecat,hiv,artv,artr,ahiv,ahivcat,tx))

#follow-up characteristics of carriage among participants (figure 2)
source('~/Rproject/Markov.Model/script/fig2.R')


#===============hidden Markov modelling of carriage dynamics within houshold and from community

#show transition frequency
phirst.fu$abx[is.na(phirst.fu$abx)] <- phirst.fu$artv[is.na(phirst.fu$artv)] <- phirst.fu$artr[is.na(phirst.fu$artr)] <- "No"
phirst.fu$ind_id <- substr(phirst.fu$visit_id,1,8)
phirst.fu <- arrange(phirst.fu,visit_id)
statetable.msm(state,ind_id,data=phirst.fu)

#initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.9),c(0.9,0.0))
rownames(matrix.Q) <- c("Clear","Carry")
colnames(matrix.Q) <- c("Clear","Carry")

#initiate emission matrix E
matrix.E <- rbind(c(1.0,0.0),c(0.1,0.9))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carry")

#fit hidden Markov models with misclassifications
p.model1 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model2 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model3 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model4 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

#comparing the AIC or BIC of the fitted models 
AIC(p.model1,p.model2,p.model3,p.model4)
AIC(p.model1,k=log(length(phirst.fu)));AIC(p.model2,k=log(length(phirst.fu)));AIC(p.model3,k=log(length(phirst.fu)));AIC(p.model4,k=log(length(phirst.fu)))

#run multiple chains to assess convergence of the selected model
j=0.05;k=2.00
for(i in 1:5){
  
matrix.Qc <- rbind(c(0.0,j), c(k,0.0))
rownames(matrix.Qc) <- c("Clear","Carry")
colnames(matrix.Qc) <- c("Clear","Carry")

matrix.Ec <- rbind(c(1.0,0.0), c(0.15,0.85))
colnames(matrix.Ec) <- c("SwabNeg","SwabPos")
rownames(matrix.Ec) <- c("Clear","Carry")

sink("~/Rproject/Markov.Model.Resources/data/convergence.txt",append=TRUE)
p.convrg <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Qc, ematrix=matrix.Ec,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                control=list(maxit=200000,trace=1,REPORT=1))
sink()
j=j+0.1;k=k-0.4
}

#observed versus predicted prevalence of selected model
p.obsexp <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

phirst.oe <- tk_tbl(prevalence.msm(p.obsexp,times=seq(0,289,14)),preserve_index=TRUE,rename_index="Time")
phirst.oe <- subset(phirst.oe,select=c(Time, Observed.State.1, Observed.State.2, Observed.percentages.State.1, Observed.percentages.State.2,
                                  Expected.Clear,Expected.Carry,Expected.percentages.Clear,Expected.percentages.Carry))
phirst.oe <- rename(phirst.oe, c("Observed.State.1"="obs.clear", "Observed.State.2"="obs.carry","Observed.percentages.State.1"="obs.p.clear","Observed.percentages.State.2"="obs.p.carry",
                           "Expected.Clear"="exp.clear","Expected.Carry"="exp.carry","Expected.percentages.Clear"="exp.p.clear","Expected.percentages.Carry"="exp.p.carry"))

phirst.oe$lci.clear=phirst.oe$exp.p.clear/100-(1.96*sqrt(phirst.oe$exp.p.clear/100*(1-phirst.oe$exp.p.clear/100)/phirst.oe$exp.clear)) 
phirst.oe$uci.clear=phirst.oe$exp.p.clear/100+(1.96*sqrt(phirst.oe$exp.p.clear/100*(1-phirst.oe$exp.p.clear/100)/phirst.oe$exp.clear))
phirst.oe$lci.carry=phirst.oe$exp.p.carry/100-(1.96*sqrt(phirst.oe$exp.p.carry/100*(1-phirst.oe$exp.p.carry/100)/phirst.oe$exp.carry)) 
phirst.oe$uci.carry=phirst.oe$exp.p.carry/100+(1.96*sqrt(phirst.oe$exp.p.carry/100*(1-phirst.oe$exp.p.carry/100)/phirst.oe$exp.carry))

#plot of model parameter convergence, and observed and predictions (supplementary figure 1)
dev.off()
source('~/Rproject/Markov.Model/script/sfig1.R')

#plot results from Viterbi algorithm (supplementary figure 2)
phirst.vi <- viterbi.msm(p.model4)
phirst.vi$fitted <- as.integer(phirst.vi$fitted)
phirst.vi$probhs1 <- as.data.frame(phirst.vi$pstate)$V1
phirst.vi$probhs2 <- as.data.frame(phirst.vi$pstate)$V2
phirst.vi$observed <- if_else(phirst.vi$observed==1L,"Clear","Carry")
phirst.vi$fitted <- if_else(phirst.vi$fitted==1L,"Clear","Carry")
dev.off()
source('~/Rproject/Markov.Model/script/sfig2.R')

#plot of within household and community acquisition rates and probabilities (figure 3)
dev.off()
source('~/Rproject/Markov.Model/script/fig3.R')

#plot of duration of carriage and carriage clearance probabilities (figure 4)
dev.off()
source('~/Rproject/Markov.Model/script/fig4.R')

##plot of sensitivity analysis of # of adults HH HIV+ & sampling times (supplementary figure 3)
dev.off()
source('~/Rproject/Markov.Model/script/sfig3.R')

#plot of acquistion rates by household size (supplementary figure 4)
dev.off()
source('~/Rproject/Markov.Model/script/sfig4.R')
