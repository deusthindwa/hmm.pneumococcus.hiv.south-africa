
#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020


#===============LOAD PACKAGES AND DATASETS IN MEMORY

phirst.packages <- c("tidyverse","plyr","msm","timetk","gridExtra","curl","minqa","table1",
                    "lubridate","magrittr","data.table","parallel","foreign","readstata13","ggpubr",
                    "wakefield","zoo","janitor","scales","msmtools","FSA","nlme","patchwork","boot")
lapply(phirst.packages, library, character.only=TRUE)

#load all phirst datasets (household-level, individual-level, follow-up & antibiotic use)
phirst.hh <- read.dta13("~/Rproject/Markov.Model/data/phirst_household.dta",generate.factors=T)
phirst.ms <- read.dta13("~/Rproject/Markov.Model/data/phirst_master2.dta",generate.factors=T)
phirst.fu <- read.dta13("~/Rproject/Markov.Model/data/phirst_follow-up1.dta",generate.factors=T)
phirst.ax <- read.dta13("~/Rproject/Markov.Model/data/phirst_follow-up2.dta",generate.factors=T)

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
phirst.ms$agecat <- as.factor(if_else(phirst.ms$age<5,"Younger child",if_else(phirst.ms$age>=5 & phirst.ms$age<18,"Older child",if_else(phirst.ms$age>=18,"Adult",NULL))))

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
                    if_else(phirst.ms$cd4<=750 & phirst.ms$agecat=="Younger child","Low",if_else(phirst.ms$cd4>750 & phirst.ms$agecat=="Younger child","High",
                    if_else(phirst.ms$cd4<=750 & phirst.ms$agecat=="Older child","Low",if_else(phirst.ms$cd4>750 & phirst.ms$agecat=="Older child","High",NULL))))))

#pcv status
phirst.ms$pcv6wks <- as.factor(if_else(phirst.ms$pcv6wks=="Yes" & phirst.ms$agecat=="Younger child","Yes",if_else(phirst.ms$pcv6wks=="No" & phirst.ms$agecat=="Younger child","No",NULL)))
phirst.ms$pcv14wks <- as.factor(if_else(phirst.ms$pcv14wks=="Yes" & phirst.ms$agecat=="Younger child","Yes",if_else(phirst.ms$pcv14wks=="No" & phirst.ms$agecat=="Younger child","No",NULL)))
phirst.ms$pcv9mth <- as.factor(if_else(phirst.ms$pcv9mth=="Yes" & phirst.ms$agecat=="Younger child","Yes",if_else(phirst.ms$pcv9mth=="No" & phirst.ms$agecat=="Younger child","No",NULL)))

#alcohol consumption
phirst.ms$alcohol <- as.factor(if_else(phirst.ms$alcohol==1 & phirst.ms$agecat=="Adult","Yes",if_else(phirst.ms$alcohol==0 & phirst.ms$agecat=="Adult","No",NULL)))

#smoking status
phirst.ms$smoke <- as.factor(if_else(phirst.ms$smoke==1 & phirst.ms$agecat=="Adult","Yes",if_else(phirst.ms$smoke==0 & phirst.ms$agecat=="Adult","No",NULL)))

#household size
phirst.ms$hhsize <- as.integer(phirst.ms$hhsize)
phirst.ms$hhsizecat <- as.factor(if_else(phirst.ms$hhsize<6,"<6",if_else(phirst.ms$hhsize>=6 & phirst.ms$hhsize<=10,"6-10","11+")))

#number of hiv+ adults in the household
phirst.hi <- subset(subset(subset(phirst.ms,select=c(hh_id,agecat,hiv)),agecat =="Adult"),!is.na(hiv))
phirst.hi$ahiv <- if_else(phirst.hi$agecat=="Adult" & phirst.hi$hiv=="Pos",1L,if_else(phirst.hi$agecat=="Adult" & phirst.hi$hiv=="Neg",0L,NULL))
phirst.hi <- setDT(phirst.hi)[,list(ahiv=sum(ahiv)),by=.(hh_id)]
phirst.hi$ahivcat <- as.factor(if_else(phirst.hi$ahiv==0,"No","Yes"))

#households with HIV+ (yes) and HIV- (no) female adults
phirst.sxf <- subset(subset(subset(phirst.ms,select=c(hh_id,sex,agecat,hiv)), sex=="Female" & agecat=="Adult"),!is.na(hiv))
phirst.sxf$sx <- if_else(phirst.sxf$hiv=="Neg",0L,1L)
phirst.sxf <- setDT(phirst.sxf)[,list(sxf=sum(sx)),by=.(hh_id)]
phirst.sxf$sxf <- as.factor(if_else(phirst.sxf$sxf==0,"No","Yes"))

#households with HIV+ (yes) and HIV- (no) male adults
phirst.sxm <- subset(subset(subset(phirst.ms,select=c(hh_id,sex,agecat,hiv)), sex=="Male" & agecat=="Adult"),!is.na(hiv))
phirst.sxm$sx <- if_else(phirst.sxm$hiv=="Neg",0L,1L)
phirst.sxm <- setDT(phirst.sxm)[,list(sxm=sum(sx)),by=.(hh_id)]
phirst.sxm$sxm <- as.factor(if_else(phirst.sxm$sxm==0,"No","Yes"))

#final individual-level dataset
phirst.ms <- merge(x=phirst.ms,y=phirst.hi,by="hh_id",all.x=TRUE)
phirst.ms <- merge(x=phirst.ms,y=phirst.sxf,by="hh_id",all.x=TRUE)
phirst.ms <- merge(x=phirst.ms,y=phirst.sxm,by="hh_id",all.x=TRUE)
phirst.ms <- subset(phirst.ms,select=c(hh_id,ind_id,year,hhsize,hhsizecat,age,agecat,site,sex,hiv,artv,artr,cd4,cd4cat,pcv6wks,pcv14wks,pcv9mth,alcohol,smoke,ahiv,ahivcat,sxf,sxm))

#baseline demographic characteristics (Table 1)
table1(~age+cd4|agecat,data=phirst.ms,topclass="Rtable1-grid Rtable1-shade Rtable1-times")
phirst.render <- function(x,name,...){if(!is.numeric(x)) return(render.categorical.default(na.omit(x)))}
table1(~site+sex+hiv+artv+artr+cd4cat+ahivcat+pcv6wks+pcv14wks+pcv9mth+smoke+alcohol|agecat, data=phirst.ms,
       topclass="Rtable1-grid Rtable1-shade Rtable1-times",render=phirst.render)


#===============FOLLOW-UP DATASET DESCRIPTION

#subset individual-level dataset that will merge follow-up dataset
phirst.ms <- arrange(subset(subset(phirst.ms,select=c(ind_id,year,site,hhsize,hhsizecat,age,sex,agecat,hiv,artv,artr,ahiv,ahivcat,sxf,sxm)),!is.na(hiv)),ind_id)

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
phirst.fu <- phirst.fu %>% mutate(obst=if_else(phirst.fu$state==9L,1L,0L)) %>% select(visit_id,ind_id,hh_id,visit_no,year,site,dys,state,obst,npdensity,abx,hhsize,hhsizecat,age,sex,agecat,hiv,artv,artr,ahiv,ahivcat,sxf,sxm)

#define community or household acquisition source
phirst.tx <- arrange(subset(subset(phirst.fu,select=c(hh_id,visit_no,state)),state !=9),visit_no)
phirst.tx$state <- if_else(phirst.tx$state==1,0L,1L)
phirst.tx <- setDT(phirst.tx)[,list(tx=sum(state)),by=.(hh_id,visit_no)]
phirst.tx$tx <- as.factor(if_else(phirst.tx$tx>=1,"hhtx","cmtx"))
phirst.fu <- subset(merge(x=phirst.fu, y=phirst.tx, by=c("hh_id","visit_no"),all.y=TRUE),select=c(visit_id,visit_no, hh_id,ind_id,year,site,dys,state,obst,npdensity,abx,hhsize,hhsizecat,age,sex,agecat,hiv,artv,artr,ahiv,ahivcat,sxf,sxm,tx))

#define infection contribution from other house members
phirst.fx <- arrange(subset(subset(phirst.fu,select=c(hh_id,visit_no,state)),state !=9),visit_no)
phirst.fx$state <- if_else(phirst.fx$state==1,0L,1L)
phirst.fx <- setDT(phirst.fx)[,list(fx=sum(state)),by=.(hh_id,visit_no)]
phirst.fu <- subset(merge(x=phirst.fu, y=phirst.fx, by=c("hh_id","visit_no"),all.y=TRUE),select=c(visit_id,ind_id,year,site,dys,state,obst,npdensity,abx,hhsize,hhsizecat,age,sex,agecat,hiv,artv,artr,ahiv,ahivcat,sxf,sxm,tx, fx))


#===============HIDDEN MARKOV MODEL WITH MISCLASSIFICATIONS

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
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat,"2-1"=~agecat+hiv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model2 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat,"2-1"=~agecat+hiv+abx),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model3 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat,"2-1"=~agecat+hiv+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

p.model4 <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

#accounting for infection probability from other household members 
p.modelfx <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+fx,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

#comparing the AIC scores of the fitted models 
AIC(p.model1,p.model2,p.model3,p.model4)

#follow-up carriage characteristics of participants (figure 2)
source('~/Rproject/Markov.Model/script/fig2_carriage_dynamics.R')

#plot within household and community acquisition probabilities (figure 3)
source('~/Rproject/Markov.Model/script/fig3_carriage_acquisition.R')

#plot duration of carriage (and by ART, ABX) and carriage clearance probabilities (figure 4)
source('~/Rproject/Markov.Model/script/fig4_carriaage_duration.R')

#plot model parameter convergence, and observed and predictions (supplementary figure 1)
source('~/Rproject/Markov.Model/script/figs1_convergence.R')

#plot results from Viterbi algorithm (supplementary figure 2)
source('~/Rproject/Markov.Model/script/figs2_viterbi.R')

##plot sensitivity analysis of # of HIV+ adults in a household & sampling times (supplementary figure 3)
source('~/Rproject/Markov.Model/script/figs3_sensitivity.R')

#plot acquistion probabilities by household size (supplementary figure 4)
source('~/Rproject/Markov.Model/script/figs4_householdsize.R')

#plot acquistion probabilities by sex (supplementary figure 5)
source('~/Rproject/Markov.Model/script/figs5_carriage_acquisition_sex.R')
