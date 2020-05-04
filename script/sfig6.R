#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 11/3/2020

#misclassification probabilities at different follow-up points
p.model4e1 <- msm(state~dys,subject=ind_id,data=subset(phirst.fu,dys>=0 & dys<=99),
qmatrix=matrix.Q, ematrix=matrix.E,
covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
opt.method="bobyqa", control=list(maxfun=100000))

p.model4e2 <- msm(state~dys,subject=ind_id,data=subset(phirst.fu,dys>=100 & dys<=139),
qmatrix=matrix.Q, ematrix=matrix.E,
covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
opt.method="bobyqa", control=list(maxfun=100000))

p.model4e3 <- msm(state~dys,subject=ind_id,data=subset(phirst.fu,dys>=140 & dys<=199),
qmatrix=matrix.Q, ematrix=matrix.E,
covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
opt.method="bobyqa", control=list(maxfun=100000))

p.model4e4 <- msm(state~dys,subject=ind_id,data=subset(phirst.fu,dys>=200 & dys<=219),
qmatrix=matrix.Q, ematrix=matrix.E,
covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
opt.method="bobyqa", control=list(maxfun=100000))

p.model4e5 <- msm(state~dys,subject=ind_id,data=subset(phirst.fu,dys>=220 & dys<=289),
qmatrix=matrix.Q, ematrix=matrix.E,
covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artv),
censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
opt.method="bobyqa", control=list(maxfun=100000))

#extract misclassification probabilities at different sampling periods
phirst.es <- data.frame("Period"=c("Overall","0-99","100-139","140-199","200-219","220-289"),
                        "Estimate"=c("NULL","NULL","NULL","NULL","NULL","NULL"), 
                        "Lbound"=c("NULL","NULL","NULL","NULL","NULL","NULL"),
                        "Ubound"=c("NULL","NULL","NULL","NULL","NULL","NULL"))

phirst.es$Estimate <- c(p.model4$Ematrices$baseline[2,1],p.model4e1$Ematrices$baseline[2,1],p.model4e2$Ematrices$baseline[2,1],
                        p.model4e3$Ematrices$baseline[2,1],p.model4e4$Ematrices$baseline[2,1],p.model4e5$Ematrices$baseline[2,1])
phirst.es$Lbound <- c(p.model4$EmatricesL$baseline[2,1],p.model4e1$EmatricesL$baseline[2,1],p.model4e2$EmatricesL$baseline[2,1],
                      p.model4e3$EmatricesL$baseline[2,1],p.model4e4$EmatricesL$baseline[2,1],p.model4e5$EmatricesL$baseline[2,1])
phirst.es$Ubound <- c(p.model4$EmatricesU$baseline[2,1],p.model4e1$EmatricesU$baseline[2,1],p.model4e2$EmatricesU$baseline[2,1],
                      p.model4e3$EmatricesU$baseline[2,1],p.model4e4$EmatricesU$baseline[2,1],p.model4e5$EmatricesU$baseline[2,1])

#plot misclassification probabilities at different sampling periods
ggplot(phirst.es) +
geom_errorbar(aes(Period,color=Period,ymin=Lbound,ymax=Ubound),width=0.1,size=0.8,position=position_dodge(width=0.5)) +
geom_point(aes(Period,Estimate,color=Period),size=2.5,position=position_dodge(width=0.5),stat="identity") +
geom_line(aes(Period,Estimate,group=1),colour="gray50",linetype="dashed",size=0.5) +
theme_bw() + 
ylim(0.10,0.18) +
labs(title="",x="",y="Misclassification probabilities") + 
theme(axis.text.y=element_text(face="bold",size=10)) + 
theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
guides(color=guide_legend(title="Sampling period")) + 
theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))
remove(phirst.es)

