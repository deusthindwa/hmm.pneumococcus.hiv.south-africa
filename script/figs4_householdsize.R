#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#age distribution of the population
A <- ggplot(phirst.fu,aes(x=age,color=hiv)) +
  geom_histogram(fill="white",binwidth=1,position=position_dodge(width=1)) +
  theme_bw() + 
  labs(title="A",x="Age (years)",y="Frequency") + 
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  theme(axis.text.y=element_text(face="bold",size=10),axis.text.x=element_text(face="bold",size=10)) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Household size")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#rerun model only for acquisition
p.hhsize <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

#compute the acquisition probabilities by household size
p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",hhsizecat="<6"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",hhsizecat="<6"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",hhsizecat="<6"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",hhsizecat="<6"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsizecat="<6"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsizecat="<6"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "status"=c("<6","<6","<6","<6","<6","<6"))

phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------

p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",hhsizecat="6-10"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",hhsizecat="6-10"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",hhsizecat="6-10"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",hhsizecat="6-10"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsizecat="6-10"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsizecat="6-10"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "status"=c("6-10","6-10","6-10","6-10","6-10","6-10"))

phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------

p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",hhsizecat="11+"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",hhsizecat="11+"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",hhsizecat="11+"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",hhsizecat="11+"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsizecat="11+"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsizecat="11+"),ci="normal",cl=0.95)

phirst.es2 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "status"=c("11+","11+","11+","11+","11+","11+"))

phirst.es2$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es2$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------

phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2)

B <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status),size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="B",x="",y="Household daily acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Household size")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

A + B
