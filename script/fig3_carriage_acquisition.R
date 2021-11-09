#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#hazard ratios for pneumococcal acquisition/clearance rates (Table 2)
hazard.msm(p.model4,hazard.scale=1,cl=0.95)

#community acquisition probability
p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="cmtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="cmtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="cmtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="cmtx"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)

phirst.es <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                        "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                        "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"))

phirst.es$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])

A<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid),shape=8,size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.01,0.08) +
  labs(title="A",x="",y="Community daily acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#format and present community acquisition probabilities (supplementary Table 2)
phirst.es <- rename(phirst.es, c("iid"="Table.1B","age"="Age","hiv"="HIV","carry.est"="Estimate","Lcarry.est"="Lower.bound","Ucarry.est"="Upper.bound"))

#within household acquisition probability
p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "hh_hiv"=c("No","No","No","No","No","No"))

phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],NA,p.modelf$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],NA,p.modelf$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],NA,p.modelf$U[1,2])

p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV-","Adult HIV+"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV-","HIV+"),
                         "hh_hiv"=c("Yes","Yes","Yes","Yes","Yes","Yes"))

phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
phirst.es <- rbind(phirst.es0,phirst.es1)

B<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv),size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.01,0.08) + 
  labs(title="B",x="",y="Household daily acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title="Age-HIV status"),shape=guide_legend(title="HH adult HIV+?"),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#format and present household acquisition probabilities (supplementary Table 2)
phirst.es <- rename(phirst.es, c("iid"="Table.1A","age"="Age","hiv"="HIV","hh_hiv"="Household.adult.HIV","carry.est"="Estimate","Lcarry.est"="Lower.bound","Ucarry.est"="Upper.bound"))

A + B


#estimate the number of episodes
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)

envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
envisits.msm(p.model4,fromt=0,tot=365.25,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)




