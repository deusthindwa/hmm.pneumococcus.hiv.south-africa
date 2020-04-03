#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 11/3/2020

#forecasted total length of time spent in each trasient state
  p.modela <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Child"),ci="normal",cl=0.95)
  p.modelb <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Child"),ci="normal",cl=0.95)
  p.modelc <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Adult"),ci="normal",cl=0.95)
  p.modeld <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Adult"),ci="normal",cl=0.95)

  p.modele <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Child",artv="Yes"),ci="normal",cl=0.95)
  p.modelf <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Child",artv="Yes"),ci="normal",cl=0.95)
  p.modelg <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Adult",artv="Yes"),ci="normal",cl=0.95)
  p.modelh <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Adult",artv="Yes"),ci="normal",cl=0.95)
  
  p.modeli <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Child",abx="Yes"),ci="normal",cl=0.95)
  p.modelj <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Child",abx="Yes"),ci="normal",cl=0.95)
  p.modelk <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)
  p.modell <- totlos.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)
  
  phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("Overall","Overall","Overall","Overall"))
  phirst.es0$carry.est <- c(p.modela[1,2],p.modelb[1,2],p.modelc[1,2],p.modeld[1,2])
  phirst.es0$Lcarry.est <- c(p.modela[2,2],p.modelb[2,2],p.modelc[2,2],p.modeld[2,2])
  phirst.es0$Ucarry.est <- c(p.modela[3,2],p.modelb[3,2],p.modelc[3,2],p.modeld[3,2])

  phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("On ART","On ART","On ART","On ART"))
  phirst.es1$carry.est <- c(p.modele[1,2],p.modelf[1,2],p.modelg[1,2],p.modelh[1,2])
  phirst.es1$Lcarry.est <- c(p.modele[2,2],p.modelf[2,2],p.modelg[2,2],p.modelh[2,2])
  phirst.es1$Ucarry.est <- c(p.modele[3,2],p.modelf[3,2],p.modelg[3,2],p.modelh[3,2])
  
  phirst.es2 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("On Abx","On Abx","On Abx","On Abx"))
  phirst.es2$carry.est <- c(p.modeli[1,2],p.modelj[1,2],p.modelk[1,2],p.modell[1,2])
  phirst.es2$Lcarry.est <- c(p.modeli[2,2],p.modelj[2,2],p.modelk[2,2],p.modell[2,2])
  phirst.es2$Ucarry.est <- c(p.modeli[3,2],p.modelj[3,2],p.modelk[3,2],p.modell[3,2])

phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2)
A<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="A",x="",y="Total carriage duration (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Status")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#expected time until Markov process first enters a carrying state (hitting time)
  p.modela <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  p.modelb <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  p.modeld <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  p.modelc <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)

  phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV- Adult","HIV+ Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"),"hh_hiv"=c("Yes","Yes","Yes","Yes"))
  phirst.es0$carry.est <- c(p.modela[1,1],p.modelb[1,1],p.modelc[1,1],p.modeld[1,1])
  phirst.es0$Lcarry.est <- c(p.modela[2,1],p.modelb[2,1],p.modelc[2,1],p.modeld[2,1])
  phirst.es0$Ucarry.est <- c(p.modela[3,1],p.modelb[3,1],p.modelc[3,1],p.modeld[3,1]) 
  
  p.modele <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
  p.modelf <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
  p.modelg <- efpt.msm(p.model4,tostate=2,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
  
  phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV- Adult"),"age"=c("Child","Child","Adult"), "hiv"=c("HIV+","HIV-","HIV-"),"hh_hiv"=c("No","No","No"))
  phirst.es1$carry.est <- c(p.modele[1,1],p.modelf[1,1],p.modelg[1,1])
  phirst.es1$Lcarry.est <- c(p.modele[2,1],p.modelf[2,1],p.modelg[2,1])
  phirst.es1$Ucarry.est <- c(p.modele[3,1],p.modelf[3,1],p.modelg[3,1])

phirst.es <- rbind(phirst.es0,phirst.es1)
B<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv),size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="B",x="",y="Hitting time (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="HH adult HIV+?")) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#expected number of visits to carriage state
  p.modela <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  p.modelb <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  p.modelc <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  p.modeld <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
  
  p.modele <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
  p.modelf <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
  p.modelg <- envisits.msm(p.model4,fromt=0,tot=289,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
  
  phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV- Adult","HIV+ Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"),"hh_hiv"=c("Yes","Yes","Yes","Yes"))
  phirst.es0$carry.est <- c(p.modela[1,2],p.modelb[1,2],p.modelc[1,2],p.modeld[1,2])
  phirst.es0$Lcarry.est <- c(p.modela[2,2],p.modelb[2,2],p.modelc[2,2],p.modeld[2,2])
  phirst.es0$Ucarry.est <- c(p.modela[3,2],p.modelb[3,2],p.modelc[3,2],p.modeld[3,2]) 
  
  phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV- Adult"),"age"=c("Child","Child","Adult"), "hiv"=c("HIV+","HIV-","HIV-"),"hh_hiv"=c("No","No","No"))
  phirst.es1$carry.est <- c(p.modele[1,2],p.modelf[1,2],p.modelg[1,2])
  phirst.es1$Lcarry.est <- c(p.modele[2,2],p.modelf[2,2],p.modelg[2,2])
  phirst.es1$Ucarry.est <- c(p.modele[3,2],p.modelf[3,2],p.modelg[3,2])
  
phirst.es <- rbind(phirst.es0,phirst.es1)
C<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="C",x="",y="Average number of acquisitions") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="HH adult HIV+?")) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

remove(phirst.es,phirst.es0,phirst.es1)
print(ggarrange(A,B,C,ncol=3,nrow=1,common.legend=FALSE,legend="right",hjust=-4))
