#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#hazard ratios for pneumococcal acquisition/clearance rates
hazard.msm(p.model4,hazard.scale=1,cl=0.95)

#household adult HIV-dependent acquisition rates
p.modela <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-"),"age"=c("Child","Child","Adult"), "hiv"=c("HIV+","HIV-","HIV-"),"hh_hiv"=c("No","No","No"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2])

p.modela <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"),"hh_hiv"=c("Yes","Yes","Yes","Yes"))
phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

phirst.es <- rbind(phirst.es0,phirst.es1)
A<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,4,2,3)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv),size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.01,0.085) +
  labs(title="A",x="",y="HH daily acquisition rate") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="HH adult HIV+?"),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#household adult HIV-dependent probability rates
p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-"),"age"=c("Child","Child","Adult"), "hiv"=c("HIV+","HIV-","HIV-"),"hh_hiv"=c("No","No","No"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2])

p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"),"hh_hiv"=c("Yes","Yes","Yes","Yes"))
phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

phirst.es <- rbind(phirst.es0,phirst.es1)
B<-ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,4,2,3)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv),size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.01,0.085) +
  labs(title="A",x="",y="Household daily acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="HH adult HIV+?"),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#community acquisition rates
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)

phirst.es <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"))
phirst.es$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

C <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid),shape=8,size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.01,0.085) +
  labs(title="C",x="",y="CM daily cquisition rate") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=FALSE,fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#community acquisition probability
p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)

phirst.es <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"))
phirst.es$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

D <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid),shape=8,size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.01,0.085) +
  labs(title="B",x="",y="Community daily acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

remove(phirst.es,phirst.es0,phirst.es1,p.modela,p.modelb,p.modelc,p.modeld)
print(ggarrange(B,D,ncol=2,nrow=1,common.legend=TRUE,legend="bottom"))
