#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#household HIV-dependent acquisition rates
p.modela <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-"),"age"=c("Child","Child","Adult"), "hiv"=c("HIV+","HIV-","HIV-"),"hh_hiv"=c("HH without Adult HIV+","HH without Adult HIV+","HH without Adult HIV+"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2])

p.modela <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"),"hh_hiv"=c("HH with Adult HIV+","HH with Adult HIV+","HH with Adult HIV+","HH with Adult HIV+"))
phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

phirst.es0 <- rbind(phirst.es0,phirst.es1)
A<-ggplot(phirst.es0) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  ylim(0.035,0.085) +
  labs(title="A",x="",y="Household acquistion per day") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#household HIV-dependent acquisition probability 
p.modela <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-"),"age"=c("Child","Child","Adult"), "hiv"=c("HIV+","HIV-","HIV-"),"hh_hiv"=c("HH without Adult HIV+","HH without Adult HIV+","HH without Adult HIV+"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2])

p.modela <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahivcat="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"),"hh_hiv"=c("HH with Adult HIV+","HH with Adult HIV+","HH with Adult HIV+","HH with Adult HIV+"))
phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

phirst.es0 <- rbind(phirst.es0,phirst.es1)
B<-ggplot(phirst.es0) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  ylim(0.35,0.85) +
  labs(title="B",x="",y="Household acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#community acquisition rates
p.modela <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4,covariates=list(hiv="Neg",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4,covariates=list(hiv="Pos",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV-","HIV+"))
phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

C<-ggplot(phirst.es1) +
  geom_point(aes(iid,carry.est,color=iid), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="C",x="",y="Community acquistion per day") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#community acquisition probability 
p.modela <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"),"hiv"=c("HIV+","HIV-","HIV-","HIV+"))
phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

D<-ggplot(phirst.es1) +
  geom_point(aes(iid,carry.est,color=iid), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="D",x="",y="Community acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

print(ggarrange(A,B,C,D,ncol=2,nrow=2,common.legend=TRUE,legend="right",vjust=-2))
remove(phirst.es0,phirst.es1,p.modela,p.modelb,p.modelc,p.modeld)
