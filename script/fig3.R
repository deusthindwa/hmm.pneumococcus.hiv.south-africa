#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#household adult HIV-dependent acquisition rates
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

phirst.es <- rbind(phirst.es0,phirst.es1)
A<-ggplot(phirst.es) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  scale_y_continuous(labels=scales::number_format(accuracy=0.001)) +
  labs(title="A",x="",y="HH acquisition per day") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#household adult HIV-dependent probability rates
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

phirst.es <- rbind(phirst.es0,phirst.es1)
B<-ggplot(phirst.es) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  scale_y_continuous(labels=scales::number_format(accuracy=0.01)) +
  labs(title="B",x="",y="HH acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

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
  geom_point(aes(iid,carry.est,color=iid),shape=1,size=2,position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est),width=0.1,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  scale_y_continuous(labels=scales::number_format(accuracy=0.001)) +
  labs(title="C",x="",y="CM acquisition per day") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#community acquisition probability
p.modela <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Child",tx="cmtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Pos",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=289,covariates=list(hiv="Neg",agecat="Adult",tx="cmtx"),ci="normal",cl=0.95)

phirst.es <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"))
phirst.es$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

D <- ggplot(phirst.es) +
  geom_point(aes(iid,carry.est,color=iid),shape=1,size=2,position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est),width=0.1,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  scale_y_continuous(labels=scales::number_format(accuracy=0.01)) +
  labs(title="D",x="",y="CM acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

remove(phirst.es,phirst.es0,phirst.es1,p.modela,p.modelb,p.modelc,p.modeld)
print(ggarrange(A,B,C,D,ncol=2,nrow=2,common.legend=TRUE,legend="bottom",vjust=-2))
