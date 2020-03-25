#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#average carriage duration by ART status
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",artv="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child",artv="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",artv="No"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",artv="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"),"hiv"=c("HIV+","HIV-","HIV-","HIV+"),"abxcat"=c("No","No","No","No"))
phirst.es0$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
phirst.es0$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
phirst.es0$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  

p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",artv="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",artv="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Adult HIV+"),"age"=c("Child","Adult"),"hiv"=c("HIV+","HIV+"),"abxcat"=c("Yes","Yes"))
phirst.es1$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1])
phirst.es1$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1])
phirst.es1$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1])  

phirst.es <- rbind(phirst.es0,phirst.es1)
A <- ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est,shape=abxcat),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,1/clear.est,color=iid,shape=abxcat),size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() +
  ylim(0,100) +
  labs(title="A",x="",y="Average duration by ART (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title=""),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#average carriage duration by ABX status
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",abx="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child",abx="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",abx="No"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",abx="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"),"hiv"=c("HIV+","HIV-","HIV-","HIV+"),"abxcat"=c("No","No","No","No"))
phirst.es0$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
phirst.es0$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
phirst.es0$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  

p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",abx="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child",abx="Yes"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"),"hiv"=c("HIV+","HIV-","HIV-","HIV+"),"abxcat"=c("Yes","Yes","Yes","Yes"))
phirst.es1$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
phirst.es1$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
phirst.es1$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  

phirst.es <- rbind(phirst.es0,phirst.es1)
B <- ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est,shape=abxcat),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,1/clear.est,color=iid,shape=abxcat),size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0,100) +
  labs(title="B",x="",y="Average duration by ABX (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title=""),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#average carriage duration
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult"),ci="normal",cl=0.95)
phirst.es <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"),"hiv"=c("HIV+","HIV-","HIV-","HIV+"))
phirst.es$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
phirst.es$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
phirst.es$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  

C <- ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,1/clear.est,color=iid),shape=8,size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0,100) +
  labs(title="C",x="",y="Average duration overall (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=FALSE,fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#probability of clearance
p.modela <- pmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult"),ci="normal",cl=0.95)
phirst.es <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-","Adult HIV+"),"age"=c("Child","Child","Adult","Adult"),"hiv"=c("HIV+","HIV-","HIV-","HIV+"))
phirst.es$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
phirst.es$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
phirst.es$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  

D <- ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=Lclear.est,ymax=Uclear.est),width=0,size=0.6,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,clear.est,color=iid),shape=8,size=1.5,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0,0.08) +
  labs(title="D",x="",y="Probability of clearance") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=FALSE,fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

remove(phirst.es,phirst.es0,phirst.es1,p.modela,p.modelb,p.modelc,p.modeld)
print(ggarrange(A,B,C,D,ncol=4,common.legend=TRUE,legend="bottom"))
