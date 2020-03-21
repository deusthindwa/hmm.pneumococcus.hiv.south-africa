#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#average carriage duration from selected model
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",abx="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child",abx="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",abx="No"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-"),"age"=c("Child","Child","Adult"),"hiv"=c("HIV+","HIV-","HIV-"),"abxcat"=c("No antibiotics","No antibiotics","No antibiotics"))
phirst.es0$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1])
phirst.es0$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1])
phirst.es0$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1])  

p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Child",abx="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Child",abx="Yes"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("Child HIV+","Child HIV-","Adult HIV-"),"age"=c("Child","Child","Adult"),"hiv"=c("HIV+","HIV-","HIV-"),"abxcat"=c("Antibiotic use","Antibiotic use","Antibiotic use"))
phirst.es1$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1])
phirst.es1$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1])
phirst.es1$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1])
phirst.es1 <- rbind(phirst.es0,phirst.es1)

C <- ggplot(phirst.es1) +
  geom_point(aes(iid,1/clear.est,color=iid,shape=abxcat), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est,shape=abxcat),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="C",x="",y="Average carriage duration (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))
