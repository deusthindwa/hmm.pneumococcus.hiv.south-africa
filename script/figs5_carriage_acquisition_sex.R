#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#rerun model with sex included as covariate for carriage acquisition
phirst.fu <- arrange(phirst.fu,visit_id)
p.acqsex <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=rbind(c(0.0,0.9),c(0.9,0.0)), ematrix=rbind(c(1.0,0.0),c(0.1,0.9)),
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat+sxf+sxm),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=200000))

#show swabbing visit vs carriage prevalence by sex
phirst.d1 <- subset(subset(phirst.fu,select=c(visit_id,state,agecat,hiv,sex)),state !=9 & agecat=="Adult")
phirst.d1$visit_id <- as.integer(substr(phirst.d1$visit_id,10,12))
phirst.d2 <- phirst.d1 %>% group_by(sex,hiv,visit_id) %>% tally(state==2)
phirst.d1 <- phirst.d1 %>% group_by(sex,hiv,visit_id) %>% tally()
phirst.d1$nPos <- phirst.d2$n; phirst.d1$prev <- phirst.d1$nPos/phirst.d1$n; remove(phirst.d2)
phirst.d1$d1group <- if_else(phirst.d1$sex=="Male" & phirst.d1$hiv=="Neg","HIV- male adults",
                     if_else(phirst.d1$sex=="Male" & phirst.d1$hiv=="Pos","HIV+ male adults",
                     if_else(phirst.d1$sex=="Female" & phirst.d1$hiv=="Neg","HIV- female adults","HIV+ female adults")))

A<-ggplot(phirst.d1, aes(visit_id,prev*100,color=d1group)) + 
  geom_line(size=1, na.rm=TRUE) + 
  scale_color_manual(values=c("black","red","orange","blue")) +
  labs(title="A",x="Sampling visit",y="Carriage prevalence (%)") + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10)) +
  guides(color=guide_legend(title="HIV-sex status")) +
  theme(legend.position="right")

#within household acquisition probability in children by HIV status in female adults
p.modela <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",sxf="No"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",sxf="No"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",sxf="No"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",sxf="No"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-"),
                         "hh_hiv"=c("No","No","No","No"))

phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

p.modela <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",sxf="Yes"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",sxf="Yes"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",sxf="Yes"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",sxf="Yes"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-"),
                         "hh_hiv"=c("Yes","Yes","Yes","Yes"))

phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])
phirst.es <- rbind(phirst.es0,phirst.es1)

B <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv),size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.04,0.22) + 
  labs(title="B",x="",y="Household daily acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title="Age-HIV status"),shape=guide_legend(title="Households with\nHIV+ female adult(s)"),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))


#within household acquisition probability in children by HIV status in male adults
p.modela <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",sxm="No"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",sxm="No"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",sxm="No"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",sxm="No"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-"),
                         "hh_hiv"=c("No","No","No","No"))

phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

p.modela <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",sxm="Yes"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",sxm="Yes"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",sxm="Yes"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.acqsex,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",sxm="Yes"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-"),
                         "hh_hiv"=c("Yes","Yes","Yes","Yes"))

phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])
phirst.es <- rbind(phirst.es0,phirst.es1)

C <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv),size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.04,0.22) + 
  labs(title="C",x="",y="Household daily acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title="Age-HIV status"),shape=guide_legend(title="Households with\nHIV+ male adult(s)"),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

A + B + C
