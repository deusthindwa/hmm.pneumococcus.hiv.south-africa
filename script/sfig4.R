#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#refit the hidden Markov model with # of adult HIV+ in the household as covariate
phirst.hz <- subset(phirst.fu,select=c(ind_id,dys,state,obst,agecat,hiv,tx,hhsize))
phirst.hz$hhsize <- as.factor(if_else(phirst.hz$hhsize<=5,"<6",if_else(phirst.hz$hhsize>=6 & phirst.hz$hhsize<=10,"6-10","11+")))

#initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.9),c(0.9,0.0))
rownames(matrix.Q) <- c("Clear","Carry")
colnames(matrix.Q) <- c("Clear","Carry")

#initiate emission matrix E
matrix.E <- rbind(c(1.0,0.0),c(0.1,0.9))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carry")

statetable.msm(state,ind_id,data=phirst.hz)
p.hhsize <- msm(state~dys,subject=ind_id,data=phirst.hz,
               qmatrix=matrix.Q, ematrix=matrix.E,
               covariates=list("1-2"=~agecat+hiv+hhsize+tx),
               censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
               opt.method="bobyqa", control=list(maxfun=100000))

#hazard ratios for pneumococcal acquisition
hazard.msm(p.hhsize,hazard.scale=1,cl=0.95)

#compute the acquisition rates by household size
hazard.msm(p.hhsize,hazard.scale=1,cl=0.95)
p.modela <- qmatrix.msm(p.hhsize,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.hhsize,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.hhsize,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.hhsize,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("<6","<6","<6","<6"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

p.modele <- qmatrix.msm(p.hhsize,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.hhsize,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
p.modelg <- qmatrix.msm(p.hhsize,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
p.modelh <- qmatrix.msm(p.hhsize,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("6-10","6-10","6-10","6-10"))
phirst.es1$carry.est <- c(p.modele$estimates[1,2],p.modelf$estimates[1,2],p.modelg$estimates[1,2],p.modelh$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modele$L[1,2],p.modelf$L[1,2],p.modelg$L[1,2],p.modelh$L[1,2])
phirst.es1$Ucarry.est <- c(p.modele$U[1,2],p.modelf$U[1,2],p.modelg$U[1,2],p.modelh$U[1,2])

p.modelj <- qmatrix.msm(p.hhsize,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",hhsize="11+"),ci="normal",cl=0.95)
p.modelk <- qmatrix.msm(p.hhsize,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsize="11+"),ci="normal",cl=0.95)
p.modell <- qmatrix.msm(p.hhsize,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsize="11+"),ci="normal",cl=0.95)
phirst.es2 <- data.frame("iid"=c("HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Adult","Adult"), "hiv"=c("HIV-","HIV+","HIV-"),"status"=c("11+","11+","11+"))
phirst.es2$carry.est <- c(p.modelj$estimates[1,2],p.modelk$estimates[1,2],p.modell$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modelj$L[1,2],p.modelk$L[1,2],p.modell$L[1,2])
phirst.es2$Ucarry.est <- c(p.modelj$U[1,2],p.modelk$U[1,2],p.modell$U[1,2])

phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2)
remove(phirst.es0,phirst.es1,phirst.es2,p.modela,p.modelb,p.modelc,p.modeld,p.modele,p.modelf,p.modelg,p.modelh,p.modelj,p.modelk,p.modell)

A <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.025,0.22) +
  labs(title="A",x="",y="HH daily acquisition rate") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Household size")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#compute the acquisition probabilities by household size
p.modela <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsize="<6"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("<6","<6","<6","<6"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

p.modele <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
p.modelg <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
p.modelh <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsize="6-10"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("6-10","6-10","6-10","6-10"))
phirst.es1$carry.est <- c(p.modele$estimates[1,2],p.modelf$estimates[1,2],p.modelg$estimates[1,2],p.modelh$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modele$L[1,2],p.modelf$L[1,2],p.modelg$L[1,2],p.modelh$L[1,2])
phirst.es1$Ucarry.est <- c(p.modele$U[1,2],p.modelf$U[1,2],p.modelg$U[1,2],p.modelh$U[1,2])

p.modelj <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",hhsize="11+"),ci="normal",cl=0.95)
p.modelk <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",hhsize="11+"),ci="normal",cl=0.95)
p.modell <- pmatrix.msm(p.hhsize,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",hhsize="11+"),ci="normal",cl=0.95)
phirst.es2 <- data.frame("iid"=c("HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Adult","Adult"), "hiv"=c("HIV-","HIV+","HIV-"),"status"=c("11+","11+","11+"))
phirst.es2$carry.est <- c(p.modelj$estimates[1,2],p.modelk$estimates[1,2],p.modell$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modelj$L[1,2],p.modelk$L[1,2],p.modell$L[1,2])
phirst.es2$Ucarry.est <- c(p.modelj$U[1,2],p.modelk$U[1,2],p.modell$U[1,2])

phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2)
remove(phirst.es0,phirst.es1,phirst.es2,p.modela,p.modelb,p.modelc,p.modeld,p.modele,p.modelf,p.modelg,p.modelh,p.modeli,p.modelj,p.modelk,p.modell)

B <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.025,0.22) +
  labs(title="B",x="",y="HH daily acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Household size")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

print(ggarrange(A,B,ncol=2,common.legend=TRUE,legend="right"))


#compute the clearance rates by reported ART status other than viral load-based ART status
p.modelART <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q, ematrix=matrix.E,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx,"2-1"=~agecat+hiv+abx+artr),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

hazard.msm(p.modelART,hazard.scale=1,cl=0.95)
