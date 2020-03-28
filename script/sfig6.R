#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#initiate transition intensity matrix Q
matrix.Q <- rbind(c(0.0,0.9),c(0.9,0.0))
rownames(matrix.Q) <- c("Clear","Carry")
colnames(matrix.Q) <- c("Clear","Carry")

#initiate emission matrix E
matrix.E <- rbind(c(1.0,0.0),c(0.1,0.9))
colnames(matrix.E) <- c("SwabNeg","SwabPos")
rownames(matrix.E) <- c("Clear","Carry")

statetable.msm(state,ind_id,data=phirst.fu)
p.misclass <- msm(state~dys,subject=ind_id,data=phirst.fu,
              qmatrix=matrix.Q, ematrix=matrix.E,
              covariates=list("1-2"=~agecat+hiv),
              misccovariates=~agecat,
              censor=9, censor.states=c(1,2), obstrue=obst,est.initprobs=T,
              opt.method="bobyqa",control=list(maxfun=100000))

#compute the acquisition rates by # of adult HIV+
p.modela <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("1","1","1","1"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

p.modele <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modelg <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modelh <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("2","2","2","2"))
phirst.es1$carry.est <- c(p.modele$estimates[1,2],p.modelf$estimates[1,2],p.modelg$estimates[1,2],p.modelh$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modele$L[1,2],p.modelf$L[1,2],p.modelg$L[1,2],p.modelh$L[1,2])
phirst.es1$Ucarry.est <- c(p.modele$U[1,2],p.modelf$U[1,2],p.modelg$U[1,2],p.modelh$U[1,2])

p.modeli <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modelj <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modelk <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modell <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
phirst.es2 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("3","3","3","3"))
phirst.es2$carry.est <- c(p.modeli$estimates[1,2],p.modelj$estimates[1,2],p.modelk$estimates[1,2],p.modell$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modeli$L[1,2],p.modelj$L[1,2],p.modelk$L[1,2],p.modell$L[1,2])
phirst.es2$Ucarry.est <- c(p.modeli$U[1,2],p.modelj$U[1,2],p.modelk$U[1,2],p.modell$U[1,2])

p.modelm <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modeln <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modelo <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modelp <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
phirst.es3 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("4","4","4","4"))
phirst.es3$carry.est <- c(p.modelm$estimates[1,2],p.modeln$estimates[1,2],p.modelo$estimates[1,2],p.modelp$estimates[1,2])
phirst.es3$Lcarry.est <- c(p.modelm$L[1,2],p.modeln$L[1,2],p.modelo$L[1,2],p.modelp$L[1,2])
phirst.es3$Ucarry.est <- c(p.modelm$U[1,2],p.modeln$U[1,2],p.modelo$U[1,2],p.modelp$U[1,2])

p.modelq <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Child",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modelr <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Child",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.models <- qmatrix.msm(p.snstvty,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modelt <- qmatrix.msm(p.snstvty,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
phirst.es4 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("5","5","5","5"))
phirst.es4$carry.est <- c(p.modelq$estimates[1,2],p.modelr$estimates[1,2],p.models$estimates[1,2],p.modelt$estimates[1,2])
phirst.es4$Lcarry.est <- c(p.modelq$L[1,2],p.modelr$L[1,2],p.models$L[1,2],p.modelt$L[1,2])
phirst.es4$Ucarry.est <- c(p.modelq$U[1,2],p.modelr$U[1,2],p.models$U[1,2],p.modelt$U[1,2])

phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2,phirst.es3,phirst.es4)
remove(phirst.es0,phirst.es1,phirst.es2,phirst.es3,phirst.es4)
remove(p.modela,p.modelb,p.modelc,p.modeld,p.modele,p.modelf,p.modelg,p.modelh,p.modeli,p.modelj,p.modelk,p.modell,p.modelm,p.modeln,p.modelo,p.modelp,p.modelq,p.modelr,p.models,p.modelt)

A <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="A",x="",y="HH acquisition per day") + 
  ylim(0,0.20) +
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="# of HH adult HIV+")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#compute the probabilities by specific sampling days
p.modela <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Pos",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Neg",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("0-99","0-99","0-99","0-99"))
phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])

p.modele <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Pos",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Neg",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelg <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelh <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("100-139","100-139","100-139","100-139"))
phirst.es1$carry.est <- c(p.modele$estimates[1,2],p.modelf$estimates[1,2],p.modelg$estimates[1,2],p.modelh$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modele$L[1,2],p.modelf$L[1,2],p.modelg$L[1,2],p.modelh$L[1,2])
phirst.es1$Ucarry.est <- c(p.modele$U[1,2],p.modelf$U[1,2],p.modelg$U[1,2],p.modelh$U[1,2])

p.modeli <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Pos",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelj <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Neg",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelk <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modell <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
phirst.es2 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("140-199","140-199","140-199","140-199"))
phirst.es2$carry.est <- c(p.modeli$estimates[1,2],p.modelj$estimates[1,2],p.modelk$estimates[1,2],p.modell$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modeli$L[1,2],p.modelj$L[1,2],p.modelk$L[1,2],p.modell$L[1,2])
phirst.es2$Ucarry.est <- c(p.modeli$U[1,2],p.modelj$U[1,2],p.modelk$U[1,2],p.modell$U[1,2])

p.modelm <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Pos",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modeln <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Neg",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelo <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelp <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
phirst.es3 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("200-219","200-219","200-219","200-219"))
phirst.es3$carry.est <- c(p.modelm$estimates[1,2],p.modeln$estimates[1,2],p.modelo$estimates[1,2],p.modelp$estimates[1,2])
phirst.es3$Lcarry.est <- c(p.modelm$L[1,2],p.modeln$L[1,2],p.modelo$L[1,2],p.modelp$L[1,2])
phirst.es3$Ucarry.est <- c(p.modelm$U[1,2],p.modeln$U[1,2],p.modelo$U[1,2],p.modelp$U[1,2])

p.modelq <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Pos",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.modelr <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Neg",agecat="Child",tx="hhtx"),ci="normal",cl=0.95)
p.models <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelt <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
phirst.es4 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"status"=c("220-289","220-289","220-289","220-289"))
phirst.es4$carry.est <- c(p.modelq$estimates[1,2],p.modelr$estimates[1,2],p.models$estimates[1,2],p.modelt$estimates[1,2])
phirst.es4$Lcarry.est <- c(p.modelq$L[1,2],p.modelr$L[1,2],p.models$L[1,2],p.modelt$L[1,2])
phirst.es4$Ucarry.est <- c(p.modelq$U[1,2],p.modelr$U[1,2],p.models$U[1,2],p.modelt$U[1,2])

phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2,phirst.es3,phirst.es4)
remove(phirst.es0,phirst.es1,phirst.es2,phirst.es3,phirst.es4)
remove(p.modela,p.modelb,p.modelc,p.modeld,p.modele,p.modelf,p.modelg,p.modelh,p.modeli,p.modelj,p.modelk,p.modell,p.modelm,p.modeln,p.modelo,p.modelp,p.modelq,p.modelr,p.models,p.modelt)

B <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=factor(iid,levels(factor(iid))[c(1,3,2,4)]),ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  ylim(0.18,0.80) +
  labs(title="B",x="",y="HH acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Sampling days")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

print(ggarrange(A,B,ncol=2,common.legend=FALSE,legend="right"))
