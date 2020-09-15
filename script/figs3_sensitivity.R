#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#refit the hidden Markov model with varying # of adults HIV+ in the household as covariate
phirst.hs <- subset(phirst.fu,select=c(ind_id,dys,state,obst,agecat,hiv,tx,ahiv,hhsizecat))
phirst.hs$ahiv <- as.factor(if_else(phirst.hs$ahiv==0,"Zero",if_else(phirst.hs$ahiv==1,"One",if_else(phirst.hs$ahiv==2,"Two",if_else(phirst.hs$ahiv==3,"Three",if_else(phirst.hs$ahiv==4,"Four","Five"))))))

statetable.msm(state,ind_id,data=phirst.hs)

p.snstvty <- msm(state~dys,subject=ind_id,data=phirst.hs,
                qmatrix=rbind(c(0.0,0.9),c(0.9,0.0)), ematrix=rbind(c(1.0,0.0),c(0.1,0.9)),
                covariates=list("1-2"=~agecat+hiv+ahiv+tx+hhsizecat),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

#compute the acquisition rates by # of adult HIV+
p.modela <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahiv="Zero"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahiv="Zero"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahiv="Zero"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahiv="Zero"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Zero"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Zero"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("0","0","0","0","0","0"))

phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],NA,p.modelf$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],NA,p.modelf$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],NA,p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="One"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("1","1","1","1","1","1"))

phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Two"),ci="normal",cl=0.95)

phirst.es2 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("2","2","2","2","2","2"))

phirst.es2$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es2$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Three"),ci="normal",cl=0.95)

phirst.es3 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("3","3","3","3","3","3"))

phirst.es3$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es3$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es3$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Four"),ci="normal",cl=0.95)

phirst.es4 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("4","4","4","4","4","4"))

phirst.es4$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es4$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es4$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t=1,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx",ahiv="Five"),ci="normal",cl=0.95)

phirst.es5 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("5","5","5","5","5","5"))

phirst.es5$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es5$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es5$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2,phirst.es3,phirst.es4,phirst.es5)

A <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="A",x="",y="Household daily acquisition probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Number of HIV+ adults\nin the household")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#compute the acquisition probabilities by specific sampling days
p.modela <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t1=0,t=99,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("0-99","0-99","0-99","0-99","0-99","0-99"))

phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t1=100,t=39,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("100-139","100-139","100-139","100-139","100-139","100-139"))

phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t1=140,t=59,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)

phirst.es2 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("140-199","140-199","140-199","140-199","140-199","140-199"))

phirst.es2$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es2$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es2$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t1=200,t=20,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)

phirst.es3 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV+"),"status"=c("200-219","200-219","200-219","200-219","200-219","200-219"))

phirst.es3$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es3$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es3$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])

#--------------------------------------------
p.modela <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Pos",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Neg",agecat="Younger child",tx="hhtx"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Pos",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Neg",agecat="Older child",tx="hhtx"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Pos",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.snstvty,t1=220,t=69,covariates=list(hiv="Neg",agecat="Adult",tx="hhtx"),ci="normal",cl=0.95)

phirst.es4 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"), 
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),"status"=c("220-289","220-289","220-289","220-289","220-289","220-289"))

phirst.es4$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2],p.modele$estimates[1,2],p.modelf$estimates[1,2])
phirst.es4$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2],p.modele$L[1,2],p.modelf$L[1,2])
phirst.es4$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2],p.modele$U[1,2],p.modelf$U[1,2])
#--------------------------------------------
phirst.es <- rbind(phirst.es0,phirst.es1,phirst.es2,phirst.es3,phirst.es4)

B <- ggplot(phirst.es) +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=status),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,carry.est,color=iid,shape=status), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="B",x="",y="Household acquisition probability\nover the sampling period") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="Sampling period, days")) + 
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

A + B

