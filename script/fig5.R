#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#---------------transition intensity estimates from selected model
set.seed(1988)
for(i in 0:1){
  p.modela <- qmatrix.msm(p.model5a, covariates=list(hiv=1,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelb <- qmatrix.msm(p.model5a, covariates=list(hiv=0,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelc <- qmatrix.msm(p.model5a, covariates=list(hiv=1,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modeld <- qmatrix.msm(p.model5a, covariates=list(hiv=0,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  if(i==0){
    phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
    phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
    phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])
  }
  else{
    phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
    phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
    phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])
  }
  if(i==0){
    phirst.es2 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es2$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
    phirst.es2$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
    phirst.es2$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  
  }
  else{
    phirst.es3 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es3$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
    phirst.es3$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
    phirst.es3$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])
  }
}

phirst.es0 <- rbind(phirst.es0,phirst.es1)
phirst.es1 <- rbind(phirst.es2,phirst.es3)

A <- ggplot(phirst.es0) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="A",x="",y="HH carriage acquistion per day") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

C <- ggplot(phirst.es1) +
  geom_point(aes(iid,1/clear.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="C",x="",y="Average carriage duration (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#---------------transition probability estimates from selected model
set.seed(1988)
for(i in 0:1){
  p.modela <- pmatrix.msm(p.model5a, covariates=list(hiv=1,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelb <- pmatrix.msm(p.model5a, covariates=list(hiv=0,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelc <- pmatrix.msm(p.model5a, covariates=list(hiv=1,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modeld <- pmatrix.msm(p.model5a, covariates=list(hiv=0,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  if(i==0){
    phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es0$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
    phirst.es0$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
    phirst.es0$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])
  }
  else{
    phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es1$carry.est <- c(p.modela$estimates[1,2],p.modelb$estimates[1,2],p.modelc$estimates[1,2],p.modeld$estimates[1,2])
    phirst.es1$Lcarry.est <- c(p.modela$L[1,2],p.modelb$L[1,2],p.modelc$L[1,2],p.modeld$L[1,2])
    phirst.es1$Ucarry.est <- c(p.modela$U[1,2],p.modelb$U[1,2],p.modelc$U[1,2],p.modeld$U[1,2])
  }
  if(i==0){
    phirst.es2 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es2$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
    phirst.es2$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
    phirst.es2$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])  
  }
  else{
    phirst.es3 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es3$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1])
    phirst.es3$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1])
    phirst.es3$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1])
  }
}
phirst.es0 <- rbind(phirst.es0,phirst.es1)
phirst.es1 <- rbind(phirst.es2,phirst.es3)

B <- ggplot(phirst.es0) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="B",x="",y="HH carriage acquistion probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

D <- ggplot(phirst.es1) +
  geom_point(aes(iid,clear.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lclear.est,ymax=Uclear.est,shape=hh_hiv),width=0.2,size=0.6,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="D",x="",y="Carriage clearance probability") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

remove(phirst.es0,phirst.es1,phirst.es2,phirst.es3)
print(ggarrange(A,B,C,D,ncol=2,nrow=2,common.legend=TRUE,legend="right",vjust=-2))
