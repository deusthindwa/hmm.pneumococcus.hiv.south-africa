#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#---------------forecasted total length of time spent in each trasient state
set.seed(1988)
for(i in 0:1){
  p.modela <- totlos.msm(p.model5,fromt=0,tot=289,covariates=list(hiv=1,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelb <- totlos.msm(p.model5,fromt=0,tot=289,covariates=list(hiv=0,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelc <- totlos.msm(p.model5,fromt=0,tot=289,covariates=list(hiv=1,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modeld <- totlos.msm(p.model5,fromt=0,tot=289,covariates=list(hiv=0,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  if(i==0){
    phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es0$carry.est <- c(p.modela[1,2],p.modelb[1,2],p.modelc[1,2],p.modeld[1,2])
    phirst.es0$Lcarry.est <- c(p.modela[2,2],p.modelb[2,2],p.modelc[2,2],p.modeld[2,2])
    phirst.es0$Ucarry.est <- c(p.modela[3,2],p.modelb[3,2],p.modelc[3,2],p.modeld[3,2])
  }
  else{
    phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es1$carry.est <- c(p.modela[1,2],p.modelb[1,2],p.modelc[1,2],p.modeld[1,2])
    phirst.es1$Lcarry.est <- c(p.modela[2,2],p.modelb[2,2],p.modelc[2,2],p.modeld[2,2])
    phirst.es1$Ucarry.est <- c(p.modela[3,2],p.modelb[3,2],p.modelc[3,2],p.modeld[3,2])
  }
} 
phirst.es <- rbind(phirst.es0,phirst.es1)
A<-ggplot(phirst.es) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="A",x="",y="Total carriage duration (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#---------------expected time until Markov process first enters a carrying state (hitting time)
set.seed(1988)
for(i in 0:1){
  p.modela <- efpt.msm(p.model5,tostate=2,covariates=list(hiv=1,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelb <- efpt.msm(p.model5,tostate=2,covariates=list(hiv=0,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelc <- efpt.msm(p.model5,tostate=2,covariates=list(hiv=1,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modeld <- efpt.msm(p.model5,tostate=2,covariates=list(hiv=0,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  if(i==0){
    phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es0$carry.est <- c(p.modela[1,1],p.modelb[1,1],p.modelc[1,1],p.modeld[1,1])
    phirst.es0$Lcarry.est <- c(p.modela[2,1],p.modelb[2,1],p.modelc[2,1],p.modeld[2,1])
    phirst.es0$Ucarry.est <- c(p.modela[3,1],p.modelb[3,1],p.modelc[3,1],p.modeld[3,1])
  }
  else{
    phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es1$carry.est <- c(p.modela[1,1],p.modelb[1,1],p.modelc[1,1],p.modeld[1,1])
    phirst.es1$Lcarry.est <- c(p.modela[2,1],p.modelb[2,1],p.modelc[2,1],p.modeld[2,1])
    phirst.es1$Ucarry.est <- c(p.modela[3,1],p.modelb[3,1],p.modelc[3,1],p.modeld[3,1]) 
  }
} 
phirst.es <- rbind(phirst.es0,phirst.es1)
B<-ggplot(phirst.es) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="B",x="",y="Hitting time (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

#---------------expected number of visits to carriage state
set.seed(1988)
for(i in 0:1){
  p.modela <- envisits.msm(p.model5, fromt=0,tot=289,covariates=list(hiv=1,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelb <- envisits.msm(p.model5, fromt=0,tot=289,covariates=list(hiv=0,age=0,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modelc <- envisits.msm(p.model5, fromt=0,tot=289,covariates=list(hiv=1,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  p.modeld <- envisits.msm(p.model5, fromt=0,tot=289,covariates=list(hiv=0,age=1,apncc=1,ahivc=i),ci="normal",cl=0.95)
  if(i==0){
    phirst.es0 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV- HH","Adult HIV- HH","Adult HIV- HH","Adult HIV- HH"))
    phirst.es0$carry.est <- c(p.modela[1,2],p.modelb[1,2],p.modelc[1,2],p.modeld[1,2])
    phirst.es0$Lcarry.est <- c(p.modela[2,2],p.modelb[2,2],p.modelc[2,2],p.modeld[2,2])
    phirst.es0$Ucarry.est <- c(p.modela[3,2],p.modelb[3,2],p.modelc[3,2],p.modeld[3,2])
  }
  else{
    phirst.es1 <- data.frame("iid"=c("HIV+ Child","HIV- Child","HIV+ Adult","HIV- Adult"),"age"=c("Child","Child","Adult","Adult"), "hiv"=c("HIV+","HIV-","HIV+","HIV-"),"hh_hiv"=c("Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH","Adult HIV+ HH"))
    phirst.es1$carry.est <- c(p.modela[1,2],p.modelb[1,2],p.modelc[1,2],p.modeld[1,2])
    phirst.es1$Lcarry.est <- c(p.modela[2,2],p.modelb[2,2],p.modelc[2,2],p.modeld[2,2])
    phirst.es1$Ucarry.est <- c(p.modela[3,2],p.modelb[3,2],p.modelc[3,2],p.modeld[3,2])
  }
}
phirst.es <- rbind(phirst.es0,phirst.es1)
C<-ggplot(phirst.es) +
  geom_point(aes(iid,carry.est,color=iid,shape=hh_hiv), size=1.5, position=position_dodge(width=0.5),stat="identity") +
  geom_errorbar(aes(iid,color=iid,ymin=Lcarry.est,ymax=Ucarry.est,shape=hh_hiv),width=0.1,size=0.5,position=position_dodge(width=0.5)) +
  theme_bw() + 
  labs(title="C",x="",y="Average number of carriage acquisitions") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) +
  theme(legend.text=element_text(size=10),legend.position="none",legend.title=element_text(face="bold",size=10))

remove(phirst.es,phirst.es0,phirst.es1)
print(ggarrange(A,B,C,ncol=3,nrow=1,common.legend=TRUE,legend="right",hjust=-4))
