#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#average carriage duration by ABX status
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Younger child",abx="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Younger child",abx="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Older child",abx="No"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Older child",abx="No"),ci="normal",cl=0.95)
p.modele <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",abx="No"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",abx="No"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"),
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "abx"=c("No","No","No","No","No","No"))

phirst.es0$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1],p.modele$estimates[2,1],p.modelf$estimates[2,1])
phirst.es0$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1],p.modele$L[2,1],p.modelf$L[2,1])
phirst.es0$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1],p.modele$U[2,1],p.modelf$U[2,1])  

p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Younger child",abx="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Younger child",abx="Yes"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Older child",abx="Yes"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Older child",abx="Yes"),ci="normal",cl=0.95)
p.modele <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",abx="Yes"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"),
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "abx"=c("Yes","Yes","Yes","Yes","Yes","Yes"))

phirst.es1$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1],p.modele$estimates[2,1],p.modelf$estimates[2,1])
phirst.es1$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1],p.modele$L[2,1],p.modelf$L[2,1])
phirst.es1$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1],p.modele$U[2,1],p.modelf$U[2,1])  
phirst.es <- rbind(phirst.es0,phirst.es1)

A<-ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est,shape=abx),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,1/clear.est,color=iid,shape=abx),size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="A",x="",y="Average carriage duration by ABX (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  guides(color=guide_legend(title="Age-HIV status"),shape=guide_legend(title="On ART/ABX?"),fill=FALSE) +
  theme(legend.text=element_text(size=10),legend.position="right",legend.title=element_text(face="bold",size=10))

#format and present average carriage duration by ABX
phirst.es$clear.est <- 1/phirst.es$clear.est
phirst.es$Lclear.est <- 1/phirst.es$Lclear.est
phirst.es$Uclear.est <- 1/phirst.es$Uclear.est
phirst.es <- subset(phirst.es, select=c(iid,age,hiv,abx,clear.est,Uclear.est,Lclear.est))
phirst.es <- rename(phirst.es, c("iid"="Table.2B","age"="Age","hiv"="HIV","abx"="ABX","clear.est"="Estimate","Uclear.est"="Lower.bound","Lclear.est"="Upper.bound"))

#average carriage duration by ART status
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Younger child",artv="No"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Younger child",artv="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Older child",artv="No"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Older child",artv="No"),ci="normal",cl=0.95)
p.modele <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",artv="No"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",artv="No"),ci="normal",cl=0.95)

phirst.es0 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"),
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "artv"=c("No","No","No","No","No","No"))

phirst.es0$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1],p.modele$estimates[2,1],p.modelf$estimates[2,1])
phirst.es0$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1],p.modele$L[2,1],p.modelf$L[2,1])
phirst.es0$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1],p.modele$U[2,1],p.modelf$U[2,1])  

p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Younger child",artv="Yes"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Younger child",artv="No"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Older child",artv="Yes"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Older child",artv="No"),ci="normal",cl=0.95)
p.modele <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult",artv="Yes"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult",artv="No"),ci="normal",cl=0.95)

phirst.es1 <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                         "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"),
                         "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"),
                         "artv"=c("Yes","Yes","Yes","Yes","Yes","Yes"))

phirst.es1$clear.est <- c(p.modela$estimates[2,1],NA,p.modelc$estimates[2,1],NA,p.modele$estimates[2,1],NA)
phirst.es1$Lclear.est <- c(p.modela$L[2,1],NA,p.modelc$L[2,1],NA,p.modele$L[2,1],NA)
phirst.es1$Uclear.est <- c(p.modela$U[2,1],NA,p.modelc$U[2,1],NA,p.modele$U[2,1],NA)  
phirst.es <- rbind(phirst.es0,phirst.es1)

B<-ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est,shape=artv),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,1/clear.est,color=iid,shape=artv),size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() +
  labs(title="B",x="",y="Average carriage duration by ART (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  theme(legend.position="none")

#format and present average carriage duration by ART
phirst.es$clear.est <- 1/phirst.es$clear.est
phirst.es$Lclear.est <- 1/phirst.es$Lclear.est
phirst.es$Uclear.est <- 1/phirst.es$Uclear.est
phirst.es <- subset(phirst.es, select=c(iid,age,hiv,artv,clear.est,Uclear.est,Lclear.est))
phirst.es <- rename(phirst.es, c("iid"="Table.2A","age"="Age","hiv"="HIV","artv"="ART","clear.est"="Estimate","Uclear.est"="Lower.bound","Lclear.est"="Upper.bound"))

#overall average carriage duration
p.modela <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Younger child"),ci="normal",cl=0.95)
p.modelb <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Younger child"),ci="normal",cl=0.95)
p.modelc <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Older child"),ci="normal",cl=0.95)
p.modeld <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Older child"),ci="normal",cl=0.95)
p.modele <- qmatrix.msm(p.model4, covariates=list(hiv="Pos",agecat="Adult"),ci="normal",cl=0.95)
p.modelf <- qmatrix.msm(p.model4, covariates=list(hiv="Neg",agecat="Adult"),ci="normal",cl=0.95)

phirst.es <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                        "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"),
                        "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"))

phirst.es$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1],p.modele$estimates[2,1],p.modelf$estimates[2,1])
phirst.es$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1],p.modele$L[2,1],p.modelf$L[2,1])
phirst.es$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1],p.modele$U[2,1],p.modelf$U[2,1])  

C<-ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=1/Lclear.est,ymax=1/Uclear.est),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,1/clear.est,color=iid),shape=8,size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="C",x="",y="Overall average carriage duration (days)") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  theme(legend.position="none")

#format and present average carriage duration
phirst.es$clear.est <- 1/phirst.es$clear.est
phirst.es$Lclear.est <- 1/phirst.es$Lclear.est
phirst.es$Uclear.est <- 1/phirst.es$Uclear.est
phirst.es <- subset(phirst.es, select=c(iid,age,hiv,clear.est,Uclear.est,Lclear.est))
phirst.es <- rename(phirst.es, c("iid"="Table.2C","age"="Age","hiv"="HIV","clear.est"="Estimate","Uclear.est"="Lower.bound","Lclear.est"="Upper.bound"))

#probability of clearance
p.modela <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Younger child"),ci="normal",cl=0.95)
p.modelb <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Younger child"),ci="normal",cl=0.95)
p.modelc <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Older child"),ci="normal",cl=0.95)
p.modeld <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Older child"),ci="normal",cl=0.95)
p.modele <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Pos",agecat="Adult"),ci="normal",cl=0.95)
p.modelf <- pmatrix.msm(p.model4,t=1,covariates=list(hiv="Neg",agecat="Adult"),ci="normal",cl=0.95)

phirst.es <- data.frame("iid"=c("Younger child HIV+","Younger child HIV-","Older child HIV+","Older child HIV-","Adult HIV+","Adult HIV-"),
                        "age"=c("Younger child","Younger child","Older child","Older child","Adult","Adult"),
                        "hiv"=c("HIV+","HIV-","HIV+","HIV-","HIV+","HIV-"))

phirst.es$clear.est <- c(p.modela$estimates[2,1],p.modelb$estimates[2,1],p.modelc$estimates[2,1],p.modeld$estimates[2,1],p.modele$estimates[2,1],p.modelf$estimates[2,1])
phirst.es$Lclear.est <- c(p.modela$L[2,1],p.modelb$L[2,1],p.modelc$L[2,1],p.modeld$L[2,1],p.modele$L[2,1],p.modelf$L[2,1])
phirst.es$Uclear.est <- c(p.modela$U[2,1],p.modelb$U[2,1],p.modelc$U[2,1],p.modeld$U[2,1],p.modele$U[2,1],p.modelf$U[2,1])  

D<-ggplot(phirst.es) + 
  geom_errorbar(aes(iid,color=iid,ymin=Lclear.est,ymax=Uclear.est),width=0,size=1,position=position_dodge(width=0.5)) +
  geom_point(aes(iid,clear.est,color=iid),shape=8,size=3,position=position_dodge(width=0.5),stat="identity") +
  theme_bw() + 
  labs(title="D",x="",y="Daily probability of carriage clearance") + 
  theme(axis.text.y=element_text(face="bold",size=10)) + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  theme(legend.position="none")

#format and present probability of carriage clearance
phirst.es <- rename(phirst.es, c("iid"="Table.2D","age"="Age","hiv"="HIV","clear.est"="Estimate","Lclear.est"="Lower.bound","Uclear.est"="Upper.bound"))

print(ggarrange(A,B,C,D,ncol=4,common.legend=TRUE,legend="right"))

