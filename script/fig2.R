#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#show swabbing visit vs carriage prevalence
phirst.d1 <- subset(subset(phirst.fu,select=c(visit_id,state,agecat,hiv)),state !=9)
phirst.d1$visit_id <- as.integer(substr(phirst.d1$visit_id,10,12))
phirst.d2 <- phirst.d1 %>% group_by(visit_id,agecat,hiv) %>% tally(state==2)
phirst.d1 <- phirst.d1 %>% group_by(visit_id,agecat,hiv) %>% tally()
phirst.d1$nPos <- phirst.d2$n; phirst.d1$prev <- phirst.d1$nPos/phirst.d1$n; remove(phirst.d2)
phirst.d1$d1group <- if_else(phirst.d1$agecat=="Child" & phirst.d1$hiv=="Neg","Child HIV-",if_else(phirst.d1$agecat=="Child" & phirst.d1$hiv=="Pos","Child HIV+",
                     if_else(phirst.d1$agecat=="Adult" & phirst.d1$hiv=="Neg","Adult HIV-",if_else(phirst.d1$agecat=="Adult" & phirst.d1$hiv=="Pos","Adult HIV+",NULL))))

A<-ggplot(phirst.d1, aes(visit_id,prev*100,color=d1group)) + 
  geom_line(size=1, na.rm=TRUE) + 
  #geom_smooth(data=subset(phirst.d1,d1group !="Child HIV+"),se=TRUE,alpha=0.6,method=loess) +
  labs(title="A",x="Sampling visit",y="Carriage prevalence (%)") + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position="none")

#show number of positive samples per person vs probability density
phirst.d2 <- arrange(subset(subset(phirst.fu,select=c(visit_id,state,agecat,hiv)),state==2),visit_id)
phirst.d2$ind_id <- substr(phirst.d2$visit_id,1,8)
phirst.d2 <- phirst.d2 %>% group_by(ind_id,agecat,hiv) %>% tally()
phirst.d2$d2group <- if_else(phirst.d2$agecat=="Child" & phirst.d2$hiv=="Neg","Child HIV-",if_else(phirst.d2$agecat=="Child" & phirst.d2$hiv=="Pos","Child HIV+",
                     if_else(phirst.d2$agecat=="Adult" & phirst.d2$hiv=="Neg","Adult HIV-",if_else(phirst.d2$agecat=="Adult" & phirst.d2$hiv=="Pos","Adult HIV+",NULL))))

B<-ggplot(phirst.d2) + 
  geom_density(aes(x=n,color=d2group,fill=d2group), alpha=0.2) +
  labs(title="B",x="Number of positive samples per person",y="Probability density") + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position=c(0.7,0.75),legend.title=element_blank(),legend.key.size=unit(0.4,"cm"),legend.key=element_rect(size=5,fill="white",colour=NA),legend.spacing.x=unit(0.2,"cm"))

#show household size vs frequency of carriage
phirst.d3 <- subset(subset(phirst.fu,select=c(hhsize,state,agecat,hiv)),state !=9)
phirst.d3$hhsize <- if_else(phirst.d3$hhsize<=5,"<6 members",if_else(phirst.d3$hhsize>=6 & phirst.d3$hhsize<=10,"6-10 members","11+ members"))
phirst.d4 <- phirst.d3 %>% group_by(hhsize,agecat,hiv) %>% tally(state==2)
phirst.d3 <- phirst.d3 %>% group_by(hhsize,agecat,hiv) %>% tally()
phirst.d3$nPos <- phirst.d4$n; phirst.d3$prev <- phirst.d3$nPos/phirst.d3$n; remove(phirst.d4)
phirst.d3 <- phirst.d3 %>% mutate(lci=prev-(1.96*sqrt(prev*(1-prev)/n)), uci=prev+(1.96*sqrt(prev*(1-prev)/n)))
phirst.d3$d3group <- if_else(phirst.d3$agecat=="Child" & phirst.d3$hiv=="Neg","Child HIV-",if_else(phirst.d3$agecat=="Child" & phirst.d3$hiv=="Pos","Child HIV+",
                     if_else(phirst.d3$agecat=="Adult" & phirst.d3$hiv=="Neg","Adult HIV-",if_else(phirst.d3$agecat=="Adult" & phirst.d3$hiv=="Pos","Adult HIV+",NULL))))
phirst.d3[nrow(phirst.d3) + 1,] = list("11+ members","Child","Pos",NA,NA,0,0,0,"Child HIV+")

C<-ggplot(phirst.d3,aes(factor(hhsize,levels(factor(hhsize))[c(1,3,2)]),prev*100,fill=d3group)) + 
  geom_bar(stat="identity",color="gray50",position=position_dodge(0.9)) + 
  geom_errorbar(aes(ymin=lci*100, ymax=uci*100), width=0.2, position=position_dodge(0.9)) +
  labs(title="C", x="Household size",y="Carriage prevalence (%)") + 
  theme_bw() +
  ylim(0,100) +
  theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position="none",legend.title=element_blank()) 

#show sampling visit vs carriage density
phirst.d4 <- arrange(subset(subset(phirst.fu,select=c(visit_id,npdensity,agecat,hiv)),!is.na(npdensity)),visit_id)
phirst.d4$visit_id <- as.integer(substr(phirst.d4$visit_id,10,12))
phirst.d4$d4group <- if_else(phirst.d4$agecat=="Child" & phirst.d4$hiv=="Neg","Child HIV-",if_else(phirst.d4$agecat=="Child" & phirst.d4$hiv=="Pos","Child HIV+",
                     if_else(phirst.d4$agecat=="Adult" & phirst.d4$hiv=="Neg","Adult HIV-",if_else(phirst.d4$agecat=="Adult" & phirst.d4$hiv=="Pos","Adult HIV+",NULL))))

D<-ggplot(phirst.d4,aes(d4group,npdensity,color=d4group)) + 
  geom_boxplot(outlier.shape=16, outlier.size=2, notch=TRUE, size=1) +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=3) +
  scale_y_log10(labels=trans_format("log10",math_format(10^.x))) +
  labs(title="D",x="",y="Carriage density (GE/ml)") + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position="none",legend.title=element_blank())

grid.arrange(A,B,C,D,nrow=2)
remove(A,B,C,D,phirst.d1,phirst.d2,phirst.d3,phirst.d4)
