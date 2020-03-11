#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

#---------------show swabbing visit vs carriage prevalence
phirst.d1 <- subset(phirst.fu,select=c(visit_id,state,age,hiv))
phirst.d1$visit_id <- as.integer(substr(phirst.d1$visit_id,10,12))
phirst.d1$state <- if_else(phirst.d1$state==2L,1L,0L)
phirst.d1 <- subset(phirst.d1, state==1)
phirst.d1 <- arrange(phirst.d1,visit_id)
phirst.d1 <- phirst.d1 %>% group_by(visit_id,age,hiv) %>% tally()
phirst.d1$d1group <- if_else(phirst.d1$age==0L & phirst.d1$hiv==0L,"HIV- child",
                             if_else(phirst.d1$age==0L & phirst.d1$hiv==1L,"HIV+ child",
                                     if_else(phirst.d1$age==1L & phirst.d1$hiv==0L,"HIV- adult",
                                             if_else(phirst.d1$age==1L & phirst.d1$hiv==1L,"HIV+ adult",NULL))))

phirst.d1 <- ddply(subset(phirst.d1,!is.na(d1group)), "visit_id", mutate, prev=n*100/sum(n))

A<-ggplot(phirst.d1, aes(visit_id, prev, color=d1group)) + 
  geom_point(size=1.5, na.rm=TRUE) + 
  geom_smooth(se=TRUE,alpha=0.8) +
  labs(title="A", x="Sampling visit", y="Carriage prevalence (%)") + 
  theme(legend.position="none") + 
  theme_bw() + 
  ylim(0,80) +
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position="none")

#---------------show number of carriage episodes per person vs carriage prevalence
phirst.d2 <- subset(phirst.fu,select=c(ind_id,state,age,hiv))
phirst.d2$state <- if_else(phirst.d2$state==2L,1L,0L)
phirst.d2 <- subset(phirst.d2, state==1)
phirst.d2 <- arrange(phirst.d2,ind_id)
phirst.d2 <- phirst.d2 %>% group_by(ind_id,age,hiv) %>% tally(state)
phirst.d2$d2group <- if_else(phirst.d2$age==0L & phirst.d2$hiv==0L,"HIV- child",
                             if_else(phirst.d2$age==0L & phirst.d2$hiv==1L,"HIV+ child",
                                     if_else(phirst.d2$age==1L & phirst.d2$hiv==0L,"HIV- adult",
                                             if_else(phirst.d2$age==1L & phirst.d2$hiv==1L,"HIV+ adult",NULL))))

B <- ggplot(subset(phirst.d2, !is.na(d2group))) + 
  geom_density(aes(x=n,color=d2group,fill=d2group), alpha=0.2) +
  labs(title="B",x="Number of positive samples per person",y="Probability density") + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position=c(0.7,0.8),legend.title=element_blank(),legend.key.size=unit(0.8,"lines")) 

#---------------show household size vs frequency of carriage
phirst.d3 <- subset(phirst.fu,select=c(hhsize,state,age,hiv))
phirst.d3$hhsize <- as.integer(phirst.d3$hhsize)
phirst.d3$hhsize <- if_else(phirst.d3$hhsize<=5,"<5 members", if_else(phirst.d3$hhsize>=6 & phirst.d3$hhsize<=10,"6-10 members","11+ members"))
phirst.d3$state <- if_else(phirst.d3$state==2L,1L,0L)
phirst.d3 <- subset(phirst.d3, state==1)
phirst.d3 <- arrange(phirst.d3,hhsize)
phirst.d3 <- phirst.d3 %>% group_by(hhsize,age,hiv) %>% tally(state)
phirst.d3$d3group <- if_else(phirst.d3$age==0L & phirst.d3$hiv==0L,"HIV- child",
                             if_else(phirst.d3$age==0L & phirst.d3$hiv==1L,"HIV+ child",
                                     if_else(phirst.d3$age==1L & phirst.d3$hiv==0L,"HIV- adult",
                                             if_else(phirst.d3$age==1L & phirst.d3$hiv==1L,"HIV+ adult",NULL))))

phirst.d3 <- ddply(subset(phirst.d3,!is.na(d3group)), "hhsize", mutate, prev=n/sum(n),
                   lci=prev-(1.96*sqrt(prev*(1-prev)/n)), uci=prev+(1.96*sqrt(prev*(1-prev)/n)))
phirst.d3$lci <- if_else(phirst.d3$lci<0,0,phirst.d3$lci)

C<-ggplot(subset(phirst.d3, !is.na(hhsize)),aes(factor(hhsize,levels(factor(hhsize))[c(1,3,2)]),prev*100,fill=d3group)) + 
  geom_bar(stat="identity",color="gray50",position=position_dodge(0.9)) + 
  geom_errorbar(aes(ymin=lci*100, ymax=uci*100), width=0.2, position=position_dodge(0.9)) +
  labs(title="C", x="Household size",y="Carriage prevalence (%)") + 
  theme_bw() +
  ylim(0,80) +
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position="none",legend.title=element_blank()) 

#---------------show sampling visit vs carriage density
phirst.d4 <- subset(phirst.fu,select=c(visit_id,npspneload,age,hiv))
phirst.d4$visit_id <- as.integer(substr(phirst.d4$visit_id,10,12))
phirst.d4 <- subset(phirst.d4, !is.na(npspneload))
phirst.d4$npspneload <- as.integer(phirst.d4$npspneload)
phirst.d4 <- arrange(phirst.d4,visit_id)
phirst.d4 <- phirst.d4 %>% group_by(visit_id,age,hiv) %>% summarise_all(mean)
phirst.d4$d4group <- if_else(phirst.d4$age==0L & phirst.d4$hiv==0L,"HIV- child",
                             if_else(phirst.d4$age==0L & phirst.d4$hiv==1L,"HIV+ child",
                                     if_else(phirst.d4$age==1L & phirst.d4$hiv==0L,"HIV- adult",
                                             if_else(phirst.d4$age==1L & phirst.d4$hiv==1L,"HIV+ adult",NULL))))

D<-ggplot(subset(phirst.d4, !is.na(d4group)),aes(x=d4group,y=npspneload,color=d4group)) + 
  geom_boxplot(outlier.shape=16, outlier.size=2, notch=TRUE, size=1) +
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=3) +
  scale_y_log10(labels=trans_format("log10",math_format(10^.x))) +
  labs(title="D",x="",y="Carriage density (GE/ml)") + 
  theme(legend.position=c(50,80)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position="none",legend.title=element_blank())

grid.arrange(A,B,C,D,nrow=2)
remove(A,B,C,D,phirst.d1,phirst.d2,phirst.d3,phirst.d4, phirst.fu.abx)
