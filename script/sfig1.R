#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#model convergence
phirst.cg <- read.csv("~/Rproject/Markov.Model.Resources/data/convergence.csv")
phirst.cg$chaincat <- if_else(phirst.cg$chain==1,"1 (q12=0.05, q21=2.00) >> (q12=0.04, q21=0.06)",
                               if_else(phirst.cg$chain==2,"2 (q12=0.15, q21=1.60) >> (q12=0.04, q21=0.06)",
                                       if_else(phirst.cg$chain==3,"3 (q12=0.25, q21=1.20) >> (q12=0.04, q21=0.06)",
                                               if_else(phirst.cg$chain==4,"4 (q12=0.35, q21=0.80) >> (q12=0.04, q21=0.06)","5 (q12=0.45, q21=0.40) >> (q12=0.04, q21=0.06)"))))

A <- ggplot(phirst.cg, aes(iter,Lik2,color=chaincat)) + 
  geom_line(size=0.8) + 
  labs(title="A", x="Iteration",y="-2Log-likelihood") + 
  xlim(0,1000) +
  scale_y_continuous(breaks=c(100000,125000,150000,175000,200000,225000,250000),labels=scientific) + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position=c(0.55,0.6),legend.key.height=unit(0.8,"line"),legend.key.width=unit(1,"line")) + 
  guides(color=guide_legend(title="Chain (initial intensity) >> (final baseline intensity)")) 


#observed versus predicted carriage
cols <- c("Observed vs predicted clearance"="#0000FF","Observed vs predicted carriage"="#FF0000")
B <- ggplot(phirst.oe,aes(Time)) + 
  geom_point(aes(Time,obs.p.clear,color="Observed vs predicted clearance"),size=2,shape=5) + 
  geom_line(aes(Time,exp.p.clear,color="Observed vs predicted clearance"),size=1) + 
  geom_ribbon(aes(ymin=lci.clear*100,ymax=uci.clear*100, color="Observed vs predicted clearance"),alpha=0.2,size=0.1) +
  geom_point(aes(Time,obs.p.carry,color="Observed vs predicted carriage"),size=2,shape=5) + 
  geom_line(aes(Time,exp.p.carry,color="Observed vs predicted carriage"),size=1) + 
  geom_ribbon(aes(ymin=lci.carry*100,ymax=uci.carry*100, color="Observed vs predicted carriage"),alpha=0.2,size=0.1) +
  labs(title="B", x="Days",y="Prevalence (%)") + 
  scale_x_continuous(breaks=c(0,40,80,120,160,200,240,280)) + 
  scale_y_continuous(breaks=c(5,30,40,50,60,70,80)) + 
  theme_bw() + 
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position=c(0.25,0.55)) + 
  guides(color=guide_legend(title="")) 
  
grid.arrange(A,B,ncol=2)
remove(A,B,i,j,k,matrix.Ec,matrix.Qc,cols)

