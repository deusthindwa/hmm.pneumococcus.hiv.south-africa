#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#11/3/2020

dev.off()  
m.converge <- read.csv("~/Rproject/Markov.Model.Resources/data/convergence.csv")
m.converge$chaincat <- if_else(m.converge$chain==1,"1 (q12=0.05, q21=2.00) >> (q12=0.04, q21=0.06)",
                               if_else(m.converge$chain==2,"2 (q12=0.15, q21=1.60) >> (q12=0.04, q21=0.06)",
                                       if_else(m.converge$chain==3,"3 (q12=0.25, q21=1.20) >> (q12=0.04, q21=0.06)",
                                               if_else(m.converge$chain==4,"4 (q12=0.35, q21=0.80) >> (q12=0.04, q21=0.06)","5 (q12=0.45, q21=0.40) >> (q12=0.04, q21=0.06)"))))

A <- ggplot(m.converge, aes(iter,Lik2,color=chaincat)) + 
  geom_line(size=0.8) + 
  labs(title="A", x="Iteration",y="-2Log-likelihood",position=position_dodge(width=0)) + 
  xlim(0,1000) +
  scale_y_continuous(breaks=c(60000,80000,100000,120000,140000,160000),labels=scientific) + 
  theme_bw() + 
  theme(legend.position=c(0.55,0.6)) + 
  guides(color=guide_legend(title="Chain (initial intensity) >> (final intensity)")) + 
  theme(legend.key.height=unit(0.8,"line")) + 
  theme(legend.key.width=unit(1,"line")) + 
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10))

cols <- c("Observed vs expected clear"="#0000FF","Observed vs expected carry"="#FF0000")
B <- ggplot(m.OEDS,aes(Time)) + 
  geom_point(aes(Time,obs.p.clear,color="Observed vs expected clear"),size=2.4,shape=5) + 
  geom_line(aes(Time,exp.p.clear,color="Observed vs expected clear"),size=1) + 
  geom_ribbon(aes(ymin=lci.clear*100, ymax=uci.clear*100, color="Observed vs expected clear"), alpha=0.2, size=0.1) +
  geom_point(aes(Time,obs.p.carry,color="Observed vs expected carry"),size=2.4,shape=5) + 
  geom_line(aes(Time,exp.p.carry,color="Observed vs expected carry"),size=1) + 
  geom_ribbon(aes(ymin=lci.carry*100, ymax=uci.carry*100, color="Observed vs expected carry"), alpha=0.2, size=0.1) +
  labs(title="B", x="Days",y="Prevalence (%)") + 
  scale_x_continuous(breaks=c(0,40,80,120,160,200,240,280)) + 
  scale_y_continuous(breaks=c(20,30,40,50,60,70,80)) + 
  theme_bw() + 
  theme(legend.position=c(0.3,0.55)) + 
  guides(color=guide_legend(title="")) +
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10))

grid.arrange(A,B,ncol=2)
