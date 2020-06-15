#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020

#plots from viterbi algorithm
A <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_point(aes(time,observed, color=observed), size=2.4, shape=20) + 
  theme_bw() + 
  labs(title="Person A",x="",y="Observed data") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.text.x=element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(5.5,5.5,-5.5,5.5),"pt"))

B <-ggplot(subset(phirst.vi,subject=="A001-001")) + 
  geom_point(aes(time,fitted,color=fitted), size=2.4, shape=20) + 
  theme_bw() + 
  labs(title="",x="",y="Viterbi states") + 
  theme(axis.text.y=element_text(face="bold",size=10),axis.text.x=element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(-5.5,5.5,-5.5,5.5),"pt"))

C <-ggplot(subset(phirst.vi,subject=="A001-001")) +
  geom_line(aes(time,probhs2),color='gray50',size=0.6) +
  theme_bw() + 
  labs(title="",x="Days",y="Probability") +
  scale_y_continuous(labels=percent) +
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(plot.margin=unit(c(-5,5.5,5.5,5.5),"pt"))

D <-ggplot(subset(phirst.vi,subject=="A234-004")) +
  geom_point(aes(time,observed, color=observed), size=2.4, shape=20) + 
  theme_bw() + 
  labs(title="Person B",x="",y="") + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(5.5,15,-5.5,-5.5),"pt"))

E <-ggplot(subset(phirst.vi,subject=="A234-004")) + 
  geom_point(aes(time,fitted,color=fitted), size=2.4, shape=20) + 
  theme_bw() + 
  labs(title="",x="",y="") + 
  theme(axis.text.x=element_blank(),axis.text.y=element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin=unit(c(-5.5,15,-5.5,-5.5),"pt"))

F <-ggplot(subset(phirst.vi,subject=="A234-004")) +
  geom_line(aes(time,probhs2),color='gray50',size=0.6) +
  theme_bw() + 
  labs(title="",x="Days",y="") +
  scale_y_continuous(labels=percent) +
  theme(axis.text.x=element_text(face="bold",size=10),axis.text.y=element_blank()) +
  theme(plot.margin=unit(c(-5,15,5.5,-5.5),"pt"))

grid.arrange(A,D,B,E,C,F,nrow=3,ncol=2)
remove(A,B,C,D,E)
