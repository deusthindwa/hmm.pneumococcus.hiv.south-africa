#Written by Deus Thindwa
#Estimating the contribution of HIV-infected adults to household pneumococcal transmission in South Africa, 2016-2018.
#Continuous-time time-homogeneous hidden Markov modelling study, PhD chapter 1.
#20/9/2019 - 10/6/2020 

#run multiple chains to assess convergence of the selected model
j=0.05;k=2.00
for(i in 1:5){
  
sink("~/Rproject/Markov.Model.Resources/data/convergence.txt",append=TRUE)
p.convrg <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=rbind(c(0.0,j), c(k,0.0)), ematrix=rbind(c(1.0,0.0), c(0.15,0.85)),
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), obstrue=obst, est.initprobs=T,
                control=list(maxit=200000,trace=1,REPORT=1))
  sink()
  j=j+0.1;k=k-0.4
}

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


#observed versus predicted prevalence of selected model
p.obsexp <- msm(state~dys,subject=ind_id,data=phirst.fu,
                qmatrix=matrix.Q,
                covariates=list("1-2"=~agecat+hiv+ahivcat+tx+hhsizecat,"2-1"=~agecat+hiv+abx+artv),
                censor=9, censor.states=c(1,2), est.initprobs=T,
                opt.method="bobyqa", control=list(maxfun=100000))

phirst.oe <- tk_tbl(prevalence.msm(p.obsexp,times=seq(0,289,14)),preserve_index=TRUE,rename_index="Time")
phirst.oe <- subset(phirst.oe,select=c(Time, 
                                       Observed.State.1, Observed.State.2, Observed.percentages.State.1, Observed.percentages.State.2,
                                       Expected.Clear,Expected.Carry,Expected.percentages.Clear,Expected.percentages.Carry))
phirst.oe <- rename(phirst.oe, c("Observed.State.1"="obs.clear", "Observed.State.2"="obs.carry","Observed.percentages.State.1"="obs.p.clear","Observed.percentages.State.2"="obs.p.carry",
                                 "Expected.Clear"="exp.clear","Expected.Carry"="exp.carry","Expected.percentages.Clear"="exp.p.clear","Expected.percentages.Carry"="exp.p.carry"))

phirst.oe$lci.clear=phirst.oe$exp.p.clear/100-(1.96*sqrt(phirst.oe$exp.p.clear/100*(1-phirst.oe$exp.p.clear/100)/phirst.oe$exp.clear)) 
phirst.oe$uci.clear=phirst.oe$exp.p.clear/100+(1.96*sqrt(phirst.oe$exp.p.clear/100*(1-phirst.oe$exp.p.clear/100)/phirst.oe$exp.clear))
phirst.oe$lci.carry=phirst.oe$exp.p.carry/100-(1.96*sqrt(phirst.oe$exp.p.carry/100*(1-phirst.oe$exp.p.carry/100)/phirst.oe$exp.carry)) 
phirst.oe$uci.carry=phirst.oe$exp.p.carry/100+(1.96*sqrt(phirst.oe$exp.p.carry/100*(1-phirst.oe$exp.p.carry/100)/phirst.oe$exp.carry))

#plot observed and predicted carriage
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

#plot observed versus predicted carriage
phirst.pcla <- subset(phirst.oe,select=c(obs.clear));phirst.pcla$state1<-"Observed clearance"
phirst.pcla <- rename(phirst.pcla, c("obs.clear" = "value1"))
phirst.pclb <- subset(phirst.oe,select=c(exp.clear));phirst.pclb$state2<-"Predicted clearance"
phirst.pclb <- rename(phirst.pclb, c("exp.clear" = "value2"))
phirst.pcl1 <-cbind(phirst.pcla,phirst.pclb)
colx <- c(phirst.pcl1$state1,phirst.pcl1$state2)

phirst.pclc <- subset(phirst.oe,select=c(obs.carry));phirst.pclc$state1<-"Observed carriage"
phirst.pclc <- rename(phirst.pclc, c("obs.carry" = "value1"))
phirst.pcld <- subset(phirst.oe,select=c(exp.carry));phirst.pcld$state2<-"Predicted carriage"
phirst.pcld <- rename(phirst.pcld, c("exp.carry" = "value2"))
phirst.pcl2 <-cbind(phirst.pclc,phirst.pcld)

ggplot() + 
  geom_point(aes(phirst.pcl1$value1,phirst.pcl1$value2),color=,size=3) + 
  geom_point(aes(phirst.pcl2$value1,phirst.pcl2$value2),size=3) + 
  geom_abline(intercept=0,slope=1,linetype=2) + 
  labs(title="C", x="Observed",y="Predicted") + 
  theme_bw() + 
  ylim(350,1150) + 
  xlim(350,1150) + 
  theme(axis.text.x=element_text(face="bold",size=10), axis.text.y=element_text(face="bold",size=10)) +
  theme(legend.position=c(0.2,0.8)) + 
  guides(color=guide_legend(title=""),shape=guide_legend(title="")) 

A + B + C
